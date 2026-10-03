"""
CRISPR Guard - backend engine
Scans a real NCBI sequence record for SpCas9 (NGG PAM) candidate sites,
scores them with the MIT mismatch-position model (Hsu et al. 2013), and
annotates them using the record's own gene/CDS/exon features.
"""
import os
import re
import ssl
import urllib.error
import urllib.request
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Seq import Seq

# Position weights from Hsu et al. 2013 (MIT score).
# Index 0 = PAM-distal end of the guide, index 19 = next to the PAM.
MIT_W = [0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079,
         0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583]

# Small demo list; only meaningful for human records.
CANCER_GENES = {"TP53", "PTEN", "RB1", "BRCA1", "BRCA2", "APC", "NF1", "VHL",
                 "MYC", "KRAS", "EGFR", "BRAF"}

# How bad is a cut in each kind of region (my own heuristic, not measured).
SEVERITY = {"Coding (cancer gene)": 10, "Coding": 6,
            "Gene (non-coding part)": 3, "Intergenic": 1}

# 20 nt protospacer + any base + GG, found on one strand.
PAM_RE = re.compile(r"(?=([ACGT]{20}[ACGT]GG))")


def load_record(accession, email="", allow_insecure_fallback=True):
    """Read a GenBank record, downloading it once and caching it in data/.

    Some school/college networks intercept HTTPS traffic with their own
    certificate, which makes the normal secure download fail with an SSL
    "self-signed certificate" error even though NCBI itself is fine. If
    that happens, and allow_insecure_fallback is True, we retry once
    without certificate verification so students on those networks aren't
    blocked. This weakens the security of that one download only -- it
    does not affect anything else on the computer.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(here, "..", "data")
    os.makedirs(data_dir, exist_ok=True)
    path = os.path.join(data_dir, accession + ".gb")
    if not os.path.exists(path):
        if not email:
            raise ValueError("Enter an email so NCBI Entrez allows the download.")
        Entrez.email = email
        try:
            with Entrez.efetch(db="nuccore", id=accession,
                                rettype="gbwithparts", retmode="text") as handle:
                text = handle.read()
        except urllib.error.URLError as err:
            is_ssl_issue = isinstance(err.reason, ssl.SSLError) or "CERTIFICATE" in str(err)
            if not (allow_insecure_fallback and is_ssl_issue):
                raise
            # Retry once with certificate verification disabled.
            insecure_ctx = ssl.create_default_context()
            insecure_ctx.check_hostname = False
            insecure_ctx.verify_mode = ssl.CERT_NONE
            old_opener = urllib.request._opener
            urllib.request.install_opener(
                urllib.request.build_opener(urllib.request.HTTPSHandler(context=insecure_ctx))
            )
            try:
                with Entrez.efetch(db="nuccore", id=accession,
                                    rettype="gbwithparts", retmode="text") as handle:
                    text = handle.read()
            finally:
                urllib.request._opener = old_opener  # restore normal secure behaviour
        with open(path, "w") as f:
            f.write(text)
    return SeqIO.read(path, "genbank")


def index_sites(seq):
    """Find every NGG-adjacent 20-mer on both strands, once.
    Returns list of (strand, start_on_forward_coords, protospacer, pam)."""
    seq = seq.upper()
    total = len(seq)
    rev = str(Seq(seq).reverse_complement())
    sites = []
    for strand, s in (("+", seq), ("-", rev)):
        for m in PAM_RE.finditer(s):
            p = m.start()
            hit = m.group(1)
            start = p if strand == "+" else total - p - 23
            sites.append((strand, start, hit[:20], hit[20:]))
    return sites


def build_regions(record):
    """Turn the record's gene/CDS/exon features into (start, end, kind, name)."""
    regions = []
    for f in record.features:
        if f.type not in ("gene", "CDS", "exon"):
            continue
        quals = f.qualifiers
        name = quals.get("gene", quals.get("locus_tag", ["unnamed"]))[0]
        for part in f.location.parts:
            regions.append((int(part.start), int(part.end), f.type, name))
    return regions


def classify(start, end, regions):
    label, gene = "Intergenic", "-"
    for s, e, kind, name in regions:
        if s < end and start < e:
            if kind in ("CDS", "exon"):
                if name in CANCER_GENES:
                    return "Coding (cancer gene)", name
                return "Coding", name
            label, gene = "Gene (non-coding part)", name
    return label, gene


def mit_score(mm):
    """Likelihood-style score 0-100 from mismatch positions (0-based)."""
    if not mm:
        return 100.0
    score = 1.0
    for i in mm:
        score *= 1 - MIT_W[i]
    n = len(mm)
    if n > 1:
        mean_gap = (mm[-1] - mm[0]) / (n - 1)
        score *= 1 / (((19 - mean_gap) / 19) * 4 + 1)
    score /= n ** 2
    return score * 100


def find_offtargets(guide, sites, max_mm=4):
    hits = []
    for strand, start, proto, pam in sites:
        mm = []
        for i in range(20):
            if guide[i] != proto[i]:
                mm.append(i)
                if len(mm) > max_mm:
                    break
        if len(mm) > max_mm:
            continue
        hits.append({"strand": strand, "start": start, "end": start + 23,
                      "mismatches": len(mm),
                      "mm_positions": ",".join(str(i + 1) for i in mm),
                      "site": proto, "pam": pam, "mit_score": mit_score(mm)})
    return hits


def score_hits(hits, regions):
    df = pd.DataFrame(hits)
    if df.empty:
        return df
    labels = [classify(r.start, r.end, regions) for r in df.itertuples()]
    df["locus"] = [l[0] for l in labels]
    df["gene"] = [l[1] for l in labels]
    df["risk"] = (df["mit_score"] / 100) * df["locus"].map(SEVERITY) * 10
    return df.sort_values("risk", ascending=False).reset_index(drop=True)


def guide_quality(guide):
    gc = 100 * (guide.count("G") + guide.count("C")) / len(guide)
    return {"gc": gc, "poly_t": "TTTT" in guide}


# ---------- Feature 1: rank every guide inside a gene ----------

def all_guides_in_region(seq, start, end):
    """Every valid 20nt+NGG guide whose protospacer falls inside [start, end)."""
    window = seq[max(0, start - 23):end + 23]
    offset = max(0, start - 23)
    out = []
    for strand, s, proto, pam in index_sites(window):
        real_start = s + offset
        if start <= real_start < end:
            out.append((strand, real_start, proto, pam))
    return out


def rank_guides(seq, sites, region_start, region_end, regions, max_mm=4):
    """For each guide inside a region, score its OWN off-targets against the
    whole indexed sequence, then rank best (fewest / lowest-risk) first."""
    candidates = all_guides_in_region(seq, region_start, region_end)
    rows = []
    for strand, start, proto, pam in candidates:
        hits = find_offtargets(proto, sites, max_mm)
        scored = score_hits(hits, regions)
        off = scored[scored["mismatches"] > 0] if not scored.empty else scored
        rows.append({
            "guide": proto, "strand": strand, "start": start,
            "gc": guide_quality(proto)["gc"],
            "poly_t": guide_quality(proto)["poly_t"],
            "num_offtargets": len(off),
            "worst_offtarget_risk": off["risk"].max() if len(off) else 0.0,
        })
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["specificity_score"] = 100 - df["worst_offtarget_risk"] - df["num_offtargets"] * 2
    df.loc[(df["gc"] < 40) | (df["gc"] > 70), "specificity_score"] -= 15
    df.loc[df["poly_t"], "specificity_score"] -= 20
    return df.sort_values("specificity_score", ascending=False).reset_index(drop=True)


# ---------- Feature 2: chromatin accessibility filter (BED upload) ----------

def load_bed_peaks(file_obj):
    """Parse a BED file (chrom, start, end, ...) into a sorted interval list.
    Assumes coordinates match the loaded record (single-contig records)."""
    peaks = []
    for line in file_obj:
        line = line.decode() if isinstance(line, bytes) else line
        if line.startswith(("track", "#", "browser")):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue
        peaks.append((int(parts[1]), int(parts[2])))
    return sorted(peaks)


def is_in_open_chromatin(pos, peaks):
    for s, e in peaks:
        if s <= pos < e:
            return True
        if s > pos:
            break
    return False


def apply_chromatin(df, peaks):
    df = df.copy()
    df["open_chromatin"] = df["start"].apply(lambda p: is_in_open_chromatin(p, peaks))
    df.loc[~df["open_chromatin"], "risk"] *= 0.1
    return df.sort_values("risk", ascending=False).reset_index(drop=True)


# ---------- Feature 3: batch scoring of many candidate guides at once ----------

def score_guides_batch(guides, sites, regions, max_mm=4):
    """Score a list of guide strings together, for comparing candidates.
    Invalid entries (wrong length/letters) are kept in the output with
    valid=False rather than silently dropped, so the user can see why."""
    rows = []
    for raw in guides:
        g = raw.strip().upper().replace("U", "T")
        if not re.fullmatch(r"[ACGT]{20}", g):
            rows.append({"guide": raw, "valid": False, "gc": None, "poly_t": None,
                         "perfect_matches": None, "num_offtargets": None,
                         "worst_offtarget_risk": None})
            continue
        hits = find_offtargets(g, sites, max_mm)
        scored = score_hits(hits, regions)
        off = scored[scored["mismatches"] > 0] if not scored.empty else scored
        q = guide_quality(g)
        rows.append({
            "guide": g, "valid": True, "gc": round(q["gc"], 1), "poly_t": q["poly_t"],
            "perfect_matches": int((scored["mismatches"] == 0).sum()) if not scored.empty else 0,
            "num_offtargets": len(off),
            "worst_offtarget_risk": round(off["risk"].max(), 1) if len(off) else 0.0,
        })
    df = pd.DataFrame(rows)
    if df.empty or not df["valid"].any():
        return df
    valid = df[df["valid"]].copy()
    valid["rank_score"] = 100 - valid["worst_offtarget_risk"].fillna(0) - valid["num_offtargets"].fillna(0) * 2
    valid.loc[(valid["gc"] < 40) | (valid["gc"] > 70), "rank_score"] -= 15
    valid.loc[valid["poly_t"] == True, "rank_score"] -= 20  # noqa: E712
    df = pd.concat([valid, df[~df["valid"]]], ignore_index=True)
    return df.sort_values("rank_score", ascending=False, na_position="last").reset_index(drop=True)


# ---------- Feature 4: accept a user's own sequence, not just an NCBI accession ----------

def parse_fasta_text(text):
    """Parse pasted FASTA (or plain raw letters) into (header, sequence).
    No gene/CDS annotation is available for a custom sequence, so downstream
    classification will label everything Intergenic -- that limitation is
    surfaced in the UI, not hidden."""
    lines = text.strip().splitlines()
    header, seq_lines = "", []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            header = line[1:].strip()
        else:
            seq_lines.append(line)
    seq = re.sub(r"[^ACGTacgt]", "", "".join(seq_lines)).upper()
    return (header or "Custom pasted sequence"), seq


# ---------- Feature 5: shareable PDF report for a single guide ----------

def build_report_pdf(guide, meta, df):
    """Build a multi-page PDF report (summary, table, chart) and return raw
    bytes, suitable for st.download_button. Uses only matplotlib, which is
    already a required dependency, so no new install is needed."""
    import io
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    buf = io.BytesIO()
    with PdfPages(buf) as pdf:
        # Page 1: title + summary
        fig = plt.figure(figsize=(8.5, 11))
        fig.text(0.08, 0.93, "CRISPR Guard -- Off-Target Report", fontsize=18, weight="bold")
        fig.text(0.08, 0.89, f"Guide RNA: {guide}", fontsize=12, family="monospace")
        y = 0.84
        for k, v in meta.items():
            fig.text(0.08, y, f"{k}: {v}", fontsize=10)
            y -= 0.035
        fig.text(0.08, y - 0.02,
                  "Method: SpCas9 NGG-PAM scan; MIT mismatch-position scoring\n"
                  "(Hsu et al. 2013); locus severity is an original heuristic,\n"
                  "not a measured value. See README for full method and limitations.",
                  fontsize=8.5)
        pdf.savefig(fig)
        plt.close(fig)

        if not df.empty:
            # Page 2: results table (first 25 rows)
            cols = [c for c in ["strand", "start", "mismatches", "mm_positions",
                                 "locus", "gene", "mit_score", "risk"] if c in df.columns]
            shown = df[cols].head(25).round(2)
            fig2, ax2 = plt.subplots(figsize=(8.5, 11))
            ax2.axis("off")
            fig2.text(0.08, 0.95, "Off-target site table", fontsize=13, weight="bold")
            tbl = ax2.table(cellText=shown.values, colLabels=shown.columns,
                             cellLoc="center", bbox=[0.0, 0.55, 1.0, 0.35])
            tbl.auto_set_font_size(False)
            tbl.set_fontsize(6.5)
            pdf.savefig(fig2)
            plt.close(fig2)

            # Page 3: risk chart
            fig3, ax3 = plt.subplots(figsize=(8.5, 6))
            for name, sub in df.groupby("locus"):
                ax3.scatter(sub["mismatches"], sub["risk"], label=name, s=40)
            ax3.set_xlabel("Mismatches")
            ax3.set_ylabel("Risk score")
            ax3.set_title("Off-target risk by mismatch count")
            ax3.legend(fontsize=8)
            pdf.savefig(fig3)
            plt.close(fig3)
    buf.seek(0)
    return buf.getvalue()


# ---------- Feature 6: allele-specific guide design around a SNP ----------

def design_allele_specific_guides(seq, snp_pos, ref_base, alt_base, window=35):
    """Find guides that can selectively cut one allele at a SNP and spare
    the other -- the real strategy used to target a single disease allele
    (e.g. a dominant pathogenic variant) without touching the healthy copy.

    Two distinct mechanisms are reported, because they behave very
    differently in practice:
      - "PAM-exclusive": the SNP creates or destroys the NGG PAM itself, so
        the guide can physically only act on one allele (near-complete
        selectivity, the gold standard for allele-specific editing).
      - "Mismatch-based": the SNP falls inside the 20nt protospacer rather
        than the PAM, so Cas9 can still bind and cut both alleles, but the
        mismatched allele is cut less efficiently. Selectivity here is
        partial and position-dependent (Hsu et al. 2013 seed-region
        weighting), not absolute.

    snp_pos is 0-based, in the same coordinate system as the loaded
    sequence. ref_base/alt_base are single letters (A/C/G/T).
    """
    region_start = max(0, snp_pos - window)
    region_end = min(len(seq), snp_pos + window)
    ref_window = seq[region_start:region_end].upper()
    rel_pos = snp_pos - region_start
    if not (0 <= rel_pos < len(ref_window)):
        raise ValueError("SNP position falls outside the loaded sequence.")
    alt_window = ref_window[:rel_pos] + alt_base.upper() + ref_window[rel_pos + 1:]

    ref_sites = {(s, st): (p, pam) for s, st, p, pam in index_sites(ref_window)}
    alt_sites = {(s, st): (p, pam) for s, st, p, pam in index_sites(alt_window)}

    rows = []
    for key in set(ref_sites) | set(alt_sites):
        strand, start = key
        in_ref, in_alt = key in ref_sites, key in alt_sites
        if in_ref and not in_alt:
            proto, _ = ref_sites[key]
            rows.append({"strand": strand, "start": start + region_start, "guide": proto,
                         "mechanism": "PAM-exclusive (cuts REF only)", "selectivity": 100.0})
        elif in_alt and not in_ref:
            proto, _ = alt_sites[key]
            rows.append({"strand": strand, "start": start + region_start, "guide": proto,
                         "mechanism": "PAM-exclusive (cuts ALT only)", "selectivity": 100.0})
        else:
            ref_proto, _ = ref_sites[key]
            alt_proto, _ = alt_sites[key]
            if ref_proto == alt_proto:
                continue  # SNP isn't inside this guide's protospacer -- not allele-specific
            mm = [i for i in range(20) if ref_proto[i] != alt_proto[i]]
            cross_cut_score = mit_score(mm)  # how well this guide still cuts the OTHER allele
            selectivity = round(100 - cross_cut_score, 1)
            rows.append({"strand": strand, "start": start + region_start, "guide": ref_proto,
                         "mechanism": "Mismatch-based (favors REF)", "selectivity": selectivity})
            rows.append({"strand": strand, "start": start + region_start, "guide": alt_proto,
                         "mechanism": "Mismatch-based (favors ALT)", "selectivity": selectivity})
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return df.sort_values("selectivity", ascending=False).reset_index(drop=True)


# ---------- Feature 7: paired-guide rearrangement risk ----------

def paired_guide_rearrangement_risk(df_a, df_b, max_distance=10_000):
    """Flag pairs of off-target sites, one from each of two guides used
    together, that fall close to each other in the genome.

    Scoring guides independently misses a real risk: when two double-strand
    breaks occur near each other at the same time, cells can repair them by
    joining the wrong ends together, producing large deletions or
    chromosomal rearrangements rather than small indels at each site
    independently (Kosicki et al. 2018, Nature Biotechnology). This matters
    whenever two guides are used in the same experiment -- paired nickases,
    excision strategies, or simple multiplexed editing.
    """
    if df_a.empty or df_b.empty:
        return pd.DataFrame()
    rows = []
    for _, a in df_a.iterrows():
        for _, b in df_b.iterrows():
            dist = abs(int(a["start"]) - int(b["start"]))
            if dist <= max_distance:
                rows.append({
                    "guide_A_site": int(a["start"]), "guide_A_locus": a["locus"],
                    "guide_B_site": int(b["start"]), "guide_B_locus": b["locus"],
                    "distance_bp": dist,
                    "combined_risk": round((a["risk"] + b["risk"]) / 2, 1),
                })
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values("distance_bp").reset_index(drop=True)
