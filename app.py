"""
CRISPR Guard - web interface
"""
import random
import re

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import streamlit as st

from pipeline import engine

st.set_page_config(page_title="CRISPR Guard", layout="wide")

COLORS = {"Coding (cancer gene)": "#c0392b", "Coding": "#e67e22",
          "Gene (non-coding part)": "#2980b9", "Intergenic": "#95a5a6"}


@st.cache_resource(show_spinner="Loading and indexing the sequence...")
def load_from_ncbi(acc, email):
    rec = engine.load_record(acc, email)
    seq = str(rec.seq)
    return rec.description, len(seq), seq, engine.index_sites(seq), engine.build_regions(rec)


@st.cache_resource(show_spinner="Indexing your sequence...")
def load_from_custom(fasta_text):
    header, seq = engine.parse_fasta_text(fasta_text)
    return header, len(seq), seq, engine.index_sites(seq), []  # no annotation available


st.title("CRISPR Guard: Cas9 off-target risk scanner")
st.caption("Scans a real sequence for NGG-adjacent sites, scores them "
           "with MIT mismatch weights, and weights them by genomic context.")

source = st.sidebar.radio("Sequence source", ["NCBI accession", "Paste my own sequence"])
max_mm = st.sidebar.slider("Max mismatches", 0, 5, 4)

if source == "NCBI accession":
    acc = st.sidebar.text_input("NCBI accession", "NC_001416.1")
    email = st.sidebar.text_input("Email (only needed for first download of this accession)")
    try:
        desc, length, seq, sites, regions = load_from_ncbi(acc, email)
    except Exception as err:
        st.error(f"Could not load {acc}: {err}")
        st.stop()
else:
    st.sidebar.caption("For a private, unpublished, or non-NCBI sequence. "
                        "No gene annotation is available for pasted sequences, "
                        "so every site will show as Intergenic.")
    fasta_text = st.sidebar.text_area("Paste FASTA (or raw letters)", height=150)
    if not fasta_text.strip():
        st.info("Paste a sequence in the sidebar to begin.")
        st.stop()
    desc, length, seq, sites, regions = load_from_custom(fasta_text)
    if length < 23:
        st.error("Sequence is too short (needs at least 23 bp).")
        st.stop()

st.write(f"**Loaded:** {desc} ({length:,} bp, {len(sites):,} candidate NGG sites)")

tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs(
    ["Single guide", "Rank guides in a gene", "Chromatin filter", "Batch compare guides",
     "Allele-specific design", "Paired-guide risk"]
)

# ---------------- TAB 1: single guide scan ----------------
with tab1:
    if "guide" not in st.session_state:
        st.session_state.guide = ""
    if st.button("Fill in a guide taken from this sequence"):
        st.session_state.guide = random.choice(sites)[2]

    guide = st.text_input("Guide RNA, 20 nt, 5'->3', without PAM", key="guide")
    guide = guide.strip().upper().replace("U", "T")

    if not re.fullmatch(r"[ACGT]{20}", guide):
        st.info("Enter exactly 20 letters from A, C, G, T (U is accepted).")
    else:
        q = engine.guide_quality(guide)
        df = engine.score_hits(engine.find_offtargets(guide, sites, max_mm), regions)
        min_risk = st.slider("Hide hits with risk below", 0, 100, 0)
        if not df.empty:
            df = df[df["risk"] >= min_risk]

        a, b, c, d = st.columns(4)
        a.metric("GC content", f"{q['gc']:.0f}%")
        b.metric("Sites found", len(df))
        c.metric("Perfect matches", int((df["mismatches"] == 0).sum()) if len(df) else 0)
        d.metric("Top risk", f"{df['risk'].max():.1f}" if len(df) else "0")

        if q["gc"] < 40 or q["gc"] > 70:
            st.warning("GC content outside the usual 40-70% window for good guides.")
        if q["poly_t"]:
            st.warning("Contains TTTT, which can terminate Pol III transcription.")

        if df.empty:
            st.success("No sites within the current settings.")
        else:
            if (df["mismatches"] == 0).sum() > 1:
                st.warning("This guide matches perfectly in more than one place. Poor specificity.")

            top = df.iloc[0]
            st.markdown(
                f"**Summary:** {len(df)} sites; highest risk is a {top['mismatches']}-mismatch "
                f"site in *{top['locus']}* ({top['gene']}) at position {top['start']:,} "
                f"on the {top['strand']} strand, risk {top['risk']:.1f}/100."
            )

            left, right = st.columns(2)
            with left:
                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 7))
                sns.stripplot(data=df, x="mismatches", y="risk", hue="locus",
                              palette=COLORS, jitter=0.2, size=7, ax=ax1)
                ax1.set_xlabel("Mismatches")
                ax1.set_ylabel("Risk score")
                ax1.legend(fontsize=7)
                for name, sub in df.groupby("locus"):
                    ax2.scatter(sub["start"], sub["risk"], label=name, color=COLORS[name], s=25)
                ax2.set_xlabel("Position on sequence (bp)")
                ax2.set_ylabel("Risk score")
                plt.tight_layout()
                st.pyplot(fig)

            with right:
                st.dataframe(
                    df[["strand", "start", "mismatches", "mm_positions", "site",
                        "pam", "locus", "gene", "mit_score", "risk"]].round(2),
                    use_container_width=True, height=430,
                )
                st.download_button("Download CSV", df.to_csv(index=False).encode(),
                                    "offtargets.csv", "text/csv")

                report_meta = {
                    "Source sequence": desc,
                    "Sequence length": f"{length:,} bp",
                    "GC content of guide": f"{q['gc']:.0f}%",
                    "Max mismatches allowed": max_mm,
                    "Sites found": len(df),
                }
                pdf_bytes = engine.build_report_pdf(guide, report_meta, df)
                st.download_button("Download PDF report", pdf_bytes,
                                    "crispr_guard_report.pdf", "application/pdf")

            st.subheader("Where do off-targets tolerate mismatches?")
            position_counts = np.zeros(20)
            for positions in df.loc[df["mismatches"] > 0, "mm_positions"]:
                for p in positions.split(","):
                    position_counts[int(p) - 1] += 1
            fig2, axh = plt.subplots(figsize=(8, 1.5))
            sns.heatmap(position_counts.reshape(1, -1), cmap="Reds", cbar=False,
                        xticklabels=range(1, 21), yticklabels=["mismatch\ncount"], ax=axh)
            axh.set_xlabel("Guide position (1 = furthest from PAM, 20 = next to PAM)")
            st.pyplot(fig2)
            st.caption("Positions near the PAM (right side) tolerating fewer mismatches "
                       "matches the known 'seed region' specificity model.")

        with st.expander("Method and limitations"):
            st.markdown(
                "- SpCas9 with NGG PAM only (NAG and bulges ignored).\n"
                "- MIT weights (Hsu et al. 2013), an older model than CFD or deep-learning scores.\n"
                "- Severity weights are my own heuristic, not experimentally measured.\n"
                "- No chromatin data by default; see the Chromatin filter tab."
            )

# ---------------- TAB 2: rank every guide in a gene ----------------
with tab2:
    st.subheader("Find and rank all guides inside a gene/CDS")
    gene_names = sorted(set(r[3] for r in regions if r[2] in ("gene", "CDS")))
    if not gene_names:
        st.info("This record has no gene/CDS features to rank guides within.")
    else:
        chosen_gene = st.selectbox("Gene", gene_names)
        matches = [r for r in regions if r[3] == chosen_gene and r[2] == "CDS"]
        if not matches:
            matches = [r for r in regions if r[3] == chosen_gene]
        g_start, g_end = min(m[0] for m in matches), max(m[1] for m in matches)
        st.write(f"Scanning **{chosen_gene}**: {g_start:,}-{g_end:,} ({g_end - g_start} bp)")

        if st.button("Rank guides in this gene"):
            with st.spinner("Scanning every candidate guide against the whole sequence..."):
                ranked = engine.rank_guides(seq, sites, g_start, g_end, regions, max_mm)
            if ranked.empty:
                st.warning("No NGG sites found in this region.")
            else:
                st.dataframe(ranked.head(30).round(2), use_container_width=True)
                st.caption("Higher specificity_score = fewer / weaker off-target sites "
                           "elsewhere in this sequence. Off-targets are only checked "
                           "against the loaded record, not a full genome.")

# ---------------- TAB 3: chromatin accessibility filter ----------------
with tab3:
    st.subheader("Filter by real accessibility data (ENCODE DNase/ATAC BED)")
    st.caption("Upload a BED/narrowPeak file whose coordinates match the loaded "
               "accession's own numbering (0 = start of this record, not the "
               "chromosome). Run a single-guide scan in Tab 1 first.")
    bed_file = st.file_uploader("BED file", type=["bed", "narrowPeak", "txt"])
    if bed_file:
        if "guide" not in st.session_state or not re.fullmatch(r"[ACGT]{20}", st.session_state.get("guide", "").upper().replace("U", "T")):
            st.warning("Enter a valid guide in Tab 1 first, then come back here.")
        else:
            g = st.session_state.guide.strip().upper().replace("U", "T")
            base_df = engine.score_hits(engine.find_offtargets(g, sites, max_mm), regions)
            if base_df.empty:
                st.info("No sites to filter for this guide.")
            else:
                peaks = engine.load_bed_peaks(bed_file)
                st.write(f"Loaded {len(peaks)} accessible regions.")
                adjusted = engine.apply_chromatin(base_df, peaks)
                st.dataframe(adjusted.round(2), use_container_width=True)
                st.caption(f"{int(adjusted['open_chromatin'].sum())} of {len(adjusted)} sites "
                           "fall in open chromatin and keep their full risk score; the rest "
                           "are scaled down 10x since closed chromatin blocks Cas9 binding.")

# ---------------- TAB 4: batch-compare multiple candidate guides ----------------
with tab4:
    st.subheader("Compare several candidate guides at once")
    st.caption("Paste one guide per line (20 nt each), or upload a CSV/text file "
               "with one guide per line. Useful when you've designed several "
               "candidates and need to pick the best one.")

    pasted = st.text_area("Guides, one per line", height=150,
                           placeholder="ATCGGTCGATCGATCGATCG\nGGGCATTAGCATTAGCATTA\n...")
    uploaded = st.file_uploader("...or upload a file (one guide per line)",
                                 type=["csv", "txt"])

    guide_list = []
    if uploaded:
        guide_list = [line.decode().strip() for line in uploaded if line.decode().strip()]
    elif pasted.strip():
        guide_list = [line.strip() for line in pasted.strip().splitlines() if line.strip()]

    if guide_list and st.button(f"Score {len(guide_list)} guides"):
        with st.spinner("Scoring every guide..."):
            batch_df = engine.score_guides_batch(guide_list, sites, regions, max_mm)
        if batch_df.empty:
            st.warning("No guides to score.")
        else:
            n_valid = int(batch_df["valid"].sum())
            n_invalid = len(batch_df) - n_valid
            st.write(f"**{n_valid} scored**, {n_invalid} skipped (not exactly 20 letters of A/C/G/T).")
            st.dataframe(batch_df.round(2), use_container_width=True)
            st.download_button("Download comparison CSV", batch_df.to_csv(index=False).encode(),
                                "batch_guide_comparison.csv", "text/csv")
            st.caption("rank_score combines off-target count/risk with GC content and "
                       "poly-T penalties -- higher is better. This is a design-time "
                       "heuristic for comparing candidates, not a validated efficiency score.")

# ---------------- TAB 5: allele-specific guide design around a SNP ----------------
with tab5:
    st.subheader("Design a guide that cuts one allele and spares the other")
    st.caption(
        "A real strategy in gene therapy: target a dominant disease allele "
        "(e.g. a point mutation) without cutting the healthy copy. Give a "
        "position and the two alleles; this finds guides overlapping that "
        "position and classifies how they'd behave on each allele."
    )
    snp_pos = st.number_input("SNP position (0-based, within the loaded sequence)",
                               min_value=0, max_value=max(length - 1, 0), value=min(1000, length - 1))
    c1, c2 = st.columns(2)
    ref_base = c1.selectbox("Reference base", ["A", "C", "G", "T"], index=2)
    alt_base = c2.selectbox("Alternate (variant) base", ["A", "C", "G", "T"], index=0)
    actual_base = seq[int(snp_pos)].upper()
    if actual_base != ref_base:
        st.warning(f"Note: the loaded sequence actually has '{actual_base}' at this position, "
                   f"not '{ref_base}'. The calculation below still runs, but double check "
                   f"your position if that's unexpected.")

    if ref_base == alt_base:
        st.info("Reference and alternate base are the same -- pick two different bases.")
    elif st.button("Find allele-specific guides"):
        allele_df = engine.design_allele_specific_guides(seq, int(snp_pos), ref_base, alt_base)
        if allele_df.empty:
            st.warning("No candidate guides overlap this position within the search window.")
        else:
            n_exclusive = (allele_df["mechanism"].str.contains("PAM-exclusive")).sum()
            st.write(f"**{len(allele_df)} candidate guide(s) found** "
                     f"({n_exclusive} are PAM-exclusive -- the strongest form of selectivity).")
            st.dataframe(allele_df, use_container_width=True)
            st.caption(
                "selectivity = 100 means the guide can essentially only act on one allele. "
                "PAM-exclusive guides get 100 because the PAM itself only exists on one "
                "allele. Mismatch-based guides get a partial score from the MIT "
                "seed-region model -- mismatches near the PAM matter far more than ones "
                "far from it, matching Hsu et al. 2013."
            )

# ---------------- TAB 6: paired-guide rearrangement risk ----------------
with tab6:
    st.subheader("Risk from using two guides together")
    st.caption(
        "Scoring guides independently misses a real risk: when two Cas9 cuts "
        "happen near each other at the same time, cells can misjoin the ends, "
        "producing large deletions or rearrangements instead of small indels "
        "at each site (Kosicki et al. 2018, Nature Biotechnology). Relevant "
        "whenever two guides are used together -- paired nickases, excision "
        "strategies, or multiplexed editing."
    )
    gc1, gc2 = st.columns(2)
    guide_a = gc1.text_input("Guide A (20 nt)", key="guide_a_pair").strip().upper().replace("U", "T")
    guide_b = gc2.text_input("Guide B (20 nt)", key="guide_b_pair").strip().upper().replace("U", "T")
    dist_threshold = st.slider("Flag pairs within this distance (bp)", 100, 100_000, 10_000, step=100)

    valid_a = bool(re.fullmatch(r"[ACGT]{20}", guide_a))
    valid_b = bool(re.fullmatch(r"[ACGT]{20}", guide_b))
    if guide_a and not valid_a:
        st.info("Guide A: enter exactly 20 letters from A, C, G, T.")
    if guide_b and not valid_b:
        st.info("Guide B: enter exactly 20 letters from A, C, G, T.")

    if valid_a and valid_b:
        if st.button("Check paired rearrangement risk"):
            df_a = engine.score_hits(engine.find_offtargets(guide_a, sites, max_mm), regions)
            df_b = engine.score_hits(engine.find_offtargets(guide_b, sites, max_mm), regions)
            pair_risk = engine.paired_guide_rearrangement_risk(df_a, df_b, dist_threshold)
            if pair_risk.empty:
                st.success("No off-target site pairs found within this distance threshold.")
            else:
                st.error(f"{len(pair_risk)} close pair(s) found -- potential rearrangement risk.")
                st.dataframe(pair_risk, use_container_width=True)
                st.caption("distance_bp is how far apart the two off-target cut sites are. "
                           "combined_risk averages each site's individual risk score.")
