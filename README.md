# CRISPR Guard

Live app: https://likhitha-tech-eng-crispr-guard-app-l1cp4s.streamlit.app/

A Streamlit tool that scans a real NCBI sequence record for SpCas9 (NGG PAM)
off-target sites, scores them with the MIT mismatch-position model
(Hsu et al. 2013), and annotates each site using the record's own
gene/CDS/exon features to rank biological risk.

## Features
- **Single guide scan** — enter a 20 nt guide, see every matching site up to
  N mismatches, a risk scatter plot, a mismatch-position heatmap, a CSV
  export, and a one-click shareable PDF report.
- **Rank guides in a gene** — scans every possible guide inside a chosen
  gene/CDS and ranks them by predicted specificity.
- **Chromatin filter** — upload a real ENCODE DNase/ATAC BED file to down-weight
  sites in closed chromatin, where Cas9 cannot physically bind.
- **Batch compare guides** — paste or upload a list of candidate guides and
  score/rank them side by side, for choosing between several designs.
- **Custom sequence support** — paste or upload your own FASTA sequence
  instead of an NCBI accession, for private or unpublished sequences. No
  gene annotation is available in this mode, so every site is labeled
  Intergenic — stated plainly in the app, not hidden.
- **Allele-specific guide design** — given a SNP position and its two
  alleles, finds guides that can selectively cut one allele and spare the
  other. This is a real strategy for targeting a single dominant disease
  allele without touching the healthy copy. Two mechanisms are reported:
  PAM-exclusive (the variant creates/destroys the PAM itself, near-total
  selectivity) and mismatch-based (the variant falls inside the
  protospacer, so selectivity depends on position, per the Hsu et al. 2013
  seed-region model). Verified against known seed-region behaviour with
  constructed test cases: a mismatch next to the PAM gives much higher
  selectivity than one far from it, as expected.
- **Paired-guide rearrangement risk** — when two guides are used together
  (paired nickases, excision strategies, multiplexed editing), flags pairs
  of off-target sites that land close together in the genome. Independent
  per-guide scoring misses this: two simultaneous double-strand breaks near
  each other can cause large deletions or chromosomal rearrangements rather
  than simple indels (Kosicki et al. 2018, *Nature Biotechnology*).

## Method
1. Every 20 nt sequence immediately followed by an NGG PAM is indexed on both
   strands (`index_sites`).
2. A candidate guide is compared against every indexed site; sites within the
   allowed mismatch count are kept.
3. Each site gets an MIT score (0–100) from the position and number of its
   mismatches — mismatches near the PAM are penalized more heavily, matching
   the known "seed region" effect.
4. Each site is classified as coding/non-coding/intergenic using the
   record's real GenBank features, and (optionally) a small hard-coded list
   of cancer-associated genes.
5. `risk = (MIT score / 100) × locus severity × 10`. Locus severity is my own
   heuristic weighting, not an experimentally measured value.

## Validation
A TP53 guide (GATAATCTGGCCGGGTGTGG) was run in this tool and in CRISPOR (human genome, SpCas9, NGG PAM).

- This tool found 2 sites in the loaded NG_017013.2 record: the guide's own perfect match, and one 4-mismatch site (mismatch positions 3, 5, 8, 12) in a non-coding part of a gene.
- CRISPOR reported more than 500 off-target sites genome-wide.

The counts differ because this tool searches only the loaded sequence record (about 33 kb), while CRISPOR searches the whole genome. This check shows the scope difference. It does not validate the scoring model against experimental data.

## Known limitations
- SpCas9 NGG PAM only — no NAG, no bulges.
- MIT score is an older model than CFD or deep-learning-based predictors.
- Severity weights are a heuristic, not measured data.
- Off-target search in "Rank guides" is limited to the loaded record, not a
  full genome.
- The chromatin filter assumes BED coordinates already match the loaded
  record's own numbering (0 = start of record). A chromosome-scale ENCODE
  BED file needs its coordinates offset first — this tool does not do that
  automatically.
- All off-target searches cover only the loaded record, not a full genome.
- The chromatin filter has not been tested with real ENCODE data. Treat it as experimental.
- The cancer-gene list is a small demo list of 12 genes.
- Paired-guide risk flags nearby site pairs by distance. It is a proximity flag, not a probability, and it includes each       guide's own target site.
- Allele-specific selectivity is a model-based estimate and has not been validated experimentally.
- Pasted sequences have no annotation, so all sites show as Intergenic.
- If the NCBI download fails with an SSL certificate error, the loader retries once without certificate verification. This     was added for a restricted network.

## Running locally
```
pip install -r requirements.txt
python -m streamlit run app.py
```

## Example accessions to try
- `NC_001416.1` — lambda phage (fast, ~50 kb, good for a first test)
- `NG_017013.2` — human TP53 RefSeqGene
- `NC_000913.3` — *E. coli* K-12 (larger, indexing takes 20–40s)
  
## Quick example
1. Open the app, set the accession to NG_017013.2 and enter your email.
2. In "Single guide", enter GATAATCTGGCCGGGTGTGG.
3. Expected: 2 sites (one perfect match, one 4-mismatch site).

## References
- Hsu et al. (2013). DNA targeting specificity of RNA-guided Cas9 nucleases. Nature Biotechnology.
- Kosicki et al. (2018). Repair of double-strand breaks induced by CRISPR-Cas9 leads to large deletions and complex rearrangements. Nature Biotechnology.
