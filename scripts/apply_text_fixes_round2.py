# -*- coding: utf-8 -*-
# Apply accepted editor-review text fixes (round 2) to Manuscript.docx.
# Every replacement is asserted to occur exactly once. Run-level edits only;
# no paragraph re-creation, so character formatting is preserved.
import docx
import os

MW = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PATH = MW + r"\Manuscript.docx"
d = docx.Document(PATH)
n_edit = 0

def replace_in_runs(p, old, new, tag):
    global n_edit
    hits = [r for r in p.runs if old in r.text]
    assert len(hits) == 1, (tag, "hits=", len(hits))
    hits[0].text = hits[0].text.replace(old, new)
    n_edit += 1
    print("ok:", tag)

def append_to_para(p, add, tag):
    global n_edit
    assert len(p.runs) >= 1, tag
    p.runs[-1].text += add
    n_edit += 1
    print("ok:", tag)

# ---- E1 abstract (in table cell) ----
abs_p = None
for t in d.tables:
    for row in t.rows:
        for cell in row.cells:
            for p in cell.paragraphs:
                if "under ten seconds" in p.text:
                    abs_p = p
assert abs_p is not None, "abstract not found"
replace_in_runs(abs_p,
    "the framework evaluates all 1.49 billion feature pairs in under ten seconds, recovers",
    "the framework counts all 1.49 billion feature pairs in 0.7 s of tensor-core kernel time "
    "(9.8 s end-to-end including I/O) and tests every emitted candidate with exact statistics "
    "in a 95-s wall-clock screen, recovers",
    "E1 abstract speed claim")

# ---- E2 para 7: ref-15 boundary (insert after superscript 15) ----
p7 = d.paragraphs[7]
replace_in_runs(p7,
    ". Here we turn that substrate into a scientific instrument.",
    ". The present work is distinct from that substrate: everything downstream of the counting "
    "wave - the deviation screen, the conservative multiplicity basis, the three-level confounding "
    "control, the planted ground-truth characterization of statistical power and approximation "
    "loss, and the end-to-end biological application - is new here. We turn that substrate into "
    "a scientific instrument.",
    "E2 ref15 boundary")

# ---- E3 para 12: router -> full-matrix bridge ----
append_to_para(d.paragraphs[12],
    " This controlled instance isolates routing; the full-matrix counting on the real cohort "
    "bitmap follows.",
    "E3 router bridge")

# ---- E4 para 26: dry-run 2.3x claim -> real-bitmap head-to-head ----
replace_in_runs(d.paragraphs[26],
    "On the same synthetic planted data used for baseline development, the pair-level recall of "
    "the framework (61/150) is 2.3 times that of Scoary (27/150) and pyseer FE/LMM (26/150 each). "
    "Figure 3a and Table 3 give the support-resolved recall and the FDR and support-floor "
    "sensitivity.",
    "Re-baselined head-to-head on this same real-cohort instance (both directions per planted pair; "
    "Bonferroni 0.05/54,590 per test; pair called when a direction passes), recall is 84/150 for "
    "the exact screen versus 80/150 for Scoary (v1.6.16, structure-naive Fisher) and 79/150 for "
    "pyseer fixed effects, while pyseer's linear mixed model recovers only 17/150 - consistent "
    "with lineage absorption of exactly the signal that the stratified permutation null preserves "
    "(Fig. 6c). The earlier synthetic dry-run instance used during baseline development had "
    "exaggerated this gap and is superseded by these real-cohort numbers. Two differences beyond "
    "recall remain: the exact Fisher layer is uniformly more powerful than the asymptotic score "
    "test used by pyseer (smaller p-value in all 46 planted pairs where both are finite, median "
    "1.05 orders of magnitude), and the screen has no candidate-generation step at all, so no "
    "interacting partner can be discarded before testing. Figure 3a and Table 3 give the "
    "support-resolved recall and the FDR and support-floor sensitivity.",
    "E4 baseline narrative")

# ---- E5 para 31: cross-validation paragraph rewrite ----
p31 = d.paragraphs[31]
assert p31.text.startswith("Cross-validation against reference microbial GWAS tools."), "p31 head"
old31_tail = "the pair test recovers pure exclusion signals to which single-variant GWAS logistic association is structurally blind."
assert p31.text.endswith(old31_tail), "p31 tail"
assert len(p31.runs) == 1, "p31 single run"
p31.runs[0].text = (
    "Cross-validation against reference microbial GWAS tools. To position the pair test against "
    "established microbial GWAS methodology, we ran pyseer (v1.3.11) fixed-effect and linear "
    "mixed-model analyses, and Scoary (v1.6.16), on the same 150-pair planted benchmark and the "
    "same real-cohort bitmap, using each planted family as the tested variant against its "
    "partner's presence/absence phenotype (Bonferroni threshold 0.05/54,590 = 9.2\u00d710\u207b\u2077 "
    "per direction; Scoary's own Bonferroni correction over 54,590 genes). Recall is "
    "statistically comparable across the structure-naive exact tests (Section above), and the "
    "exact Fisher layer is uniformly more powerful than pyseer's asymptotic score test (smaller "
    "p in all 46 pairs with finite values on both sides, median 1.05 orders of magnitude). The "
    "practical difference is what the baseline protocols cost at cohort scale: both baselines "
    "require one model fit per candidate pair, so an exhaustive C(54,590,2) screen is out of reach "
    "for them - reproducing this 150-pair benchmark alone took 300 Scoary runs (about 14.6 "
    "CPU-hours, after a semantics-preserving parser speed patch) against a single 95-s screen "
    "that tests every pair at once. Conversely, lineage control moves in opposite directions: "
    "pyseer's LMM, run on all 150 planted pairs with the cohort kinship matrix, recovers only "
    "17/150, demonstrating that the component absorbed by kinship models is precisely the "
    "lineage-recurrent signal that the stratified permutation null preserves and flags (Fig. 6c). "
    "The two views are therefore complementary: pyseer and Scoary answer phenotype-first "
    "questions under population-structure control, whereas exact pair mining answers "
    "genotype-first questions without a candidate set, and the pair test recovers pure "
    "exclusion signals to which single-variant GWAS association is structurally blind.")
n_edit += 1
print("ok: E5 cross-validation rewrite")

# ---- E6 para 46: B2 LMM-for-calling rebuttal ----
append_to_para(d.paragraphs[46],
    " A natural question is why the calling layer itself does not use a structure-corrected "
    "mixed model. Exhaustive pair-level mixed-model fitting is computationally prohibitive - "
    "each pair fit carries a cubic-complexity kinship solve - whereas the stratified "
    "permutation null controls confounding at cohort scale with a single GPU recomputation "
    "and is fully auditable. We therefore separate the two problems: lineage-stratified "
    "permutation for cohort-scale calling, and mixed models for targeted confirmation, with "
    "agreement between them (Fig. 8b) as the validation of this split.",
    "E6 B2 LMM rebuttal")

# ---- E7 para 47: A5 CPU statistical phase ----
append_to_para(d.paragraphs[47],
    " The statistical phase, which dominates wall-clock time at cohort scale (67-71%), is "
    "deliberately executed on the CPU in this work: per-pair Fisher tests are embarrassingly "
    "parallel across candidate pairs, and mapping them onto the same wave contract as counting "
    "is a direct extension that we leave to future work.",
    "E7 A5 CPU phase")

# ---- E8 para 56: Fig 3d caption ----
replace_in_runs(d.paragraphs[56],
    "d, Planted-interaction recall at BH-FDR < 0.05: 61/150 for the exact screen versus 27 "
    "(Scoary) and 26 (pyseer FE and LMM); against bitmap ground truth the exact screen recalls "
    "46/150, all misses lying below the detectability floor (Fig. 5).",
    "d, Planted-interaction recall on the real K. pneumoniae bitmap at matched Bonferroni "
    "0.05/54,590 (both directions per planted pair): 84/150 for the exact screen versus 80 "
    "(Scoary 1.6.16), 79 (pyseer FE) and 17 (pyseer LMM); all exact-screen misses lie below "
    "the detectability floor (Fig. 5).",
    "E8 Fig3d caption")

# ---- E9 para 109: A3 memory-bandwidth bound ----
append_to_para(d.paragraphs[109],
    " Aggregate memory traffic is quadratic in M at fixed N (approximately 2 \u00d7 C(M,2) "
    "\u00d7 N/8 bytes, ~0.74 TB for the K. pneumoniae matrix), so counting throughput is "
    "memory-bandwidth-bound: a single core already approaches the streaming limit of a "
    "16-thread host (32.0 s versus 28.94 s), which is why thread count improves counting time "
    "only marginally.",
    "E9 A3 bandwidth")

# ---- E10 para 115: positive control pre-declaration + Scoary protocol ----
append_to_para(d.paragraphs[115],
    " For the positive control, the seven AMR co-carriage combinations were fixed from the "
    "literature before the pipeline was run; none was added or removed on the basis of the "
    "results. Scoary (v1.6.16) was run with --no_pairwise (structure-naive Fisher, matching the "
    "fixed-effect comparison); pair-level calls required both directions to pass Bonferroni "
    "0.05; its gene-matrix parser was locally patched with a semantics-preserving fast path "
    "(identical data structures and statistics, parse speed only), and the patch ships with "
    "the reproduction code.",
    "E10 Methods additions")

cp = d.core_properties
cp.author = "DAPR"
cp.last_modified_by = "DAPR"
try:
    cp.comments = ""
except Exception:
    pass
d.save(PATH)
print("ALL OK, edits:", n_edit)
