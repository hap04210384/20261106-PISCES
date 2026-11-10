# -*- coding: utf-8 -*-
# B1: add references 21-29 (all bibliographic details verified by web search
# on 2026-11-02) with superscript in-text citations at first-appearance order.
import copy
import docx
from docx.text.run import Run
from docx.oxml.ns import qn
import os

MW = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PATH = MW + r"\Manuscript.docx"
d = docx.Document(PATH)

def add_sup_after(p, anchor, num, tag):
    runs = p.runs
    full = "".join(r.text for r in runs)
    pos = full.find(anchor)
    assert pos >= 0, (tag, "anchor not found")
    end = pos + len(anchor)
    acc = 0
    for r in runs:
        ln = len(r.text)
        if acc + ln >= end:
            off = end - acc
            tail = r.text[off:]
            r.text = r.text[:off]
            # superscript citation run
            sup_el = copy.deepcopy(r._element)
            r._element.addnext(sup_el)
            sup_run = Run(sup_el, p)
            for t in sup_el.findall(qn("w:t")):
                sup_el.remove(t)
            t = sup_el.makeelement(qn("w:t"), {})
            t.text = num
            sup_el.append(t)
            rpr = sup_el.find(qn("w:rPr"))
            if rpr is None:
                rpr = sup_el.makeelement(qn("w:rPr"), {})
                sup_el.insert(0, rpr)
            va = rpr.find(qn("w:vertAlign"))
            if va is None:
                va = rpr.makeelement(qn("w:vertAlign"), {})
                rpr.append(va)
            va.set(qn("w:val"), "superscript")
            # tail run
            if tail:
                tail_el = copy.deepcopy(r._element)
                sup_el.addnext(tail_el)
                tail_run = Run(tail_el, p)
                tail_run.text = tail
            print("ok cite:", tag)
            return
        acc += ln
    raise AssertionError(tag + " run walk failed")

# para 31 (single run after rewrite)
p31 = d.paragraphs[31]
add_sup_after(p31, "established microbial GWAS methodology", "21,22", "Voichek+Sheppard")
# refresh runs after split
add_sup_after(d.paragraphs[31], "with the cohort kinship matrix", "23", "Lippert")
assert d.paragraphs[31].runs[-1].text.endswith("structurally blind.")
d.paragraphs[31].runs[-1].text += (
    " Prior GPU epistasis implementations exist24, but they target regression "
    "likelihoods over human SNP pairs rather than exact discrete counting over "
    "genome-by-feature bitmaps.")

# para 46: PopPUNK
add_sup_after(d.paragraphs[46], "complementary to, not competing with, microbial GWAS",
              "25", "PopPUNK")

# para 107: CheckM2, Panaroo
add_sup_after(d.paragraphs[107], "passing CheckM2 quality filters", "26", "CheckM2")
add_sup_after(d.paragraphs[107], "PanGBank Panaroo pipeline", "27", "Panaroo")

# para 113: BH
add_sup_after(d.paragraphs[113], "controlled by Benjamini-Hochberg", "28", "BH")

# para 115: LSH
add_sup_after(d.paragraphs[115], "MinHash-LSH used 128 hash functions", "29", "Gionis")

# ---- reference list entries (order = first appearance) ----
refs = [
 "21. Voichek, Y. & Weigel, D. Identifying genetic variants underlying phenotypic "
 "variation in plants without complete genomes. Nat. Genet. 52, 534\u2013540 (2020).",
 "22. Sheppard, S. K. et al. Genome-wide association study identifies vitamin B5 "
 "biosynthesis as a host specificity factor in Campylobacter. Proc. Natl Acad. Sci. "
 "USA 110, 11923\u201311927 (2013).",
 "23. Lippert, C. et al. FaST linear mixed models for genome-wide association studies. "
 "Nat. Methods 8, 833\u2013835 (2011).",
 "24. Hemani, G., Theocharidis, A., Wei, W. & Haley, C. EpiGPU: exhaustive pairwise "
 "epistasis scans parallelized on consumer level graphics cards. Bioinformatics 27, "
 "1462\u20131465 (2011).",
 "25. Lees, J. A. et al. Fast and flexible bacterial genomic epidemiology with PopPUNK. "
 "Genome Res. 29, 304\u2013316 (2019).",
 "26. Chklovski, A., Parks, D. H., Woodcroft, B. J. & Tyson, G. W. CheckM2: a rapid, "
 "scalable and accurate tool for assessing microbial genome quality using machine "
 "learning. Nat. Methods 20, 1203\u20131212 (2023).",
 "27. Tonkin-Hill, G. et al. Producing polished prokaryotic pangenomes with the Panaroo "
 "pipeline. Genome Biol. 21, 180 (2020).",
 "28. Benjamini, Y. & Hochberg, Y. Controlling the false discovery rate: a practical and "
 "powerful approach to multiple testing. J. R. Stat. Soc. Series B 57, 289\u2013300 "
 "(1995).",
 "29. Gionis, A., Indyk, P. & Motwani, R. Similarity search in high dimensions via "
 "hashing. In Proc. 25th International Conference on Very Large Data Bases (VLDB), "
 "518\u2013529 (1999).",
]
anchor = None
for i, p in enumerate(d.paragraphs):
    if p.text.strip() == "Methods":
        anchor = p
        break
assert anchor is not None, "Methods heading not found"
style = d.paragraphs[104].style
for txt in refs:
    anchor.insert_paragraph_before(txt, style=style)
    print("ok ref:", txt[:20])

cp = d.core_properties
cp.author = "DAPR"
cp.last_modified_by = "DAPR"
d.save(PATH)
print("ALL REFS OK")
