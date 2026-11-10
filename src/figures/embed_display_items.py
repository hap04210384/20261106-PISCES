# -*- coding: utf-8 -*-
"""Embed Figures 2-4 (EMF) and journal-style Tables 2-6 into Manuscript.docx.
Display items go into the single-column stretch that already holds Table 1
(before the sectPr paragraph preceding 'Methods'). Prose is updated to cite them.
"""
import docx
from docx.oxml.ns import qn, nsdecls
from docx.oxml import parse_xml, OxmlElement
from docx.opc.part import Part
from docx.opc.packuri import PackURI
from docx.opc.constants import RELATIONSHIP_TYPE as RT
import os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DOC = ROOT + r'\Manuscript.docx'
FIG = ROOT + r'\figures\emf'

XMLNS = ('xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" '
         'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" '
         'xmlns:wp="http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing" '
         'xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" '
         'xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture"')
EMU_MM = 36000
TW_PER_MM = 56.6929

doc = docx.Document(DOC)
_DOC_ID = [1000]


def new_para_xml(runs, center=False, after_pt=3.0, before_pt=0.0):
    """runs: list of (text, bold). Returns w:p element, Times New Roman 7pt."""
    jc = '<w:jc w:val="center"/>' if center else ''
    rxml = ''
    for text, bold in runs:
        rpr = ('<w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"/>'
               '<w:sz w:val="14"/><w:szCs w:val="14"/>' + ('<w:b/>' if bold else '') + '</w:rPr>')
        rxml += '<w:r>%s<w:t xml:space="preserve">%s</w:t></w:r>' % (rpr, text.replace('&', '&amp;').replace('<', '&lt;'))
    return parse_xml(
        '<w:p %s><w:pPr><w:spacing w:after="%d" w:before="%d"/>%s</w:pPr>%s</w:p>'
        % (nsdecls('w'), int(after_pt * 20), int(before_pt * 20), jc, rxml))


def add_emf_para(path, w_mm, h_mm, center=False):
    blob = open(path, 'rb').read()
    existing = {str(p.partname) for p in doc.part.package.iter_parts()}
    i = 1
    while ('/word/media/image%d.emf' % i) in existing:
        i += 1
    part = Part(PackURI('/word/media/image%d.emf' % i), 'image/x-emf', blob, doc.part.package)
    rId = doc.part.relate_to(part, RT.IMAGE)
    cx, cy = int(w_mm * EMU_MM), int(h_mm * EMU_MM)
    _DOC_ID[0] += 1
    xml = (
        '<w:p %s><w:pPr><w:spacing w:after="40"/></w:pPr><w:r><w:drawing>'
        '<wp:inline distT="0" distB="0" distL="0" distR="0">'
        '<wp:extent cx="%d" cy="%d"/><wp:effectExtent l="0" t="0" r="0" b="0"/>'
        '<wp:docPr id="%d" name="fig%d"/><wp:cNvGraphicFramePr/>'
        '<a:graphic><a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture">'
        '<pic:pic><pic:nvPicPr><pic:cNvPr id="0" name="fig%d.emf"/><pic:cNvPicPr/></pic:nvPicPr>'
        '<pic:blipFill><a:blip r:embed="%s"/><a:stretch><a:fillRect/></a:stretch></pic:blipFill>'
        '<pic:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="%d" cy="%d"/></a:xfrm>'
        '<a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr></pic:pic>'
        '</a:graphicData></a:graphic></wp:inline></w:drawing></w:r></w:p>'
    ) % (XMLNS, cx, cy, _DOC_ID[0], i, i, rId, cx, cy)
    return parse_xml(xml)


def make_nm_table(header, rows, widths_mm):
    n = len(header)
    total_dxa = int(sum(widths_mm) * TW_PER_MM)
    tbl = OxmlElement('w:tbl')
    tblPr = OxmlElement('w:tblPr')
    tblPr.append(parse_xml('<w:tblW %s w:w="%d" w:type="dxa"/>' % (nsdecls('w'), total_dxa)))
    jc = OxmlElement('w:jc'); jc.set(qn('w:val'), 'center'); tblPr.append(jc)
    borders = parse_xml(
        '<w:tblBorders %s>'
        '<w:top w:val="single" w:sz="12" w:space="0" w:color="000000"/>'
        '<w:left w:val="none" w:sz="0" w:space="0" w:color="auto"/>'
        '<w:bottom w:val="single" w:sz="12" w:space="0" w:color="000000"/>'
        '<w:right w:val="none" w:sz="0" w:space="0" w:color="auto"/>'
        '<w:insideH w:val="none" w:sz="0" w:space="0" w:color="auto"/>'
        '<w:insideV w:val="none" w:sz="0" w:space="0" w:color="auto"/>'
        '</w:tblBorders>' % nsdecls('w'))
    tblPr.append(borders)
    layout = OxmlElement('w:tblLayout'); layout.set(qn('w:type'), 'fixed'); tblPr.append(layout)
    cellmar = parse_xml(
        '<w:tblCellMar %s><w:top w:w="28" w:type="dxa"/><w:left w:w="57" w:type="dxa"/>'
        '<w:bottom w:w="28" w:type="dxa"/><w:right w:w="57" w:type="dxa"/></w:tblCellMar>' % nsdecls('w'))
    tblPr.append(cellmar)
    tbl.append(tblPr)
    grid = OxmlElement('w:tblGrid')
    for w in widths_mm:
        gc = OxmlElement('w:gridCol'); gc.set(qn('w:w'), str(int(w * TW_PER_MM)))
        grid.append(gc)
    tbl.append(grid)

    def row_el(texts, bold, header_row=False):
        tr = OxmlElement('w:tr')
        if header_row:
            trPr = OxmlElement('w:trPr')
            th = OxmlElement('w:tblHeader'); trPr.append(th)
            tr.append(trPr)
        for text, w in zip(texts, widths_mm):
            tc = OxmlElement('w:tc')
            tcPr = OxmlElement('w:tcPr')
            tcW = OxmlElement('w:tcW'); tcW.set(qn('w:w'), str(int(w * TW_PER_MM))); tcW.set(qn('w:type'), 'dxa')
            tcPr.append(tcW)
            if header_row:
                tcPr.append(parse_xml('<w:tcBorders %s><w:bottom w:val="single" w:sz="6" w:space="0" w:color="000000"/></w:tcBorders>' % nsdecls('w')))
            tc.append(tcPr)
            tc.append(new_para_xml([(text, bold)], after_pt=1.0, before_pt=1.0))
            tr.append(tc)
        return tr

    tbl.append(row_el(header, True, header_row=True))
    for r in rows:
        tbl.append(row_el(r, False))
    return tbl


# ---------------- locate insertion anchor ----------------
anchor = None
for p in doc.paragraphs:
    if 'sectPr' in p._p.xml:
        nxt = p._p.getnext()
        while nxt is not None and nxt.tag == qn('w:p') and not nxt.findall('.//' + qn('w:t')):
            nxt = nxt.getnext()
        if nxt is not None:
            txt = ''.join(t.text or '' for t in nxt.findall('.//' + qn('w:t')))
            if txt.strip() == 'Methods':
                anchor = p
                break
assert anchor is not None, 'anchor sectPr paragraph not found'
blocks = []


def caption(label, text):
    return new_para_xml([(label + ' ', True), (text, False)], after_pt=3.0, before_pt=6.0)


def note(text):
    return new_para_xml([(text, False)], after_pt=8.0)

# ---------------- Figures ----------------
blocks.append(add_emf_para(FIG + r'\fig2_scaling.emf', 183, 62))
blocks.append(caption('Fig. 2 |', 'Scaling behaviour of exact pair counting. a, Wall time (log scale) for the full '
    'M x M pair count (M = 54,590 gene families; 1.49 billion pairs) on one AVX2 CPU core and on the GPU BMMA wave '
    'kernel, at N = 1,992 (real K. pneumoniae pangenome), 4,000 and 8,000 genomes; GPU output is byte-identical to '
    'the CPU reference at every point. End-to-end GPU time including screened-pair I/O is 9.9 s, 1.4 s and 1.4 s. '
    'b, Throughput (pairs per second, log scale): CPU throughput halves as N doubles while the GPU kernel holds '
    '~2.2 billion pairs per second.'))
blocks.append(add_emf_para(FIG + r'\fig3_planted.emf', 183, 60))
blocks.append(caption('Fig. 3 |', 'Recall and power limits on planted ground truth. a, Recall at conservative '
    'BH-FDR 0.05 versus planted support bin for 1,000 interactions planted into the real K. pneumoniae bitmap '
    '(400 co-occurrence, 300 conditional, 300 exclusion). b, Overall recall is insensitive to the FDR level '
    '(0.288 to 0.272 from FDR 0.1 to 0.001) but drops under support floors (c11 >= 10 or minsup >= 20), showing '
    'that discovery is support-limited rather than threshold-limited.'))
blocks.append(add_emf_para(FIG + r'\fig4_baselines.emf', 183, 58))
blocks.append(caption('Fig. 4 |', 'Approximate pipelines miss rare, high-effect interactions. a, Recovery of the '
    'z-screened pair set by the exact pipeline and by MinHash-LSH (128 hash functions, 32 x 4 banding) on the real '
    'K. pneumoniae bitmap; exclusion pairs are nearly invisible to Jaccard-LSH. b, Marginal-support pre-screening '
    'trades tested pairs against retained signal: even the most permissive floor (minsup >= 5) must test 16.7% of '
    'all pairs yet loses 61% of the deviating pairs.'))

# ---------------- Tables ----------------
W3 = [52, 45, 85]
blocks.append(caption('Table 2 |', 'Statistical calibration and multiplicity basis.'))
blocks.append(make_nm_table(
    ['Quantity', 'Result', 'Note'],
    [['Screen rate on iid null bitmap', '1.71M emitted vs 1.49M expected', 'chance level; nominal rate 1e-3 x C(M,2)'],
     ['BH-FDR within screened set (m = 22.2M)', 'empirical FPR = 1.000 (0.1 / 0.05 / 0.01)', 'no false-positive control'],
     ['BH-FDR on C(M,2) (conservative)', 'empirical FPR = 0.000 (0.1 / 0.05 / 0.01 / 0.001)', 'valid; basis for all reported significance'],
     ['Significant pairs, conservative BH', '2.97M / 2.69M / 2.21M / 1.77M', 'FDR 0.1 / 0.05 / 0.01 / 0.001'],
     ['Multiplicity by unique 2x2 patterns (18.8M)', 'screened-set BH still flags 100% at 0.05', 'screen is data-dependent; pattern count does not rescue']],
    W3))
blocks.append(note('FPR: fraction of null (iid Bernoulli) instances called significant; BH inside the screened set '
    'has no control because the screen itself is a function of the data, so all significance is reported against '
    'the full pair count C(M,2).'))

blocks.append(caption('Table 3 |', 'Planted-benchmark recall (1,000 interactions; conservative BH).'))
blocks.append(make_nm_table(
    ['Type', 'n', 'Recall @0.05', '@0.01', '@0.001', 'Median OR', 'Detection behaviour'],
    [['co-occurrence', '400', '46.7%', '46.8%', '-', '1.46e5', 'support >= 6 essentially always recovered'],
     ['conditional', '300', '28.7%', '27.7%', '-', '124', 'recall rises with planted support'],
     ['exclusion', '300', '4.7%', '3.0%', '-', '0', 'detectable only when s^2/N >= 5 (power floor)']],
    [34, 12, 26, 20, 20, 24, 46]))
blocks.append(note('Overall recall 287/1,000 (28.7%); interactions planted at log-uniform support (1e-3 to 1e-1) into '
    'the real bitmap. Recall varies by at most 1.6 percentage points across FDR 0.1-0.001 but falls to 0.209 (c11 >= 10) '
    'and 0.184 (minsup >= 20) under support floors (Fig. 3b).'))

W4 = [34, 40, 44, 40, 24]
blocks.append(caption('Table 4 |', 'Approximate baselines on the real K. pneumoniae bitmap.'))
blocks.append(make_nm_table(
    ['Method', 'Configuration', 'Cost', 'Recovery of z-screened pairs', 'Consequence'],
    [['MinHash-LSH', '128 hashes, 32 x 4 bands', '14.2M candidates in 40 s', '12.75% (co-occ. 15.1%, exclusion 0.69%)', 'blind to exclusion structure'],
     ['Marginal pre-screen', 'minsup >= 5', '249M pairs tested (16.7%)', 'retains 38.5% (loses 61%)', 'dominated on both axes'],
     ['Marginal pre-screen', 'minsup >= 100', '28.5M pairs tested (1.9%)', 'retains 8.1% (loses 92%)', 'dominated on both axes'],
     ['Split-half replication', 'two 996-genome cohorts', '-', 'Jaccard 0.504; top-100k overlap 83.0%', 'stable structure, not cohort noise']],
    W4))
blocks.append(note('Exact bitmap counting with a deviation screen is both faster-to-candidate and more complete than '
    'either approximate route in this regime (Fig. 4).'))

blocks.append(caption('Table 5 |', 'K. pneumoniae pangenome dataset.'))
blocks.append(make_nm_table(
    ['Item', 'Value'],
    [['Source', 'PanGBank pangenome 11756, GTDB_refseq collection 2.0.0'],
     ['Genomes', '2,001 total; 1,992 post-QC (completeness >= 90%, contamination <= 5%)'],
     ['Assembly quality (CheckM2)', 'completeness 99.41% mean; contamination 0.17% mean'],
     ['Gene families', '54,590 (Panaroo)'],
     ['Presences', '10,452,880'],
     ['All feature pairs', '1.49 x 10^9'],
     ['Screened pairs (|z| >= 3.29)', '22,155,477'],
     ['Significant pairs (BH < 0.05)', '2,685,180 (2.21M at FDR 0.01)']],
    [70, 112]))
blocks.append(note('QC filters applied per genome before bitmap construction; see Methods for the manifest and '
    'checksums.'))

blocks.append(caption('Table 6 |', 'Flagship interacting locus: the yersiniabactin (ybt) siderophore system.'))
blocks.append(make_nm_table(
    ['Family A', 'Family B', 'c11', 'supA', 'supB', 'OR', 'padj'],
    [['AB8Z72_RS09615 (ybtT, type II thioesterase)', 'AB8Z72_RS09630 (irp2, yersiniabactin biosynthesis)', '476', '512', '483', '2.56e3', '<1e-300'],
     ['AB8Z72_RS09615 (ybtT)', 'AB8Z72_RS09620 (ybtU, yersiniabactin biosynthesis)', '511', '512', '513', '2.02e5', '<1e-300'],
     ['AB8Z72_RS09615 (ybtT)', 'AB8Z72_RS09640 (ybtP, ybt ABC transporter)', '512', '512', '513', '1.01e6', '<1e-300']],
    [48, 52, 14, 14, 14, 16, 24]))
blocks.append(note('The locus anchor family (ybtT; present in 512 of 1,992 genomes) is recovered together with its '
    'biosynthetic and transport partners at odds ratios up to 10^6; VFDB annotations in parentheses.'))

for b in blocks:
    anchor._p.addprevious(b)

# ---------------- prose updates ----------------
def para_by_prefix(prefix):
    for p in doc.paragraphs:
        if p.text.strip().startswith(prefix):
            return p
    raise KeyError(prefix)

def append_to_para(p, sentence):
    r = p.add_run(' ' + sentence)
    r.font.name = 'Times New Roman'
    r.font.size = docx.shared.Pt(8.5)

append_to_para(para_by_prefix('A calibration experiment on pure null data'),
    'Table 2 compiles these calibration results and the multiplicity analysis.')

# rewrite the scaling paragraph with the GPU ladder numbers
p = para_by_para = para_by_prefix('Scaling was measured on synthetic null bitmaps')
for r in list(p.runs):
    r.text = ''
p.runs[0].text = ('Scaling was measured at fixed M = 54,590 families (1.49 billion pairs) with N = 1,992 (real '
    'K. pneumoniae pangenome), 4,000 and 8,000 genomes. The GPU BMMA kernel counts all pairs in 0.68 s '
    'independently of N (CPU AVX2: 28.9 s, 46.8 s, 90.6 s; 43x to 130x kernel speed-up, byte-identical output), and '
    'end-to-end GPU time including screened-pair I/O is 9.9 s, 1.4 s and 1.4 s (Fig. 2). CPU throughput halves as '
    'N doubles (51.5 to 16.5 Mpair/s) while the GPU kernel holds ~2.2 billion pairs per second, so at fixed feature '
    'count the exact pipeline is nearly free in the number of genomes on the matrix engine; the real, heavily '
    'skewed bitmap is the fastest CPU instance, indicating that output size and support structure, not nominal pair '
    'count, dominate cost at this scale. The full GTDB ladder (1e5 to 1e6+ genomes) and datacenter-GPU energy/cost '
    'figures remain pending and will be added; the commodity-GPU reference already covers the biologically relevant '
    'single-species regime in under a minute.')

append_to_para(para_by_prefix('Ground truth was planted directly'),
    'Figure 3a and Table 3 give the support-resolved recall and the FDR and support-floor sensitivity.')
append_to_para(para_by_prefix('Two approximate alternatives were quantified'),
    'Figure 4 and Table 4 summarize both comparisons.')
append_to_para(para_by_prefix('Applying the pipeline to the quality-filtered'),
    'Table 5 summarizes the dataset. Among the top-ranked pairs, the complete yersiniabactin (ybt) siderophore '
    'locus emerges as a coherent module: its anchor family is recovered together with its biosynthetic and '
    'transport partners at odds ratios up to 10^6 (Table 6).')

doc.core_properties.author = 'DAPR'
doc.core_properties.last_modified_by = 'DAPR'
doc.save(DOC)
print('embedded: %d blocks' % len(blocks))
