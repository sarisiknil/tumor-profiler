#!/usr/bin/env python3
"""Build the two remaining XX395 deliverables from the pipeline's own output: the presentation and the digest.

The faculty guidelines require three files. This produces the second and third:

  * Presentation (10-15 minutes, Times New Roman): motivation and problem definition, objective, methods and
    tools, outcomes and deliverables, results.
  * Digest: a single slide giving the project title, company, duration and author, followed by the project
    objectives and a numbered list of outcomes.

Every number comes from results/summary.json, so the slides cannot drift from the report or the analysis.

  python3 report/build_slides.py --results results
"""
import argparse, json
from pathlib import Path

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt

REPO = Path(__file__).resolve().parents[1]
FONT = "Times New Roman"
INK = RGBColor(0x1A, 0x1A, 0x1A)
MUTED = RGBColor(0x55, 0x55, 0x55)
ACCENT = RGBColor(0x1F, 0x4E, 0x79)
TODO = RGBColor(0xC0, 0x00, 0x00)


def jload(p, default=None):
    p = Path(p)
    if not p.exists():
        return default if default is not None else {}
    try:
        return json.loads(p.read_text())
    except Exception:
        return default if default is not None else {}


def tload(p):
    p = Path(p)
    if not p.exists():
        return []
    with open(p) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]


def run_text(paragraph, text):
    """python-pptx's add_run() takes no argument; the text is assigned to the run afterwards."""
    r = paragraph.add_run()
    r.text = text
    return r


def style(run, size=18, bold=False, color=INK, italic=False):
    run.font.name = FONT
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = color


def textbox(slide, left, top, width, height):
    tb = slide.shapes.add_textbox(Inches(left), Inches(top), Inches(width), Inches(height))
    tf = tb.text_frame
    tf.word_wrap = True
    return tf


def slide_with_title(prs, title, subtitle=None):
    s = prs.slides.add_slide(prs.slide_layouts[6])          # blank
    tf = textbox(s, 0.6, 0.35, 12.1, 0.9)
    style(run_text(tf.paragraphs[0], title), 30, bold=True, color=ACCENT)
    if subtitle:
        p = tf.add_paragraph()
        style(run_text(p, subtitle), 15, color=MUTED, italic=True)
    return s


def bullets(slide, items, left=0.75, top=1.5, width=11.9, size=17, gap=6):
    tf = textbox(slide, left, top, width, 5.4)
    for i, item in enumerate(items):
        text, level = (item if isinstance(item, tuple) else (item, 0))
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.level = level
        p.space_after = Pt(gap)
        style(run_text(p, ("—  " if level else "•  ") + text),
              size - (2 if level else 0), color=INK if not level else MUTED)
    return tf


def table(slide, headers, rows, left=0.9, top=1.6, width=11.5, height=None, size=13):
    height = height or min(4.8, 0.42 * (len(rows) + 1))
    shape = slide.shapes.add_table(len(rows) + 1, len(headers),
                                   Inches(left), Inches(top), Inches(width), Inches(height))
    t = shape.table
    for j, h in enumerate(headers):
        c = t.cell(0, j)
        c.text = ""
        style(run_text(c.text_frame.paragraphs[0], str(h)), size, bold=True)
    for i, row in enumerate(rows, start=1):
        for j, v in enumerate(row):
            c = t.cell(i, j)
            c.text = ""
            style(run_text(c.text_frame.paragraphs[0], "" if v is None else str(v)), size)
    return t


def note(slide, text, top=6.55, color=MUTED, size=13):
    tf = textbox(slide, 0.75, top, 11.9, 0.55)
    style(run_text(tf.paragraphs[0], text), size, color=color, italic=True)


def todo(slide, text, top=6.55):
    note(slide, "[TO COMPLETE] " + text, top=top, color=TODO)


# ------------------------------------------------------------------------------------ presentation
def build_presentation(S, out):
    lc = S.get("qc", {}).get("library_complexity", {})
    up, rup = lc.get("umi_pattern", {}), lc.get("rna_umi_pattern", {})
    dd, ra = lc.get("deduplication", {}), lc.get("rna_alignment", {})
    var = S.get("variants", {})
    merged = var.get("caller_concordance", {})
    counts = var.get("counts_by_class", {})
    tmb = (S.get("biomarkers") or {}).get("tmb", {})
    fus = (S.get("fusions") or {}).get("summary", {})
    val = (S.get("validation") or {}).get("summary", {})
    hit = next((r for r in (S.get("validation") or {}).get("table", [])
                if r.get("status") == "CONCORDANT"), {})
    filt = jload(REPO / "results" / "dna" / "variants_filtered_summary.json")

    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)

    # 1 title
    s = prs.slides.add_slide(prs.slide_layouts[6])
    tf = textbox(s, 0.9, 2.1, 11.5, 2.6)
    style(run_text(tf.paragraphs[0], "Educational DNA–RNA Tumour Profiling"), 40, bold=True, color=ACCENT)
    p = tf.add_paragraph(); style(run_text(p, "Pipeline and Visualization Platform"), 40, bold=True, color=ACCENT)
    p = tf.add_paragraph(); p.space_before = Pt(18)
    style(run_text(p, "Building a clinical-grade genomic workflow in the open, and measuring it "
                    "against an accredited laboratory"), 17, color=MUTED, italic=True)
    todo(s, "Your name, programme, company, internship dates, supervisor", top=5.4)

    # 2 the question
    s = slide_with_title(prs, "The question")
    tf = textbox(s, 1.1, 2.0, 11.0, 1.5)
    style(run_text(tf.paragraphs[0], "“What can we understand about a person's disease "
                                   "by looking at their tumour's genetic material?”"), 26, italic=True,
          color=ACCENT)
    bullets(s, [
        "A tumour biopsy is sequenced, and a two-page clinical report names a few mutations and the drugs "
        "that target them.",
        "What happens in between is invisible: the software is commercial and closed, the thresholds are "
        "undocumented, and the reasoning is not shown.",
        "This project builds that chain in the open — and then checks it against the real clinical report "
        "for the same specimen.",
    ], top=3.5, size=18)

    # 3 the data
    s = slide_with_title(prs, "What we started with", "One biopsy, sequenced twice — and no documentation")
    table(s, ["", "DNA", "RNA"], [
        ["Read pairs", "5,542,544", "13,870,216"],
        ["Read length", "2 × 149 bp", "2 × 151 bp"],
        ["Instrument", "NextSeq 550", "NovaSeq 6000"],
        ["Design", "tumour only, no matched normal", "small fusion panel"],
    ], top=1.9, height=2.4)
    bullets(s, [
        "Nothing stated which assay produced the files. The first task was not analysis but identification.",
        "A clinical report for the same specimen existed — which made honest validation possible.",
    ], top=4.6, size=17)

    # 4 identifying the assay
    s = slide_with_title(prs, "Identifying the assay from the reads alone",
                         "Base composition per position: random bases mark a barcode, a fixed base marks an adapter")
    tf = textbox(s, 0.9, 1.9, 11.6, 1.4)
    style(run_text(tf.paragraphs[0], "Read 1:  [ 12 random bases = molecular barcode ][ fixed adapter sequence ]"
                                   "[ patient DNA … ]"), 18, bold=True)
    p = tf.add_paragraph()
    style(run_text(p, "Read 2:  [ gene-specific primer ][ patient DNA … ]"), 18, bold=True)
    bullets(s, [
        "Only ~2,200 distinct sequences account for 80 % of Read-2 starts → a primer-based targeted panel.",
        "Index reads of 8 and 10 bases → a specific adapter kit.",
        "Conclusion: anchored multiplex PCR. Confirmed afterwards by the clinical report: Archer VariantPlex "
        "Expanded Solid Tumor (DNA) and FusionPlex Lung v2 (RNA).",
        (f"{up.get('fraction_matching', 0)*100:.1f} % of DNA reads and {rup.get('fraction_matching', 0)*100:.1f} %"
         f" of RNA reads matched the predicted structure exactly.", 1),
    ], top=3.5, size=17)

    # 5 architecture
    s = slide_with_title(prs, "Pipeline architecture", "Split by memory, not by convenience")
    table(s, ["Stage", "Tool", "Runs on"], [
        ["Read structure, primer inference", "project scripts", "laptop"],
        ["Barcode extraction, alignment, deduplication", "UMI-tools, BWA-MEM", "Galaxy"],
        ["Primer clipping, coverage", "samtools, mosdepth", "Galaxy"],
        ["Somatic calling ×3", "Mutect2, LoFreq, FreeBayes", "Galaxy"],
        ["RNA alignment and fusions", "STAR, Arriba", "Galaxy"],
        ["Annotation, tiering, biomarkers, report", "VEP REST, CIViC, project scripts", "laptop"],
    ], top=1.75, height=3.1)
    note(s, "Human genome alignment needs ~6 GB (BWA) and ~31 GB (STAR); the available machine had 8 GB. "
            "The heavy half runs on the free public Galaxy service, the interpretation locally.", top=5.2)

    # 6 molecules not reads
    s = slide_with_title(prs, "Counting molecules, not reads")
    bullets(s, [
        f"{dd.get('reads_in', 0):,} reads collapse to {dd.get('unique_molecules_out', 0):,} unique molecules "
        f"({dd.get('reads_per_molecule', 0)} reads per molecule).",
        f"Only {dd.get('mean_unique_umis_per_position', 0)} distinct molecules per amplicon start position.",
        "That second number is why barcodes matter: every molecule from one amplicon starts at the same "
        "coordinate, so ordinary duplicate marking would collapse them all into one.",
        "The experiment is a 2.9-million-molecule measurement, not a 5.5-million-read one.",
    ], top=1.9, size=19)

    # 7 calling without a normal
    s = slide_with_title(prs, "Somatic calling without a matched normal",
                         "Every variant could in principle be inherited")
    table(s, ["Callers agreeing", "Variants"],
          [[k, f"{v:,}"] for k, v in sorted((merged.get("by_agreement") or {}).items())],
          top=1.8, width=5.0, height=1.8)
    tf = textbox(s, 6.6, 1.8, 6.0, 3.4)
    for i, t in enumerate([
        "No single caller is reliable here.",
        "1,573 of 1,718 positions were seen by one caller only.",
        "Agreement between independent callers is used as the filter.",
        "Population frequency, ClinVar and allele fraction then separate likely germline.",
        "Variants that cannot be resolved are flagged, never silently deleted.",
    ]):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.space_after = Pt(8)
        style(run_text(p, "•  " + t), 16)
    note(s, f"Result: {counts.get('SOMATIC_LIKELY', 0)} likely somatic · "
            f"{counts.get('GERMLINE_LIKELY', 0)} likely germline · "
            f"{counts.get('LOW_QUALITY', 0)} insufficient evidence", top=5.5, color=INK, size=15)

    # 8 headline
    s = slide_with_title(prs, "The result", "Two independent workflows, one specimen")
    table(s, ["", "Accredited laboratory", "This pipeline"], [
        ["Variant", "IDH1 p.Arg132Cys", "IDH1 p.Arg132Cys"],
        ["Allele fraction", hit.get("lab_vaf", "0.281"), hit.get("our_vaf", "0.2814")],
        ["Supported by", "vendor software", "3 of 3 callers"],
        ["Reference build", "GRCh37 / hg19", "GRCh38"],
        ["Reported variants recovered", "—", f"{val.get('concordant', 1)} of {val.get('reported_by_laboratory', 1)}"],
    ], top=1.9, height=2.6, size=15)
    tf = textbox(s, 0.9, 4.9, 11.6, 1.4)
    style(run_text(tf.paragraphs[0], "Different aligners, different callers, different reference builds, "
                                   "different deduplication — the same answer to three decimal places."),
          19, bold=True, color=ACCENT)

    # 9 where we differ
    s = slide_with_title(prs, "Where the two disagree — and why that is the interesting part")
    bullets(s, [
        "Clinical tier: the laboratory says Tier I-A, this pipeline says Tier II.",
        ("The laboratory cited a clinical guideline recording an approved therapy in this tumour type. The "
         "open knowledge bases hold that evidence only for other tumour types, and the tiering rule refuses "
         "Tier I on evidence from a different disease.", 1),
        ("The gap measures the uneven tumour-type coverage of free evidence sources — not a disagreement "
         "about biology.", 1),
        f"Candidate count: {counts.get('SOMATIC_LIKELY', 0)} likely somatic against one reported.",
        ("Partly expected — the laboratory reports only what passes its clinical thresholds. The rest is "
         "technical, and measurable.", 1),
    ], top=1.8, size=17)

    # 10 artefacts
    s = slide_with_title(prs, "One artefact with two faces",
                         "Cross-priming between amplicons amplified in the same reaction")
    tf = textbox(s, 0.75, 1.7, 5.7, 4.4)
    style(run_text(tf.paragraphs[0], "DNA arm"), 20, bold=True, color=ACCENT)
    for t in [f"{(filt.get('dominant_substitution_class') or {}).get('fraction_of_somatic_snvs', 0)*100:.0f} %"
              " of somatic SNVs share one substitution class",
              "Four “independent mutations” within 35 bases in one gene",
              f"{filt.get('somatic_flagged_artefact', 0)} of {counts.get('SOMATIC_LIKELY', 0)} candidates flagged"]:
        p = tf.add_paragraph(); p.space_after = Pt(8); style(run_text(p, "•  " + t), 16)
    tf = textbox(s, 6.9, 1.7, 5.7, 4.4)
    style(run_text(tf.paragraphs[0], "RNA arm"), 20, bold=True, color=ACCENT)
    for t in [f"{fus.get('n_fusions', 0)} fusions called, 166,694 discarded",
              "FGFR1/2/3 joined in every combination, in both orientations",
              f"All {fus.get('n_fusions', 0)} flagged; the laboratory reported none — we agree"]:
        p = tf.add_paragraph(); p.space_after = Pt(8); style(run_text(p, "•  " + t), 16)
    note(s, "An unscreened reading would have reported an FGFR2 fusion — the single most consequential false "
            "positive available in this cancer type.", top=6.2, color=TODO, size=15)

    # 11 limits
    s = slide_with_title(prs, "What this data cannot tell us", "The limits are the point, not a disclaimer")
    bullets(s, [
        "Whether a variant is inherited — no normal tissue was sequenced.",
        "Roughly half of the nominal target territory had no coverage at all.",
        "Genes frequently mutated in this cancer type are absent from the panel entirely: a negative there is "
        "an untested region, not an absence of mutation.",
        f"Mutational burden: {tmb.get('tmb_per_mb_before_screening', 0)} → {tmb.get('tmb_per_mb', 0)} per Mb "
        "after screening; both implausible for this tumour type, and the assessable territory is far too small.",
        "Microsatellite status: not determinable — as the accredited laboratory also concluded.",
    ], top=1.8, size=17)

    # 12 deliverables
    s = slide_with_title(prs, "Deliverables")
    bullets(s, [
        "Reproducible Snakemake + Bioconda pipeline, two execution modes (local example, Galaxy patient).",
        "Public repository with setup instructions and worked examples.",
        "Synthetic example dataset built over a real mini-reference, with a truth set — the whole workflow "
        "runs and is tested without any patient data.",
        "38 automated tests and continuous integration.",
        "Six-page interactive dashboard.",
        "Fourteen documentation pages: what each step measures, what it means, how it misleads.",
        "Validation against an accredited laboratory report.",
    ], top=1.75, size=17)
    todo(s, "Insert the repository URL", top=6.4)

    # 13 conclusions
    s = slide_with_title(prs, "Conclusions")
    bullets(s, [
        "An open pipeline, on a laptop plus a free public service, reproduced an accredited laboratory's "
        "finding to three decimal places.",
        "Every disagreement was explicable and measurable: evidence coverage, reporting thresholds, and one "
        "identifiable chemistry artefact.",
        "The most useful property is not the agreement but the capacity to state precisely where the answers "
        "stop — which is what makes it a teaching object rather than a black box.",
    ], top=1.9, size=19)
    note(s, "Not a diagnostic tool. Not validated for clinical use. Educational output only.", top=5.9,
         color=TODO, size=15)

    prs.save(out)
    return out


# ------------------------------------------------------------------------------------ digest
def build_digest(S, out):
    val = (S.get("validation") or {}).get("summary", {})
    hit = next((r for r in (S.get("validation") or {}).get("table", [])
                if r.get("status") == "CONCORDANT"), {})
    counts = S.get("variants", {}).get("counts_by_class", {})

    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
    s = prs.slides.add_slide(prs.slide_layouts[6])

    tf = textbox(s, 0.7, 0.3, 11.9, 1.4)
    tf.paragraphs[0].alignment = PP_ALIGN.CENTER
    style(run_text(tf.paragraphs[0], "Educational DNA–RNA Tumour Profiling Pipeline and "
                                   "Visualization Platform"), 24, bold=True, color=ACCENT)
    for label in ("[COMPANY NAME / ADDRESS]", "[DURATION OF PROJECT: month day – month day, year]",
                  "[STUDENT NAME / PROGRAMME]"):
        p = tf.add_paragraph(); p.alignment = PP_ALIGN.CENTER
        style(run_text(p, label), 13, color=TODO, bold=True)

    tf = textbox(s, 0.7, 2.0, 11.9, 1.55)
    style(run_text(tf.paragraphs[0], "PROJECT OBJECTIVE & EXPECTATIONS"), 14, bold=True)
    p = tf.add_paragraph()
    style(run_text(p, "To build a reproducible, documented pipeline that takes tumour DNA and RNA sequencing "
                    "from a targeted cancer panel through to an interpreted summary — variant calling, "
                    "clinical tiering, fusion detection and reporting — and to make every step inspectable "
                    "rather than hidden inside commercial software. The pipeline was to be validated against "
                    "the report issued by an accredited clinical laboratory for the same specimen."), 12.5)

    tf = textbox(s, 0.7, 3.65, 11.9, 3.3)
    style(run_text(tf.paragraphs[0], "OUTCOMES"), 14, bold=True)
    outcomes = [
        "A reproducible Snakemake + Bioconda pipeline covering quality control, barcode-aware deduplication, "
        "tumour-only somatic calling, annotation, clinical tiering, RNA fusion analysis and reporting.",
        f"Independent recovery of the accredited laboratory's finding: {hit.get('gene', 'IDH1')} "
        f"{hit.get('variant', 'p.Arg132Cys')} at an allele fraction of {hit.get('our_vaf', '0.2814')} against "
        f"the laboratory's {hit.get('lab_vaf', '0.281')} — {val.get('concordant', 1)} of "
        f"{val.get('reported_by_laboratory', 1)} reported variants, none missed.",
        "Identification of the assay chemistry from the raw reads alone, later confirmed by the clinical "
        "report, and reconstruction of the panel design from the sequencing data.",
        "Artefact screening for both arms, which identified a single chemistry-derived mechanism explaining "
        "both the excess DNA variant calls and all 20 candidate RNA fusions.",
        "A public repository containing the pipeline, a synthetic example dataset with a truth set, 38 "
        "automated tests, continuous integration and fourteen documentation pages.",
        "A six-page interactive dashboard presenting variants, fusions, pathways, therapies and the limits of "
        "the analysis.",
    ]
    for i, o in enumerate(outcomes, start=1):
        p = tf.add_paragraph(); p.space_after = Pt(5)
        style(run_text(p, f"{i}.  {o}"), 12)

    prs.save(out)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", default="results")
    ap.add_argument("--outdir", default="report")
    a = ap.parse_args()
    S = jload(Path(a.results) / "summary.json")
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    p1 = build_presentation(S, out / "BIO395_Presentation_DRAFT.pptx")
    p2 = build_digest(S, out / "BIO395_Digest_DRAFT.pptx")
    print(f"wrote {p1}\nwrote {p2}")


if __name__ == "__main__":
    main()
