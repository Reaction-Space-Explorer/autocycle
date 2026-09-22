"""MANUSCRIPT.md to a submission .docx, laid out like the group's other papers.

A4, Roboto with a serif fallback, first-line indents on body paragraphs, bold
run-in headings rather than Word heading styles, and continuous line numbers for
review. Figures are placed inline where the text first calls them, each image
followed immediately by its bold caption, which is how the nucleoside-analogue
manuscript is built. Markdown pipe tables become real Word tables.
"""
import re
import sys
from pathlib import Path

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt

ROOT = Path(__file__).resolve().parents[1]
FIGS = {                      # numbered by order of first citation, not by filename
    1: ("fig1_ladder.png", 6.3),
    2: ("fig3_paired.png", 6.5),
    3: ("fig6_rule_removal.png", 5.8),
    4: ("fig_shared.png", 6.4),
    5: ("fig_bound.png", 6.1),
    6: ("fig4_triage_thermo.png", 6.3),
    7: ("fig5_depth.png", 6.3),
    8: ("fig_coresize.png", 5.6),
    9: ("fig7_calvin_core.png", 4.2),
    10: ("fig2_formose_core.png", 3.6),
}
BODY, SIZE = "Roboto", Pt(11)
LINE = 1.03                   # as the nucleoside-analogue manuscript is set


def line_numbers(section):
    ln = OxmlElement("w:lnNumType")
    ln.set(qn("w:countBy"), "1")
    ln.set(qn("w:restart"), "continuous")
    ln.set(qn("w:distance"), "360")
    section._sectPr.append(ln)


def style_run(r, bold=False, italic=False, mono=False, sup=False):
    r.bold, r.italic = bold, italic
    if sup:
        r.font.superscript = True
    r.font.name = "Courier New" if mono else BODY
    r.font.size = Pt(9.5) if mono else SIZE
    if not mono:                       # so the font survives on machines without Roboto
        r._element.rPr.rFonts.set(qn("w:cs"), BODY)


def runs(par, text, bold=False):
    for piece in re.split(r"(\*\*[^*]+\*\*|\*[^*]+\*|`[^`]+`|\^[^\s^]+\^)", text):
        if not piece:
            continue
        piece = piece.replace("\\*", "*")      # an escaped asterisk is a literal one
        if piece.startswith("^") and piece.endswith("^") and len(piece) > 2:
            style_run(par.add_run(piece[1:-1]), sup=True)
        elif piece.startswith("**"):
            style_run(par.add_run(piece[2:-2]), bold=True)
        elif piece.startswith("*"):
            style_run(par.add_run(piece[1:-1]), italic=True)
        elif piece.startswith("`"):
            style_run(par.add_run(piece[1:-1]), mono=True)
        else:
            style_run(par.add_run(piece), bold=bold)


def main(src="MANUSCRIPT.md", out=None):
    src = Path(src)
    if not src.is_absolute() and not src.exists():
        src = ROOT / src
    root = src.parent          # figures sit beside the manuscript
    out = Path(out) if out else src.with_suffix(".docx")
    text = src.read_text()
    doc = Document()
    n = doc.styles["Normal"]
    n.font.name, n.font.size = BODY, SIZE
    for s in doc.sections:
        s.page_width, s.page_height = Inches(8.27), Inches(11.69)
        s.left_margin = s.right_margin = Inches(1)
        s.top_margin = s.bottom_margin = Inches(1)
        line_numbers(s)

    captions = {}
    for m in re.finditer(r"\*\*Figure (\d+)\.\*\* (.+?)(?=\n\n)", text, re.S):
        captions[int(m.group(1))] = " ".join(m.group(0).split())
    placed = set()

    def para(t, bold=False, before=0, align=None, gap=False):
        """One paragraph, laid out as the nucleoside-analogue manuscript is.

        Justified, no first-line indent, line spacing 1.03, and paragraphs
        separated by a blank paragraph rather than by space after.
        """
        p = doc.add_paragraph()
        pf = p.paragraph_format
        pf.line_spacing = LINE
        pf.space_after = Pt(0)
        pf.space_before = Pt(before)
        p.alignment = align if align is not None else WD_ALIGN_PARAGRAPH.JUSTIFY
        runs(p, t, bold=bold)
        if gap:
            g = doc.add_paragraph()
            g.paragraph_format.line_spacing = LINE
            g.paragraph_format.space_after = Pt(0)
        return p

    def place_figure(k):
        """Image, then its caption, at the point the text first calls it."""
        name, width = FIGS[k]
        path = root / "figures" / name
        if not path.exists() or k in placed:
            return
        placed.add(k)
        holder = doc.add_paragraph()
        holder.alignment = WD_ALIGN_PARAGRAPH.CENTER
        holder.paragraph_format.space_before = Pt(10)
        holder.add_run().add_picture(str(path), width=Inches(width))
        cap = doc.add_paragraph()
        cap.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
        cap.paragraph_format.space_after = Pt(10)
        runs(cap, captions.get(k, f"**Figure {k}.**"))
        for r in cap.runs:
            r.font.size = Pt(9.5)

    blocks = re.split(r"\n\s*\n", text)
    in_captions = False
    for b in blocks:
        b = b.strip()
        if not b:
            continue
        if b.startswith("## Figure captions"):
            in_captions = True
            continue
        if in_captions and b.startswith("**Figure "):
            continue                    # captions travel with their figures
        if b.startswith("## "):
            in_captions = False
        if b.startswith("|"):
            rows = [r for r in b.split("\n") if r.strip().startswith("|")]
            cells = [[c.strip() for c in r.strip().strip("|").split("|")] for r in rows]
            cells = [c for c in cells if not set("".join(c)) <= set("-: ")]
            tab = doc.add_table(rows=0, cols=max(len(c) for c in cells))
            tab.style = "Table Grid"
            for j, row in enumerate(cells):
                wr = tab.add_row().cells
                for i, c in enumerate(row):
                    cp = wr[i].paragraphs[0]
                    cp.paragraph_format.space_after = Pt(0)
                    runs(cp, c, bold=(j == 0))
                    for r in cp.runs:
                        r.font.size = Pt(9)
            doc.add_paragraph()
            continue
        m = re.match(r"^(#{1,3}) (.+)$", b)
        if m:
            para(m.group(2), bold=True, before=12 if len(m.group(1)) > 1 else 0, gap=True)
            continue
        flat = b.replace("\n", " ")
        para(flat, gap=True)
        for k in sorted(FIGS):
            if re.search(rf"Figure {k}\b", flat) and k not in placed:
                place_figure(k)

    for k in sorted(FIGS):               # anything the text never called
        place_figure(k)
    doc.save(out)
    print(f"  {out}: {len(doc.paragraphs)} paragraphs, {len(doc.tables)} tables, "
          f"{len(doc.inline_shapes)} figures")


if __name__ == "__main__":
    main(*sys.argv[1:])
