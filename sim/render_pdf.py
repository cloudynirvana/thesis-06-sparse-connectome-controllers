#!/usr/bin/env python3
"""Render THESIS.md to a citeable PDF. Research manuscript only."""

from __future__ import annotations

import re
from pathlib import Path

from fpdf import FPDF

ROOT = Path(__file__).resolve().parents[1]
MD = ROOT / "THESIS.md"
PDF = ROOT / "THESIS.pdf"

SERIF = "/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf"
SERIF_B = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-Bold.ttf"
MONO = "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf"


class ThesisPDF(FPDF):
    def header(self):
        if self.page_no() == 1:
            return
        self.set_font("Serif", "I", 8)
        self.set_text_color(80, 80, 80)
        self.cell(
            0,
            8,
            "Ogbonna KE. Thesis #6 — sparse connectome-style controllers (research only; not a medical device)",
            new_x="LMARGIN",
            new_y="NEXT",
        )
        self.set_text_color(0, 0, 0)
        self.ln(2)

    def footer(self):
        self.set_y(-14)
        self.set_font("Serif", "I", 8)
        self.set_text_color(80, 80, 80)
        self.cell(0, 8, str(self.page_no()), align="C")
        self.set_text_color(0, 0, 0)


def tidy_math(text: str) -> str:
    text = text.replace(r"\bigl(", "(").replace(r"\bigr)", ")")
    text = text.replace(r"\,", " ")
    text = text.replace(r"\mathrm{", "").replace(r"\mathrm", "")
    text = text.replace(r"\mathrm{clip}", "clip")
    text = re.sub(r"\\mathrm\{([^}]+)\}", r"\1", text)
    text = text.replace(r"\mathrm", "")
    replacements = {
        r"\frac{dS}{dt}": "dS/dt",
        r"\frac{dR}{dt}": "dR/dt",
        r"\frac{S + \alpha_{SR} R}{K}": "(S + α_SR R)/K",
        r"\frac{R + \alpha_{RS} S}{K}": "(R + α_RS S)/K",
        r"\alpha_{SR}": "α_SR",
        r"\alpha_{RS}": "α_RS",
        r"\delta_S": "δ_S",
        r"\delta_R": "δ_R",
        r"r_S": "r_S",
        r"r_R": "r_R",
        r"T_{\mathrm{on}}": "T_on",
        r"T_{\mathrm{off}}": "T_off",
        r"u_{\mathrm{prev}}": "u_prev",
        r"n_{\mathrm{claw}}": "n_claw",
        r"\ell_2": "l2",
        r"\lambda": "λ",
        r"\ge": "≥",
        r"\le": "≤",
        r"\in": "∈",
        r"\times": "×",
        r"\top": "T",
        r"\Delta": "Δ",
        r"\left(": "(",
        r"\right)": ")",
        r"\left[": "[",
        r"\right]": "]",
        r"\[": "",
        r"\]": "",
        r"\(": "",
        r"\)": "",
        "{": "",
        "}": "",
    }
    for a, b in replacements.items():
        text = text.replace(a, b)
    text = text.replace("$", "")
    text = re.sub(r" +", " ", text)
    return text


def strip_md_row(row: str) -> list[str]:
    cells = [c.strip() for c in row.strip().strip("|").split("|")]
    return [re.sub(r"\*\*(.+?)\*\*", r"\1", c) for c in cells]


def write_inline(pdf: ThesisPDF, text: str, size: float = 11):
    text = tidy_math(text)
    parts = re.split(r"(\*\*[^*]+\*\*|`[^`]+`|\*[^*]+\*)", text)
    pdf.set_font("Serif", "", size)
    for part in parts:
        if not part:
            continue
        if part.startswith("**") and part.endswith("**"):
            pdf.set_font("Serif", "B", size)
            pdf.write(5.4, part[2:-2])
            pdf.set_font("Serif", "", size)
        elif part.startswith("`") and part.endswith("`"):
            pdf.set_font("Mono", "", size - 1)
            pdf.write(5.4, part[1:-1])
            pdf.set_font("Serif", "", size)
        elif part.startswith("*") and part.endswith("*") and len(part) > 2:
            pdf.set_font("Serif", "I", size)
            pdf.write(5.4, part[1:-1])
            pdf.set_font("Serif", "", size)
        else:
            pdf.set_font("Serif", "", size)
            pdf.write(5.4, part)
    pdf.ln(6.2)


def render():
    md = MD.read_text(encoding="utf-8")
    pdf = ThesisPDF(format="A4")
    pdf.set_auto_page_break(auto=True, margin=18)
    pdf.set_margins(18, 18, 18)
    pdf.add_font("Serif", "", SERIF)
    pdf.add_font("Serif", "B", SERIF_B)
    pdf.add_font("Serif", "I", SERIF)  # italic face not installed; regular fallback
    pdf.add_font("Serif", "BI", SERIF_B)
    pdf.add_font("Mono", "", MONO)
    pdf.add_page()
    pdf.set_title("Sparse Connectome-Style Controllers as In-Silico Policy Classes")
    pdf.set_author("Kelechi Emeka Ogbonna")

    lines = md.splitlines()
    i = 0
    in_code = False
    code_buf: list[str] = []
    table_buf: list[str] = []

    def flush_table():
        nonlocal table_buf
        if not table_buf:
            return
        rows = [strip_md_row(r) for r in table_buf if not re.match(r"^\|?\s*:?--", r)]
        table_buf = []
        if not rows:
            return
        ncol = max(len(r) for r in rows)
        for r in rows:
            while len(r) < ncol:
                r.append("")
        usable = pdf.epw
        col_w = usable / ncol
        pdf.set_font("Serif", "", 8.5)
        for ri, row in enumerate(rows):
            if ri == 0:
                pdf.set_font("Serif", "B", 8.5)
            else:
                pdf.set_font("Serif", "", 8.5)
            # estimate height
            line_h = 4.2
            max_lines = 1
            for cell in row:
                nlines = max(1, pdf.multi_cell(col_w, line_h, tidy_math(cell), dry_run=True, output="LINES").__len__() if False else 1)
            # simpler: multi_cell with ln
            x0 = pdf.l_margin
            y0 = pdf.get_y()
            heights = []
            for cell in row:
                # measure
                h = 4.2 * max(1, int(pdf.get_string_width(tidy_math(cell)) / max(col_w - 1, 1) + 1))
                heights.append(h)
            row_h = min(max(heights + [4.2]), 22)
            if y0 + row_h > pdf.h - pdf.b_margin:
                pdf.add_page()
                y0 = pdf.get_y()
            for ci, cell in enumerate(row):
                pdf.set_xy(x0 + ci * col_w, y0)
                pdf.multi_cell(col_w, 4.2, tidy_math(cell), border=1, max_line_height=4.2)
            pdf.set_y(y0 + row_h)
        pdf.ln(3)

    while i < len(lines):
        line = lines[i]
        if line.strip().startswith("```"):
            if in_code:
                pdf.set_font("Mono", "", 8)
                pdf.set_fill_color(245, 245, 245)
                for cl in code_buf:
                    pdf.multi_cell(0, 4, cl, fill=True)
                pdf.ln(2)
                code_buf = []
                in_code = False
            else:
                flush_table()
                in_code = True
            i += 1
            continue
        if in_code:
            code_buf.append(line)
            i += 1
            continue
        if line.strip().startswith("|"):
            table_buf.append(line)
            i += 1
            continue
        else:
            flush_table()

        if line.strip() == "---":
            pdf.ln(1)
            pdf.set_draw_color(160, 160, 160)
            pdf.line(pdf.l_margin, pdf.get_y(), pdf.w - pdf.r_margin, pdf.get_y())
            pdf.ln(4)
            i += 1
            continue
        if line.startswith("# "):
            pdf.ln(2)
            pdf.set_font("Serif", "B", 16)
            pdf.multi_cell(0, 7.5, line[2:].strip())
            pdf.ln(3)
            i += 1
            continue
        if line.startswith("## "):
            pdf.ln(3)
            pdf.set_font("Serif", "B", 13)
            pdf.multi_cell(0, 6.5, line[3:].strip())
            pdf.ln(2)
            i += 1
            continue
        if line.startswith("### "):
            pdf.ln(2)
            pdf.set_font("Serif", "B", 11.5)
            pdf.multi_cell(0, 6, line[4:].strip())
            pdf.ln(1)
            i += 1
            continue
        if line.startswith("- "):
            pdf.set_x(pdf.l_margin + 4)
            write_inline(pdf, "• " + line[2:].strip(), 11)
            i += 1
            continue
        if re.match(r"^\d+\. ", line):
            write_inline(pdf, line.strip(), 11)
            i += 1
            continue
        if line.strip() == "":
            pdf.ln(1.5)
            i += 1
            continue
        write_inline(pdf, line.strip(), 11)
        i += 1

    flush_table()
    pdf.output(str(PDF))
    print(f"Wrote {PDF} ({PDF.stat().st_size} bytes)")


if __name__ == "__main__":
    render()
