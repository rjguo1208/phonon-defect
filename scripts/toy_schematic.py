#!/usr/bin/env python3
"""Schematic of the toy model for the toy-validation page (run: npm run schematic).

Compiles scripts/toy_schematic.tex (standalone class, TikZ) with pdflatex in a temporary
directory and writes site/results/toy-schematic.{pdf,svg,png}. The SVG and PNG are
rendered from the PDF with PyMuPDF; the SVG has its text converted to paths. Needs
pdflatex on the PATH (on Banff: TinyTeX in ~/.TinyTeX) and PyMuPDF.
"""
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pymupdf

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scripts" / "toy_schematic.tex"
OUT = ROOT / "site" / "results"
NAME = "toy-schematic"
DISPLAY_WIDTH = 1000                     # width attribute of the <img> on the Toy page


def main():
    pdflatex = shutil.which("pdflatex")
    if pdflatex is None:
        sys.exit("pdflatex not found")
    OUT.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        run = subprocess.run([pdflatex, "-interaction=nonstopmode", "-halt-on-error",
                              "-output-directory", tmp, str(SOURCE)],
                             capture_output=True, text=True)
        log = (Path(tmp) / (SOURCE.stem + ".log")).read_text(errors="replace")
        if run.returncode != 0:
            sys.exit(log[-3000:])
        problems = re.findall(r"^(?:LaTeX|Package \w+) Warning.*$|^(?:Overfull|Underfull) .*$", log, re.M)
        if problems:
            sys.exit("pdflatex warnings:\n" + "\n".join(problems))
        shutil.copyfile(Path(tmp) / (SOURCE.stem + ".pdf"), OUT / f"{NAME}.pdf")

    with pymupdf.open(OUT / f"{NAME}.pdf") as doc:
        page = doc[0]
        (OUT / f"{NAME}.svg").write_text(page.get_svg_image(text_as_path=True))
        page.get_pixmap(dpi=300).save(OUT / f"{NAME}.png")
        width, height = page.rect.width, page.rect.height
    print(f"wrote site/results/{NAME}.pdf, .svg, .png: {width / 72 * 2.54:.2f} x {height / 72 * 2.54:.2f} cm; "
          f'<img width="{DISPLAY_WIDTH}" height="{round(DISPLAY_WIDTH * height / width)}">')


if __name__ == "__main__":
    main()
