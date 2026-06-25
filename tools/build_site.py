#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_site.py  —  Static-site generator for the phonon-defect repo.

Converts Markdown notes in content/ into MathJax-rendered sub-pages under
docs/pages/ and regenerates the landing page (docs/index.html) with the
catalog table. Converter copied from the claude-sternheimer repo (math-
protected, stdlib only); the driver below is registry-based: to add a page,
drop a Markdown file in content/ and register it in PAGES + CATALOG.

Usage
-----
    python3 tools/build_site.py
"""
import os, re, html

ROOT      = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCS      = os.path.join(ROOT, "docs")
PAGES_DIR = os.path.join(DOCS, "pages")

GEN_DATE   = "2026-06-10"
SITE_TITLE = "Phonon-Defect"
BRAND      = "Phonon&nbsp;Defect"

# ======================================================================
#  Markdown -> HTML  (math-protected, stdlib only)
# ======================================================================
NUL = "\x00"

def _protect(md, store):
    """Replace fenced code, display math, inline math, inline code with
    null-delimited placeholders so Markdown processing cannot mangle them."""
    def fence(m):
        store["c"].append((m.group(1) or "", m.group(2)))
        return "%sC%d%s" % (NUL, len(store["c"]) - 1, NUL)
    md = re.sub(r"```[ \t]*([A-Za-z0-9_+-]*)[ \t]*\n(.*?)\n```", fence, md, flags=re.DOTALL)
    def disp(m):
        store["d"].append(m.group(1))
        return "%sD%d%s" % (NUL, len(store["d"]) - 1, NUL)
    md = re.sub(r"\$\$(.*?)\$\$", disp, md, flags=re.DOTALL)
    def inl(m):
        store["i"].append(m.group(1))
        return "%sI%d%s" % (NUL, len(store["i"]) - 1, NUL)
    md = re.sub(r"\$([^$\n]+?)\$", inl, md)
    def ic(m):
        store["k"].append(m.group(1))
        return "%sK%d%s" % (NUL, len(store["k"]) - 1, NUL)
    md = re.sub(r"`([^`]+?)`", ic, md)
    return md

def _restore(text, store):
    """Restore placeholders. Math is emitted RAW (for MathJax); code is escaped."""
    for n, (lang, code) in enumerate(store["c"]):
        cls = ' class="language-%s"' % lang if lang else ""
        repl = "<pre><code%s>%s</code></pre>" % (cls, html.escape(code))
        text = text.replace("%sC%d%s" % (NUL, n, NUL), repl)
    for n, tex in enumerate(store["d"]):
        repl = '<div class="math">\n$$%s$$\n</div>' % tex
        text = text.replace("%sD%d%s" % (NUL, n, NUL), repl)
    for n, tex in enumerate(store["i"]):
        text = text.replace("%sI%d%s" % (NUL, n, NUL), "$" + tex + "$")
    for n, code in enumerate(store["k"]):
        text = text.replace("%sK%d%s" % (NUL, n, NUL), "<code>" + html.escape(code) + "</code>")
    return text

def _inline(text):
    """Escape HTML, then apply images / links / bold / italic. Placeholders untouched."""
    text = html.escape(text, quote=False)
    text = re.sub(r"!\[([^\]]*)\]\(([^)]+)\)",
                  r'<img src="\2" alt="\1" style="max-width:100%;height:auto;display:block;'
                  r'margin:1.2rem auto;border:1px solid #d0d7de;border-radius:6px">', text)
    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2">\1</a>', text)
    text = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", text)
    text = re.sub(r"(?<!\*)\*(?!\*)(.+?)(?<!\*)\*(?!\*)", r"<em>\1</em>", text)
    return text

def _slug(headtext):
    m = re.match(r"\s*(\d+)\.", headtext)
    if m:
        return "sec-" + m.group(1)
    s = re.sub(r"%s[A-Z]\d+%s" % (NUL, NUL), "", headtext)
    s = re.sub(r"[^A-Za-z0-9]+", "-", s).strip("-").lower()
    return s or "sec"

def _render_table(header, body_rows):
    def cells(row):
        row = row.strip()
        if row.startswith("|"): row = row[1:]
        if row.endswith("|"):   row = row[:-1]
        return [c.strip() for c in row.split("|")]
    h = "".join("<th>%s</th>" % _inline(c) for c in cells(header))
    out = ['<div class="table-wrap"><table><thead><tr>%s</tr></thead><tbody>' % h]
    for r in body_rows:
        tds = "".join("<td>%s</td>" % _inline(c) for c in cells(r))
        out.append("<tr>%s</tr>" % tds)
    out.append("</tbody></table></div>")
    return "".join(out)

def _render_list(lines):
    ordered = bool(re.match(r"\s*\d+\.\s+", lines[0]))
    tag = "ol" if ordered else "ul"
    items, cur = [], None
    for ln in lines:
        m = re.match(r"\s*(?:[-*+]|\d+\.)\s+(.*)$", ln)
        if m:
            if cur is not None: items.append(cur)
            cur = m.group(1)
        else:
            cur = (cur + " " + ln.strip()) if cur is not None else ln.strip()
    if cur is not None: items.append(cur)
    has_task = any(re.match(r"\[[ xX]\]\s+", it) for it in items)
    lis = []
    for it in items:
        mt = re.match(r"\[([ xX])\]\s+(.*)$", it)
        if mt:
            chk = " checked" if mt.group(1) in ("x", "X") else ""
            lis.append('<li class="task"><input type="checkbox" disabled%s> %s</li>'
                       % (chk, _inline(mt.group(2))))
        else:
            lis.append("<li>%s</li>" % _inline(it))
    cls = ' class="tasklist"' if has_task else ""
    return "<%s%s>%s</%s>" % (tag, cls, "".join(lis), tag)

def _is_block_start(line):
    s = line.strip()
    if re.match(r"^#{1,6}\s", line): return True
    if re.match(r"^%s[DC]\d+%s$" % (NUL, NUL), s): return True
    if s.startswith(">"): return True
    if re.match(r"\s*(?:[-*+]|\d+\.)\s+", line): return True
    return False

def _parse_blocks(md):
    lines = md.split("\n")
    out, i = [], 0
    while i < len(lines):
        line = lines[i]
        if line.strip() == "":
            i += 1; continue
        m = re.match(r"^(#{1,6})\s+(.*)$", line)
        if m:
            lv = len(m.group(1)); out.append("<h%d>%s</h%d>" % (lv, _inline(m.group(2).strip()), lv))
            i += 1; continue
        if re.match(r"^%s[DC]\d+%s$" % (NUL, NUL), line.strip()):
            out.append(line.strip()); i += 1; continue
        if line.lstrip().startswith(">"):
            buf = []
            while i < len(lines) and lines[i].lstrip().startswith(">"):
                buf.append(re.sub(r"^\s*>\s?", "", lines[i])); i += 1
            out.append("<blockquote>%s</blockquote>" % _parse_blocks("\n".join(buf)))
            continue
        if "|" in line and i + 1 < len(lines) and re.match(r"^\s*\|?[\s:|-]*-[\s:|-]*\|?\s*$", lines[i+1]):
            header = line; i += 2; body = []
            while i < len(lines) and lines[i].strip() != "" and "|" in lines[i]:
                body.append(lines[i]); i += 1
            out.append(_render_table(header, body)); continue
        if re.match(r"^\s*(?:[-*+]|\d+\.)\s+", line):
            buf = []
            while i < len(lines) and lines[i].strip() != "" and \
                  (re.match(r"^\s*(?:[-*+]|\d+\.)\s+", lines[i]) or lines[i].startswith("  ")):
                buf.append(lines[i]); i += 1
            out.append(_render_list(buf)); continue
        buf = [line]; i += 1
        while i < len(lines) and lines[i].strip() != "" and not _is_block_start(lines[i]):
            buf.append(lines[i]); i += 1
        out.append("<p>%s</p>" % _inline(" ".join(b.strip() for b in buf)))
    return "\n".join(out)

def convert_doc(md, want_subtitle=False):
    """Convert a Markdown doc to HTML pieces: dict(title, subtitle, preamble, body, toc)."""
    lines = md.split("\n")
    title_md = lines[0][2:].strip() if lines[0].startswith("# ") else SITE_TITLE
    start, subtitle_md = 1, ""
    if want_subtitle:
        for j in range(1, min(8, len(lines))):
            if lines[j].startswith("### "):
                subtitle_md = lines[j][4:].strip(); start = j + 1; break
            if lines[j].startswith("## "):
                break
    rest_md = "\n".join(lines[start:])

    store = {"c": [], "d": [], "i": [], "k": []}
    rest_md = _protect(rest_md, store)

    pre_lines, sections, cur = [], [], None
    for ln in rest_md.split("\n"):
        if ln.startswith("## "):
            if cur: sections.append(cur)
            cur = {"head": ln[3:].strip(), "body": []}
        elif re.match(r"^---+\s*$", ln.strip()):
            continue
        elif cur is None:
            pre_lines.append(ln)
        else:
            cur["body"].append(ln)
    if cur: sections.append(cur)

    preamble_html = _parse_blocks("\n".join(pre_lines)) if any(l.strip() for l in pre_lines) else ""

    toc, html_sections = [], []
    for s in sections:
        slug = _slug(s["head"]); head_html = _inline(s["head"])
        toc.append((slug, head_html))
        inner = _parse_blocks("\n".join(s["body"]))
        html_sections.append(
            '<section id="%s"><h2>%s</h2>\n%s\n</section>' % (slug, head_html, inner))
    body_html = "\n".join(html_sections)

    def inline_only(t):
        st = {"c": [], "d": [], "i": [], "k": []}
        return _restore(_inline(_protect(t, st)), st)

    return {
        "title": inline_only(title_md),
        "subtitle": inline_only(subtitle_md),
        "preamble": _restore(preamble_html, store),
        "body": _restore(body_html, store),
        "toc": [(sl, _restore(tx, store)) for sl, tx in toc],
    }

# ======================================================================
#  HTML templates
# ======================================================================
MATHJAX = (
    '<script>window.MathJax={tex:{inlineMath:[[\'$\',\'$\'],[\'\\\\(\',\'\\\\)\']],'
    'displayMath:[[\'$$\',\'$$\'],[\'\\\\[\',\'\\\\]\']]},'
    'options:{skipHtmlTags:[\'script\',\'noscript\',\'style\',\'textarea\',\'pre\',\'code\']}};</script>\n'
    '<script id="MathJax-script" async '
    'src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>'
)

# Site navigation: (href-from-docs-root, label, key)
NAV = [
    ("index.html", "Home", "home"),
    ("pages/fm-phonon-defect.html", "FM self-energy &amp; defect phonons", "fm"),
    ("pages/mos2-vs-phonon-plan.html", "V<sub>S</sub>@MoS&#8322; phonons: plan", "vsplan"),
    ("pages/mos2-vs-phonon-results.html", "V<sub>S</sub>@MoS&#8322; phonons: results", "vsres"),
    ("pages/mos2-os-phonon-results.html", "O<sub>S</sub>@MoS&#8322; phonons: results", "osres"),
    ("pages/mos2-pristine-baseline.html", "Pristine baseline", "prist"),
    ("pages/supercell-bands.html", "Electronic bands", "bands"),
    ("pages/vs-defect-dos-diag.html", "V<sub>S</sub> defect DOS (diag)", "dosdiag"),
    ("pages/tmatrix.html", "T-matrix &amp; spectral function", "tmat"),
    ("pages/static-downfolding.html", "Static vs dynamic downfolding", "statdf"),
    ("pages/self-consistent-tmatrix.html", "Self-consistent T-matrix", "sctma"),
]

def _topnav(active, prefix=""):
    links = "".join(
        '<a href="%s%s"%s>%s</a>' % (
            prefix, href,
            ' style="color:#fff;text-decoration:underline"' if key == active else "",
            label)
        for href, label, key in NAV)
    return ('<nav class="topnav"><div class="inner">'
            '<span class="brand">%s</span>%s</div></nav>' % (BRAND, links))

def page_shell(title, head_html, nav_html, body_html, css_href):
    return """<!doctype html><html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<link rel="stylesheet" href="{css}">
{mathjax}
</head><body>
{nav}
{head}
<main>
{body}
</main>
<footer>Generated {date}. Static HTML; equations rendered client-side with MathJax v3.
No raw wavefunctions, cubes, binaries, or logs are published.</footer>
</body></html>""".format(title=html.escape(title), css=css_href, mathjax=MATHJAX,
                          nav=nav_html, head=head_html, body=body_html, date=GEN_DATE)

# ======================================================================
#  Page registry — to add a page: register here (and in NAV / CATALOG)
# ======================================================================
PAGES = [
    dict(
        key="fm",
        md=os.path.join(ROOT, "content", "note_fm_phonon_defect.md"),
        out=os.path.join(PAGES_DIR, "fm-phonon-defect.html"),
        html_title=SITE_TITLE + " — Fan–Migdal self-energy & defect phonon spectral function",
        subtitle=(r"Two linked derivations: the Fan&ndash;Migdal electron self-energy as the lowest-order "
                  r"electron&ndash;phonon term, built explicitly from the <em>bare</em> $G^0$ and $D^0$ "
                  r"(S-matrix, Wick, Matsubara sum, analytic continuation); then the phonon spectral function "
                  r"$B_{\mathbf q\nu}(\omega)$ of a crystal with a dilute concentration of point defects, "
                  r"via the exact defect $T$-matrix with $V(z)=\Delta\mathcal D - z\,\varepsilon$."),
        pill="Theory / derivation",
    ),
    dict(
        key="vsplan",
        md=os.path.join(ROOT, "content", "note_mos2_vs_phonon_plan.md"),
        out=os.path.join(PAGES_DIR, "mos2-vs-phonon-plan.html"),
        html_title=SITE_TITLE + " — MoS2 S-vacancy phonon spectral function: plan",
        subtitle=(r"Computational plan + TODO for $B_{\mathbf q\nu}(\omega)$ of monolayer MoS&#8322; with "
                  r"dilute S vacancies: host DFPT on a $6\times6\times1$ grid, finite-displacement "
                  r"$\Delta\Phi$ from the same $6\times6$ supercells used in the EDI electron-defect study, "
                  r"cluster T-matrix, validation gates V0&ndash;V5, and a Kestrel cost estimate."),
        pill="Plan / TODO",
    ),
    dict(
        key="vsres",
        md=os.path.join(ROOT, "content", "note_mos2_vs_phonon_results.md"),
        out=os.path.join(PAGES_DIR, "mos2-vs-phonon-results.html"),
        html_title=SITE_TITLE + " — MoS2 S-vacancy phonon spectral function: results",
        subtitle=(r"Executed plan P0&ndash;P5 on Kestrel (~13 node-h): all six gates passed. "
                  r"Vacancy-induced resonances (a&#8321; 40.9 meV, e 42.2 / 46.7 meV, band-bottom 34.5 meV) "
                  r"confirmed against the direct 321-DOF supercell diagonalization to &le;0.5 meV; "
                  r"defect-limited lifetimes $\tau_{\min}\approx37$ ps at $n_d=10^{12}\,{\rm cm^{-2}}$. "
                  r"No in-gap modes &mdash; the S vacancy scatters resonantly inside the bands."),
        pill="Numerical results",
    ),
    dict(
        key="osres",
        md=os.path.join(ROOT, "content", "note_mos2_os_phonon_results.md"),
        out=os.path.join(PAGES_DIR, "mos2-os-phonon-results.html"),
        html_title=SITE_TITLE + " — MoS2 O-substitution phonon spectral function: results",
        subtitle=(r"Same T-matrix pipeline with the frequency-dependent $V(z)=\Delta\mathcal D-z\,\varepsilon$ "
                  r"($\Delta M/M_S=-0.50$): O$_S$ produces true localized modes above the host spectrum "
                  r"&mdash; e 59.2 + a&#8321; 66.2 meV &mdash; matching the real-O-mass supercell "
                  r"diagonalization to 0.12&ndash;0.15 meV. Gentler in-band scattering than V$_S$ "
                  r"($\tau_{\min}\approx49$ ps at $10^{12}$ cm$^{-2}$). Marginal cost ~12 node-h."),
        pill="Numerical results",
    ),
    dict(
        key="prist",
        md=os.path.join(ROOT, "content", "note_mos2_pristine_baseline.md"),
        out=os.path.join(PAGES_DIR, "mos2-pristine-baseline.html"),
        html_title=SITE_TITLE + " — Pristine MoS2 phonon spectral function (baseline)",
        subtitle=(r"Null reference: with $\pi_{q\nu}=0$ the spectral function collapses to "
                  r"$\delta(\omega-\omega_{q\nu})$ ($\eta=0.05$ meV broadening only). Produced by the "
                  r"identical code path with $V=0$ &mdash; doubling as a machinery null test "
                  r"($\max|t|=0$, no resonances, $\Delta$DOS $\equiv0$). Zero DFT cost."),
        pill="Baseline",
    ),
    dict(
        key="bands",
        md=os.path.join(ROOT, "content", "note_supercell_bands.md"),
        out=os.path.join(PAGES_DIR, "supercell-bands.html"),
        html_title=SITE_TITLE + " — Electronic band structure of the three MoS2 supercells",
        subtitle=(r"Companion electronic-structure comparison along $\Gamma$&ndash;M&ndash;K&ndash;$\Gamma$ "
                  r"(coarse 26-$k$ path). V$_S$ shows the classic empty $e$ defect doublet deep in the gap "
                  r"(+1.15 eV, width &lt;10 meV); O$_S$ is isovalent &mdash; no gap state; pristine ~1.66 eV gap. "
                  r"The electronic mirror of the phonon story."),
        pill="Electronic structure",
    ),
    dict(
        key="dosdiag",
        md=os.path.join(ROOT, "content", "note_vs_defect_dos_diag.md"),
        out=os.path.join(PAGES_DIR, "vs-defect-dos-diag.html"),
        html_title=SITE_TITLE + " — V_S defect DOS by single-defect Hamiltonian diagonalization",
        subtitle=(r"Literature single-defect Hamiltonian $H=\mathrm{diag}(\varepsilon_{n\mathbf k})+M\,\mathrm{Ry}/N_{\mathbf k}$ "
                  r"from EDI (direct mode). Two fixes (found vs an independent EDT code): the "
                  r"$\mathrm{Ry}/N_{\mathbf k}$ normalization, and a stale-bra bug that made $M$ non-Hermitian "
                  r"and split the $e$-doublet. Corrected V$_S$: $a_1$ + a DEGENERATE $e$-doublet at "
                  r"$+1.196$ eV (EDT: $+1.205$); $C_3$ restored (143/144). Earlier numbers retracted."),
        pill="Method / electronic structure",
    ),
    dict(
        key="tmat",
        md=os.path.join(ROOT, "content", "note_tmatrix.md"),
        out=os.path.join(PAGES_DIR, "tmatrix.html"),
        html_title=SITE_TITLE + " — Electron-defect T-matrix, self-energy & spectral function (MoS2 V_S, O_S)",
        subtitle=(r"Beyond-Born electron&ndash;defect scattering by explicit-summation downfolding + Wannier "
                  r"interpolation (Anvil EDT path, ported to Kestrel, validated to the meV): the on-shell "
                  r"diagonal $T(\mathbf k,\mathbf k;\omega)$, defect self-energy $n_d T$, and spectral function "
                  r"$A(\mathbf k,\omega)$ along $\Gamma$&ndash;M&ndash;K&ndash;$\Gamma$. O$_S$ scatters the VBM "
                  r"~3&ndash;4&times; harder than V$_S$; V$_S$ carries the stronger rest-space dressing."),
        pill="Beyond-Born / T-matrix",
    ),
    dict(
        key="statdf",
        md=os.path.join(ROOT, "content", "note_static_downfolding.md"),
        out=os.path.join(PAGES_DIR, "static-downfolding.html"),
        html_title=SITE_TITLE + " — Static vs dynamic rest-space in the downfolded T-matrix",
        subtitle=(r"Derivation: the rest dressing $T^R=[1-gG^R(\omega_0)]^{-1}g$, its static limit at a single "
                  r"reference $\omega_0$, the exact frequency-dependent form $T^R(\omega)$, and the proof that the "
                  r"two-layer split $T=[1-T^R G^A]^{-1}T^R$ is exact iff $T^R$ runs with $\omega$. Static error "
                  r"$\delta T^R\!\approx\!(\omega-\omega_0)\,T^R(\partial_\omega G^R)T^R$ &rArr; valid only when "
                  r"$W_A/\Delta\ll1$, so it degrades as the active space grows."),
        pill="Theory / downfolding",
    ),
    dict(
        key="sctma",
        md=os.path.join(ROOT, "content", "note_self_consistent_tmatrix.md"),
        out=os.path.join(PAGES_DIR, "self-consistent-tmatrix.html"),
        html_title=SITE_TITLE + " — Born, self-consistent Born, and the self-consistent T-matrix",
        subtitle=(r"Derivation of the disorder self-energy hierarchy for dilute defects — 1st Born &rarr; SCBA "
                  r"&rarr; $T$-matrix &rarr; self-consistent $T$-matrix (SCTMA), $\Sigma=n_d\,v[1-\mathcal G v]^{-1}$ "
                  r"with $\mathcal G=[(G^A)^{-1}-\Sigma]^{-1}$. Shown to be the dilute limit of CPA, mapped onto the "
                  r"downfolded Koster&ndash;Slater pipeline, with what it can fix (the V$_S$ $a_1$ under-binding) "
                  r"and cannot (inter-defect interference, full DFT self-consistency)."),
        pill="Theory / beyond-Born",
    ),
]

# Landing-page catalog rows: (item, type, date, badge_class, badge_label, summary, link_html)
CATALOG = [
    (r"Fan&ndash;Migdal self-energy &amp; defect phonon spectral function", "Theory", GEN_DATE, "ok", "Complete",
     r"Lowest-order e&ndash;ph self-energy from bare $G^0$ and $D^0$: S-matrix expansion, Wick contractions "
     r"(first order and tadpoles vanish at the relaxed geometry), the boson Matsubara sum done explicitly "
     r"($n_B(\xi-i\omega_j)=-n_F(\xi)$), retarded $\Sigma^{\rm FM}$ with absorption/emission golden-rule structure. "
     r"Then phonons with defects: Lifshitz perturbation $V(z)=\Delta\mathcal D-z\varepsilon$, exact small-block "
     r"$T$-matrix (phonon Koster&ndash;Slater), dilute average $\pi_{q\nu}=c\,t_{q\nu}$, spectral function "
     r"$B_{q\nu}(\omega)$ with Tamura/Born and resonant-mode limits, plus the non-adiabatic e&ndash;ph bubble "
     r"with defect-dressed electrons.",
     '<a href="pages/fm-phonon-defect.html">Open derivation &rarr;</a>'),
    (r"MoS&#8322; S-vacancy phonon spectral function &mdash; results", "Result", GEN_DATE, "prod", "Complete",
     r"Full pipeline executed (~13 node-h): host DFPT + 124-displacement defect FD (P3m1 preserved), "
     r"32-atom cluster T-matrix. Gates V0&ndash;V5 all passed (Born&rarr;Tamura scaling 10.5%&rarr;0.40% as "
     r"$\epsilon$ 0.06&rarr;0.01; $\int B\,d\omega=0.9985$; resonances match supercell diagonalization "
     r"&le;0.5 meV). V$_S$ gives in-band resonances (a&#8321; 40.9, e 42.2/46.7, 34.5 meV), no gap states; "
     r"$\Gamma_{\max}\approx9\times10^{-3}$ meV at 43.4/34.3 meV &rArr; $\tau_{\min}\approx37$ ps at "
     r"$n_d=10^{12}$ cm$^{-2}$.",
     '<a href="pages/mos2-vs-phonon-results.html">Open results &rarr;</a>'),
    (r"MoS&#8322; O-substitution (O$_S$) phonon spectral function &mdash; results", "Result", GEN_DATE, "prod", "Complete",
     r"First exercise of the frequency-dependent $V(z)=\Delta\mathcal D-z\varepsilon$ (mass term "
     r"$\varepsilon=-0.50$): three local modes ejected above the host spectrum &mdash; e doublet 59.2 meV + "
     r"a&#8321; 66.2 meV (PR 0.004, cluster weight 0.999) &mdash; T-matrix poles vs real-O-mass supercell "
     r"diagonalization agree to 0.12&ndash;0.15 meV. In-band scattering gentler than V$_S$: "
     r"$\bar\Gamma=1.1\times10^{-3}$ meV, $\tau_{\min}\approx49$ ps at $10^{12}$ cm$^{-2}$. "
     r"Reused host DFPT + pristine FCs + cached $\rho_{ab}(\lambda)$; ~12 node-h marginal.",
     '<a href="pages/mos2-os-phonon-results.html">Open results &rarr;</a>'),
    (r"Pristine MoS&#8322; phonon spectral function (baseline)", "Result", GEN_DATE, "ok", "Complete",
     r"Null reference for the defect runs: $B_{q\nu}(\omega)\to\delta(\omega-\omega_{q\nu})$ at "
     r"$\pi_{q\nu}=0$, same path/mesh/$\eta$ and code path with $V=0$ (null test: $\max|t|=0$, "
     r"$\int B=0.9987$, no resonances, $\Delta$DOS$\,\equiv0$). Includes the three-way "
     r"pristine / V$_S$ / O$_S$ comparison table. Zero DFT cost (~5 min post-processing).",
     '<a href="pages/mos2-pristine-baseline.html">Open baseline &rarr;</a>'),
    (r"Electronic band structure of the supercells (6×6 &amp; 9×9)", "Result", "2026-06-17", "ok", "Complete",
     r"Companion electronic comparison ($\Gamma$&ndash;M&ndash;K&ndash;$\Gamma$, sparse $k$, standard "
     r"whole-node): V$_S$ shows the empty $e$ defect doublet deep in the gap; O$_S$ isovalent (no in-gap "
     r"state, gap 1.58 vs pristine 1.66 eV). Plus a 9×9 (242-atom) V$_S$ size-convergence run: the $e$ "
     r"doublet dispersion collapses ~10× (9.1&rarr;0.9 meV) and drops to +1.06 eV &mdash; direct evidence "
     r"the defect states reach the isolated limit, validating the dilute assumption of the phonon work.",
     '<a href="pages/supercell-bands.html">Open band structure &rarr;</a>'),
    (r"Defect DOS by single-defect Hamiltonian diagonalization (MoS$_2$ V$_S$)", "Result", "2026-06-20", "prod", "Complete",
     r"Diagonalize $H=\mathrm{diag}(\varepsilon_{n\mathbf k})+M\,\mathrm{Ry}/N_{\mathbf k}$ from EDI direct mode. "
     r"Two corrections, both found by cross-checking an independent EDT code: (1) the $\mathrm{Ry}/N_{\mathbf k}$ "
     r"normalization; (2) a stale-bra bug in `ed_direct_from_files` that made the same-$\mathbf k$ block "
     r"**non-Hermitian** ($\lVert M-M^\dagger\rVert/\lVert M\rVert\!\approx\!2$), breaking $C_3$ and splitting the "
     r"$e$-doublet. After the fix: Hermitian to $10^{-13}$, $C_3$ triplets 143/144 (= EDT), and V$_S$ shows an "
     r"$a_1$ + a **degenerate $e$-doublet at $+1.196$ eV** (EDT $+1.205$; DFT 1.06&ndash;1.15) &mdash; the correct "
     r"defect structure. Earlier tables retracted.",
     '<a href="pages/vs-defect-dos-diag.html">Open defect DOS &rarr;</a>'),
    (r"Electron&ndash;defect $T$-matrix, self-energy &amp; spectral function (MoS$_2$ V$_S$, O$_S$)", "Result", "2026-06-24", "prod", "Complete",
     r"Beyond-Born: explicit-summation downfolding (11 active Wannier bands + exact static rest "
     r"$\Sigma(\omega_0)$, 60$\to$11) + Wannier interpolation, ported from the independent Anvil EDT code and "
     r"validated to the meV (V$_S$ $+1.21$, O$_S$ $+0.73$ eV). On-shell diagonal $T(\mathbf k,\mathbf k;\omega)$, "
     r"self-energy $n_d T$, and $A(\mathbf k,\omega)$ along $\Gamma$&ndash;M&ndash;K&ndash;$\Gamma$: O$_S$ "
     r"scatters the VBM ~3.7&times; harder than V$_S$ (Im $T$ at K $-0.24$ vs $-0.065$ Ry; defect broadening "
     r"~90 vs ~31 meV at $n_d=2.78\%$); rest-space dressing strong for V$_S$, mild for O$_S$.",
     '<a href="pages/tmatrix.html">Open T-matrix &rarr;</a>'),
    (r"Static vs dynamic rest-space in the downfolded $T$-matrix", "Theory", "2026-06-24", "ok", "Complete",
     r"Derivation note. The two-layer $T$-matrix $T=[1-T^R G^A(\omega)]^{-1}T^R$ with rest dressing "
     r"$T^R=[1-gG^R(\omega_0)]^{-1}g$: shown identical to the Feshbach self-energy $\tilde V=W_{PP}+\Sigma_P$, "
     r"and the split proven exact iff $T^R$ is taken at the running $\omega$ &mdash; so freezing $T^R(\omega_0)$ "
     r"is the only frequency approximation in Layer 1. Static error "
     r"$\delta T^R\approx(\omega-\omega_0)T^R(\partial_\omega G^R)T^R$, equivalently a $Z=(1+\beta)^{-1}<1$ "
     r"renormalization, valid only for $W_A/\Delta\ll1$ (degrades as the active space grows). Pole positions are "
     r"weakly affected (V$_S$: static $+1.209$ vs dynamic $+1.192$ eV, $Z\approx0.99$); spectral weight and "
     r"satellites are not. Includes the exact dynamic $T^R(\omega)$ form (cheap: one GEMM/$\omega$) and a "
     r"comparison plan. Derivation only &mdash; not yet run.",
     '<a href="pages/static-downfolding.html">Open derivation &rarr;</a>'),
    (r"Born, self-consistent Born &amp; the self-consistent $T$-matrix (SCTMA)", "Theory", "2026-06-24", "ok", "Complete",
     r"Derivation of the dilute-defect self-energy hierarchy: 1st Born ($n_d vG_0v$, = golden-rule EDI) &rarr; SCBA "
     r"($n_d v\bar G v$, rainbow diagrams, no bound states) &rarr; $T$-matrix ($n_d v[1-G_0v]^{-1}$, our current "
     r"beyond-Born result) &rarr; SCTMA ($n_d v[1-\mathcal G v]^{-1}$ with $\mathcal G=[(G^A)^{-1}-\Sigma]^{-1}$, "
     r"self-consistent). Proven to be the $O(n_d)$ limit of CPA; mapped onto the Koster&ndash;Slater pipeline "
     r"(one extra self-consistency loop, $G^A\!\to\!\mathcal G$). Numerically tested on MoS$_2$ V$_S$: at the "
     r"physical $n_d{=}2.78\%$ SCTMA moves the $a_1$ by only $-9$ meV (does NOT lift it toward DFT $+0.14$); its "
     r"actual effect is a disorder broadening of the deep $e$ level. The band-edge $a_1$ needs DFT self-consistency.",
     '<a href="pages/self-consistent-tmatrix.html">Open derivation &rarr;</a>'),
    (r"MoS&#8322; S-vacancy phonon spectral function &mdash; plan", "Plan", "2026-06-09", "ok", "Executed",
     r"Workflow P0&ndash;P6 reusing the EDI structural model ($a=3.185$ &#8491;, $6\times6$ supercells, "
     r"NC stringent pseudos, $E_{\rm cut}=100$ Ry, 2D Coulomb cutoff): host DFPT ($6\times6\times1$ q-grid, "
     r"2D LO&ndash;TO), defect-supercell relaxation + phonopy finite displacements, $\Delta\Phi$ cluster "
     r"truncation with ASR, dense-grid host resolvent, T-matrix $\to\pi_{q\nu}=c\,t_{q\nu}$ at "
     r"$n_d=10^{12}\,{\rm cm^{-2}}$ ($c\approx8.8\times10^{-4}$). Six validation gates (incl. Tamura Born "
     r"check and supercell-diagonalization fingerprint); est. 15&ndash;30 node-h.",
     '<a href="pages/mos2-vs-phonon-plan.html">Open plan &amp; TODO &rarr;</a>'),
]

def build_page(p):
    with open(p["md"], encoding="utf-8") as f:
        md = f.read()
    r = convert_doc(md, want_subtitle=False)
    toc_links = "".join('<a href="#%s">%s</a>' % (sl, tx) for sl, tx in r["toc"])
    header = ('<header><div class="header-inner"><h1>{t}</h1>'
              '<p class="subtitle">{s}</p>'
              '<div class="meta"><span class="pill">{p}</span>'
              '<span class="pill">{n} sections</span>'
              '<span class="pill">MathJax v3</span>'
              '<span class="pill">Generated {d}</span></div></div></header>'
             ).format(t=r["title"], s=p["subtitle"], p=p["pill"], n=len(r["toc"]), d=GEN_DATE)
    toc_section = ('<section id="contents"><h2>Contents</h2>'
                   '<div class="toc">%s</div></section>' % toc_links)
    body = toc_section + "\n" + r["preamble"] + "\n" + r["body"]
    out = page_shell(p["html_title"], header, _topnav(p["key"], prefix="../"),
                     body, "../assets/style.css")
    with open(p["out"], "w", encoding="utf-8") as f:
        f.write(out)

def build_index():
    rows = ""
    for item, typ, date, bc, bl, summ, link in CATALOG:
        planned = ' class="planned"' if bc == "plan" else ""
        rows += ("<tr%s><td><strong>%s</strong></td><td>%s</td><td>%s</td>"
                 "<td><span class=\"badge %s\">%s</span></td><td>%s</td><td>%s</td></tr>\n"
                 % (planned, item, typ, date, bc, bl, summ, link))
    catalog = (
        '<section id="catalog"><h2>Catalog</h2>'
        '<p>One row per piece of work: what it is, when, status, the key result, and a link.</p>'
        '<div class="table-wrap"><table><thead><tr>'
        '<th>Item</th><th>Type</th><th>Date</th><th>Status</th><th>Key result / summary</th><th>Link</th>'
        '</tr></thead><tbody>%s</tbody></table></div>'
        '<p class="small">Legend: '
        '<span class="badge ok">Complete</span> done &nbsp; '
        '<span class="badge plan">Planned</span> not yet run &nbsp; '
        '<span class="badge prod">Production</span> headline result.</p></section>'
        % rows)

    warnings = (
        '<section id="not-published" class="warning"><h2>Warnings: Files Not Published</h2>'
        '<p>This GitHub Pages artifact contains only the static report under <code>docs/</code>. '
        'It intentionally excludes raw or large research data so the published site stays small '
        'and contains no credentials.</p>'
        '<div class="table-wrap"><table><thead><tr><th>Class</th><th>Policy</th></tr></thead><tbody>'
        '<tr><td>Wavefunctions, QE <code>*.save/</code></td><td>Excluded (large binary).</td></tr>'
        '<tr><td>Cube / volumetric grids <code>*.cube</code></td><td>Excluded (often &gt;100&nbsp;MB).</td></tr>'
        '<tr><td>Numerical arrays <code>*.npy/*.npz/*.bin/*.h5/*.dat</code></td><td>Excluded; summarized in tables only.</td></tr>'
        '<tr><td>Scheduler logs <code>*.out/*.err/*.log</code></td><td>Excluded.</td></tr>'
        '<tr><td>Local/private config <code>.claude/</code>, keys/tokens</td><td>Excluded.</td></tr>'
        '</tbody></table></div></section>')

    header = ('<header><div class="header-inner"><h1>{t}</h1>'
              '<p class="subtitle">Phonon&ndash;defect interaction project &mdash; theory notes, '
              'derivations, and (forthcoming) numerical results.</p>'
              '<div class="meta">'
              '<span class="pill">Generated {d}</span>'
              '<span class="pill">GitHub Pages: docs/</span>'
              '<span class="pill">Branch: main</span>'
              '<span class="pill">MathJax v3</span></div></div></header>'
             ).format(t=SITE_TITLE, d=GEN_DATE)

    body = catalog + warnings
    out = page_shell(SITE_TITLE, header, _topnav("home"), body, "assets/style.css")
    with open(os.path.join(DOCS, "index.html"), "w", encoding="utf-8") as f:
        f.write(out)

# ======================================================================
#  main + self-checks
# ======================================================================
def main():
    os.makedirs(PAGES_DIR, exist_ok=True)
    for p in PAGES:
        build_page(p)
    build_index()

    def check(txt):
        problems = []
        if NUL in txt: problems.append("UNRESTORED placeholder (\\x00) present!")
        leftover = re.findall(r"%s[A-Z]\d+%s" % (NUL, NUL), txt)
        if leftover: problems.append("leftover tokens: %s" % leftover[:5])
        return problems

    print("=== build_site.py ===")
    files = [(os.path.basename(p["out"]), p["out"]) for p in PAGES]
    files.append(("index.html", os.path.join(DOCS, "index.html")))
    for nm, path in files:
        txt = open(path, encoding="utf-8").read()
        print("%-28s: %d bytes, %d sections, %d display-eq, %d tables" % (
            nm, len(txt), txt.count("<section id="),
            txt.count('class="math"'), txt.count("<table>")))
        p = check(txt)
        print("  [%s] %s" % (nm, "OK" if not p else " ; ".join(p)))

if __name__ == "__main__":
    main()
