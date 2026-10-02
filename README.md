# phonon-defect

Website: <https://rjguo1208.github.io/phonon-defect/> (published by the manual
"Deploy GitHub Pages" workflow)

Electrons scattered by point defects and phonons at once. The defect potential is
kept to all orders (electron–defect T-matrix) and the electron–phonon coupling to
lowest order; the goal is the cross terms of order `n_d g^2` that Matthiessen's
rule leaves out: resonance-enhanced phonon emission and absorption, scattering from
defect-localized vibrations, the defect-modified electron–phonon vertex, and
defect–phonon interference.

The method borrows the upfolding of Lihm and Park
([arXiv:2409.07393](https://arxiv.org/abs/2409.07393)): a frequency-dependent
self-energy is replaced by static couplings to auxiliary modes, so one static
matrix replaces a frequency loop and electrons couple through static vertices. For a
point defect the auxiliary modes are the sites of a block Lanczos chain of the host
phonon bath, attached to the defect cluster.

## Status

- Method: designed (`src/method.html`).
- Upfolding core `phdef`: implemented; 8 unit tests pass; on a 121 x 121 toy
  lattice a local mode is exact with 10 chain blocks and an in-band resonance
  converges to a relative error of 5e-6 with 160 blocks, with one chain shared by
  all defects (`src/toy.html`).
- Plasmon–phonon upfolding (Lihm and Park): re-derived and checked numerically, with two notes on the
  printed equations (`src/plasmon-phonon.html`).
- First-principles MoS2 (O_S, V_S): planned, not run (`src/plan.html`).

## Layout

| Path | Content |
| --- | --- |
| `phdef/` | block Lanczos chain, upfolded matrix, cluster Green function (mode sum and continued fraction), damped companion matrix, toy models |
| `examples/square_lattice_toy.py` | toy validation (about one minute) |
| `tests/` | unit tests: `python -m pytest tests` |
| `reference/plasmon_phonon_upfolding/` | numerical checks of the plasmon–phonon upfolding identities (Lihm and Park) |
| `src/`, `site/` | website source and generated website |
| `scripts/` | website build (`build.mjs`) and checks (`check_site.py`) |
| `docs/verification.md` | what was checked for the current version |

## Website

```sh
npm ci --ignore-scripts
npm run figures  # toy-validation figures and CSV into site/results/ and site/data/ (about one minute)
npm run schematic  # toy-model schematic from scripts/toy_schematic.tex (TikZ; needs pdflatex and PyMuPDF)
npm run build    # KaTeX renders all LaTeX into site/ (HTML + MathML, no client JavaScript)
npm run check    # rendered math, navigation, links, anchors, fonts, English-only text
```

Commit the regenerated `site/` with every source change; CI fails if it is stale.
The website is English only; see `AGENTS.md`.

## Earlier content

Until 2026-07-01 this repository held the phonon-side T-matrix of MoS2 with V_S and
O_S and notes on the electron–defect T-matrix. That work is preserved at the tag
[`archive/2026-07-01`](https://github.com/rjguo1208/phonon-defect/tree/archive/2026-07-01).

Raw research data (wavefunctions, `*.save/`, cubes, binary arrays, scheduler logs)
are never committed; see `.gitignore`.
