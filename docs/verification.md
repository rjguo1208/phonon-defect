# Verification record: repository restart (2026-10-01)

What was run for the first version of the repurposed repository, and what it shows.
Environment: Banff login node (tiny toy problems only), Python 3.12.12, NumPy 2.4.0,
SciPy 1.16.3, Node 22.23.3, KaTeX 0.18.7.

## Unit tests

`python -m pytest -q tests` → 8 passed in 7.9 s (41 x 41 lattice):
Lanczos orthonormality; defect-independent coupling block; local mode exact with
12 chain blocks (1e-9); continued fraction = mode sum (1e-10); pristine, light and
heavy defect reproduce the exact cluster Green function at eta = 0.3 (1e-6);
damped companion matrix reproduces [(w + i G)^2 - C]^-1 (1e-10) with all
eigenvalues in the lower half plane.

## Toy validation (`python examples/square_lattice_toy.py`)

121 x 121 lattice, cluster of 13 sites, bath of 14628 sites, block width 8,
160 blocks built once in 15.5 s and shared by both defects.

| Defect | Blocks | Top mode (upfolded / exact) | Max rel. error of G_00 |
| --- | --- | --- | --- |
| light (M'=0.3, K'=1.3) | 10 | 4.348826 / 4.348826 | 6.0e-1 |
| light | 40 | 4.348826 / 4.348826 | 4.7e-2 |
| light | 160 | 4.348826 / 4.348826 | 4.9e-6 |
| heavy (M'=8, K'=0.6) | 10 | 2.812116 / 2.828119 | 4.8e-1 |
| heavy | 40 | 2.826994 / 2.828119 | 5.0e-2 |
| heavy | 160 | 2.828119 / 2.828119 | 4.9e-6 |

Continued fraction vs mode sum at 40 blocks: 1.3e-15 (light), 4.7e-14 (heavy).

## Plasmon–phonon upfolding checks (`reference/plasmon_phonon_upfolding/`)

`check_upfold.py`: S14->S15 5.62e-16; Eq. 9 eigen-expansion 1.46e-13; pairing
7.14e-14 with max Im(lambda) = -0.302; S59 vs S58 with static conjugation 8.33e-16
(literal conjugation 3.80e-2 in this weak-coupling model); multipole with
pole-resolved couplings 7.17e-16 (pole-independent 1.85e-2).

`check_loss.py` (strong coupling): S75 vs RPA closed form 4.53e-15; literal S74
1.64; loss-function integral 34.8712 (literal) vs 33.2308 (S75); peak 39.91 vs
40.33.

## Toy figures (`npm run figures`, added 2026-10-01)

`scripts/toy_figures.py` reuses the toy lattice and chain and writes
`site/results/toy-local-dos.{svg,pdf,png}`, `site/results/toy-convergence.{svg,pdf,png}`
and the plotted numbers as CSV in `site/data/`. Full scan of the maximum relative
error of G_00 (light / heavy): m = 5: 9.3e-1 / 9.5e-1; 20: 2.3e-1 / 2.1e-1;
60: 1.1e-2 / 1.2e-2; 100: 6.1e-4 / 6.4e-4; 140: 2.4e-5 / 2.5e-5; 160: 4.9e-6 / 4.9e-6.
Local-mode frequency error of the light defect: m = 1: 1.2e-6; 2: 2.0e-8; 3: 3.5e-10;
4: 6.0e-12; 5: 1.1e-13; 6: 1.8e-15; 7-10: 8.9e-16 (double-precision floor).

Colours were checked with the data-viz palette validator against the white page
surface: the ordinal blue ramp `#86b6ef,#2a78d6,#104281` (m = 10, 40, 160) passes all
ramp checks (light end 2.11:1); the categorical pair `#2a78d6,#eb6834` passes all
categorical checks (worst CVD Delta E 24.7). The SVGs have a viewBox, no raster
images and text converted to paths. The PNG renders were inspected for label
collisions and clipping.

## Website

`npm run build` → 5 pages, 195 LaTeX expressions rendered with KaTeX in strict
mode. `npm run check` → "OK: 5 English pages; 195 LaTeX expressions with MathML;
0 scientific plots; page navigation, local links, anchors and fonts." No browser
check was run.

## Not verified

Nothing first-principles has been run in this version; every gate on the Plan
page except U0 is open.
