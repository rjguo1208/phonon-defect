# phonon-defect

Phonon–defect interaction project.

The GitHub Pages site lives under `docs/` (published from the `main` branch,
`/docs` folder) — same layout as the `claude-sternheimer` repo. Markdown
sources live under `content/` and are rendered into `docs/` by
`python3 tools/build_site.py` (stdlib-only, math-protected Markdown→MathJax
HTML). To add a page: drop a note in `content/`, register it in the
`PAGES`/`NAV`/`CATALOG` lists of `tools/build_site.py`, and re-run.

Publish policy: only the static report under `docs/` is published. Raw or
large research data (wavefunctions, `*.save/`, cubes, binary arrays,
scheduler logs) are never committed — see `.gitignore`.
