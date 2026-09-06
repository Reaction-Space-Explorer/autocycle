# How every committed figure was made. `make figures` rebuilds the ones that need only
# what is in this repository. PNGs need libcairo; on macOS with homebrew add
# CAIRO_LIB=/opt/homebrew/lib, since macOS strips DYLD_* from the shell make spawns.
# The style swatches need obabel on PATH.

PY ?= python -m autocycle.cli
# obabel draws small molecules more legibly (O = CH2, not a bare =O), but it cannot
# label an R-group, so the two cycles carrying the [CoA] stub are drawn by RDKit.
CANON_OB := formose_core krebs_tca
CANON_RD := acetyl_coa_sol0 malyl_coa_arm
STYLES := paper annotated rich
CAIRO_LIB ?=

.PHONY: figures svg png gallery panel test lint

figures: svg png

svg:
	@for c in $(CANON_OB); do \
	  $(PY) draw examples/canonical/$$c.yaml --style annotated --backend obabel \
	    -o examples/canonical/$$c.svg; \
	done
	@for c in $(CANON_RD); do \
	  $(PY) draw examples/canonical/$$c.yaml --style annotated \
	    -o examples/canonical/$$c.svg; \
	done
	@for s in $(STYLES); do \
	  $(PY) draw examples/formose_gain.yaml --style $$s --backend obabel \
	    -o examples/styles/cycle_$$s.svg; \
	  $(PY) route examples/ribose_route.yaml --style $$s --backend obabel \
	    -o examples/styles/route_$$s.svg; \
	done
	$(PY) linear examples/canonical/formose_core.yaml --backend obabel \
	  -o examples/styles/linear_paper.svg
	$(PY) from-crs examples/crs/example-01.crs --style annotated \
	  -o examples/crs/example-01.svg
	$(PY) bench-routes examples/traced_sample \
	  --seeds examples/traced_sample/products.tsv \
	  --rels examples/traced_sample/rels.tsv \
	  --sample 3 --sample-dir examples/routes

png:
	@DYLD_FALLBACK_LIBRARY_PATH=$(CAIRO_LIB) python -c "import pathlib; from autocycle.export import write; \
	[write(p.read_text(), p.with_suffix('.png')) for p in \
	 sorted(pathlib.Path('examples').rglob('*.svg'))]"

# The search figures come from the full glucose corpus of Arya et al. 2022, which is not
# redistributed here. Point CORPUS at a directory of Cypher result CSVs to rebuild them.
gallery:
	@test -n "$(CORPUS)" || (echo "set CORPUS=/path/to/cycle/csvs"; false)
	$(PY) bench $(CORPUS) --sample 5 --sample-dir examples/gallery --style annotated

panel:
	@test -n "$(CORPUS)" || (echo "set CORPUS=/path/to/cycle/csvs"; false)
	$(PY) panel $(CORPUS) --sample 2 -o examples/figure_panel.svg

test:
	pytest -q

lint:
	ruff check src tests
