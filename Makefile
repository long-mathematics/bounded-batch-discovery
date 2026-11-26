PYTHON ?= python3
LATEXMK ?= latexmk
PAPER := papers/bounded_batch_discovery
CHECK_ENV := OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

.PHONY: all pdf check check-finite check-revision validate clean

all: validate

# latexmk tracks source changes and repeats pdflatex until references settle.
pdf:
	$(LATEXMK) -cd -pdf -interaction=nonstopmode -halt-on-error -file-line-error $(PAPER).tex

check: check-finite check-revision

check-finite:
	$(CHECK_ENV) $(PYTHON) scripts/verify_finite_models.py --output reports/verification_results.json

check-revision:
	$(CHECK_ENV) $(PYTHON) scripts/verify_revision.py --output reports/revision_verification_results.json

# Depend on both steps so validation sees the current PDF and verification scripts.
validate: pdf check
	$(PYTHON) scripts/validate_manuscript.py --output reports/manuscript_validation.json

# Preserve the committed paper PDF and JSON reports. Use -c, never latexmk -C.
clean:
	$(LATEXMK) -cd -c $(PAPER).tex
	rm -f $(PAPER).synctex.gz
