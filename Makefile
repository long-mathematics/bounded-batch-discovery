PYTHON ?= python3
LATEXMK ?= latexmk
PAPERS := papers/bounded_batch_discovery papers/optimal_berry_esseen papers/quantitative_extensions
CHECK_ENV := OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

.PHONY: all pdf check check-finite check-revision check-optimal check-extensions validate clean

all: validate

# latexmk tracks source changes and repeats pdflatex until references settle.
pdf:
	$(LATEXMK) -cd -pdf -interaction=nonstopmode -halt-on-error -file-line-error papers/bounded_batch_discovery.tex
	$(LATEXMK) -cd -pdf -interaction=nonstopmode -halt-on-error -file-line-error papers/optimal_berry_esseen.tex
	$(LATEXMK) -cd -pdf -interaction=nonstopmode -halt-on-error -file-line-error papers/quantitative_extensions.tex

check: check-finite check-revision check-optimal check-extensions

check-finite:
	$(CHECK_ENV) $(PYTHON) scripts/verify_finite_models.py --output reports/verification_results.json

check-revision:
	$(CHECK_ENV) $(PYTHON) scripts/verify_revision.py --output reports/revision_verification_results.json

check-optimal:
	$(CHECK_ENV) $(PYTHON) scripts/verify_optimal_rate.py --output reports/optimal_rate_verification.json

check-extensions:
	$(CHECK_ENV) $(PYTHON) scripts/verify_quantitative_extensions.py --output reports/quantitative_extensions_verification.json

# Depend on both steps so validation sees the current PDF and verification scripts.
validate: pdf check
	$(PYTHON) scripts/validate_manuscript.py --output reports/manuscript_validation.json
	$(PYTHON) scripts/validate_manuscript.py --source papers/optimal_berry_esseen.tex --pdf papers/optimal_berry_esseen.pdf --log papers/optimal_berry_esseen.log --header-text "Optimal bounds for bounded-batch discovery" --output reports/optimal_manuscript_validation.json
	$(PYTHON) scripts/validate_manuscript.py --source papers/quantitative_extensions.tex --pdf papers/quantitative_extensions.pdf --log papers/quantitative_extensions.log --header-text "Further bounds for batch discovery" --output reports/quantitative_extensions_validation.json

# Preserve the committed paper PDF and JSON reports. Use -c, never latexmk -C.
clean:
	$(LATEXMK) -cd -c papers/bounded_batch_discovery.tex
	$(LATEXMK) -cd -c papers/optimal_berry_esseen.tex
	$(LATEXMK) -cd -c papers/quantitative_extensions.tex
	rm -f $(addsuffix .synctex.gz,$(PAPERS))
