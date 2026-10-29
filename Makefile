CONFIG ?= config/run.yaml
CORES ?= 1

.PHONY: lint validate report-test dry-run run

lint:
	python3 -m py_compile pipeline/scripts/*.py pipeline/validation/*.py

validate:
	python3 pipeline/validation/validate_config.py --config $(CONFIG)

report-test:
	python3 tests/test_report.py

dry-run:
	bash pipeline/run_snakemake.sh --config $(CONFIG) --cores $(CORES) --dry-run

run:
	bash pipeline/run_snakemake.sh --config $(CONFIG) --cores $(CORES)
