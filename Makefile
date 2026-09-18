.PHONY: test validate syntax predict-test predict-demo

PYTHON ?= python3

test:
	$(PYTHON) -m pytest tests -q

validate:
	$(PYTHON) tools/validate_public_artifacts.py

syntax:
	$(PYTHON) -c 'import ast; from pathlib import Path; files=list(Path("src").rglob("*.py"))+list(Path("scripts").rglob("*.py"))+list(Path("tools").rglob("*.py"))+list(Path("tests").rglob("*.py"))+list(Path("predict").rglob("*.py")); [ast.parse(p.read_text(encoding="utf-8"), filename=str(p)) for p in files]; print(f"parsed {len(files)} Python files")'

# Requires: pip install -r predict/requirements-inference.txt
predict-test:
	$(PYTHON) -m pytest tests/test_predict.py tests/test_render_variant_report.py -q

predict-demo:
	rm -rf artifacts/predict_demo
	$(PYTHON) predict/predict_variants.py --input predict/examples/input.csv --bundle predict/bundle --output artifacts/predict_demo
	$(PYTHON) predict/render_variant_report.py --predictions artifacts/predict_demo/predictions.jsonl --out artifacts/predict_demo/report.md
	@echo "Report written to artifacts/predict_demo/report.md"
