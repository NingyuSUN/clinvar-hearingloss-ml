.PHONY: test validate syntax

PYTHON ?= python3

test:
	$(PYTHON) -m pytest tests -q

validate:
	$(PYTHON) tools/validate_public_artifacts.py

syntax:
	$(PYTHON) -c 'import ast; from pathlib import Path; files=list(Path("src").rglob("*.py"))+list(Path("scripts").rglob("*.py"))+list(Path("tools").rglob("*.py"))+list(Path("tests").rglob("*.py")); [ast.parse(p.read_text(encoding="utf-8"), filename=str(p)) for p in files]; print(f"parsed {len(files)} Python files")'
