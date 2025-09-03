UV     := uv
PIP    := pip
UVX    := $(UV) run
PYTHON := $(UVX) python


.PHONY: build format install publish


build:
	$(UV) build

format:
	$(UVX) ruff format

install:
	$(PIP) install -e .

publish:
	$(UV) publish
