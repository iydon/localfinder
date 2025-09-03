UV     := uv
PIP    := pip
PYTHON := $(UVX) python


.PHONY: \
	build install publish  \
	mypy ruff


build:
	$(UV) build

install:
	# conda create --name lf python=3.9
	# conda install -c bioconda -c conda-forge samtools bedtools ucsc-bigwigtobedgraph
	# conda install -c bioconda -c conda-forge igvtools=2.16.2
	$(PIP) install -e .

publish:
	$(UV) publish


mypy:
	# uv tool install mypy --with "pip,matplotlib"
	# uv tool run mypy localfinder --install-types
	$(UV) tool run mypy localfinder

ruff:
	# uv tool install ruff
	$(UV) tool run ruff format
