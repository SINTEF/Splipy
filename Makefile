.PHONY: install mypy lint fmt fmtcheck doc

install:
	poetry install --with=dev

mypy:
	poetry run mypy splipy

lint:
	poetry run ruff splipy

bench:
	poetry run pytest --benchmark-only

fmt:
	poetry run black splipy
	poetry run isort splipy

fmtcheck:
	poetry run black splipy --check
	poetry run isort splipy --check


doc:
	$(MAKE) -C doc html


# Test targets

.PHONY: pytest
pytest:
	poetry run pytest --benchmark-skip

.PHONY: examples
examples:
	poetry run python examples/circle_animation.py --ci
	poetry run python examples/lissajous.py --ci
	poetry run python examples/loft.py
	poetry run python examples/read.py
	poetry run python examples/reuleaux.py --ci
	poetry run python examples/trefoil.py
	poetry run python examples/write.py

.PHONY: test  # most common test commands for everyday development
test: pytest

.PHONY: test-all  # run from CI: the whole kitchen sink
test-all: test examples


# Build targets (used from CI)

.PHONY: sdist
sdist:
	poetry build -f sdist

.PHONY: wheel
wheel:
	poetry build -f wheel

.PHONY: build
build: sdist wheel
