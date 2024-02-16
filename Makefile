package := splipy


# Convenience targets

.PHONY: install
install:
	poetry install --with=dev

.PHONY: doc
doc:
	$(MAKE) -C doc html


# Linting targets

.PHONY: format
format:
	poetry run ruff format $(package)

.PHONY: lint
lint:
	poetry run ruff check --fix $(package)


# Test targets

.PHONY:
benchmark:
	poetry run pytest --benchmark-only

.PHONY: pytest
pytest:
	poetry run pytest --benchmark-skip

.PHONY: mypy
mypy:
	poetry run mypy $(package)

.PHONY: lint-check
lint-check:
	poetry run ruff check $(package)
	poetry run ruff format --check $(package)

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
test: pytest mypy lint-check

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
