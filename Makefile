# Convenience targets

.PHONY: sync
sync:
	uv sync --all-packages --group dev

.PHONY: doc
doc:
	$(MAKE) -C doc html


# Other targets

.PHONY: bench
bench:
	uv run pytest --benchmark-only


# Linting targets

.PHONY: format
format:
	uv run ruff format

.PHONY: lint
lint:
	uv run ruff check --fix


# Test targets

.PHONY: pytest
pytest:
	uv run pytest --benchmark-skip

.PHONY: mypy
mypy:
	uv run mypy

.PHONY: ruff
ruff:
	uv run ruff check
	uv run ruff format --check

.PHONY: examples
examples:
	uv run python examples/circle_animation.py --ci
	uv run python examples/lissajous.py --ci
	uv run python examples/loft.py
	uv run python examples/read.py
	uv run python examples/reuleaux.py --ci
	uv run python examples/trefoil.py
	uv run python examples/write.py

.PHONY: test
test: pytest ruff examples


# Build targets (used from CI)

.PHONY: sdist
sdist:
	uv build --all-packages --sdist

.PHONY: wheel
wheel:
	uv build --all-packages --wheel

.PHONY: build
build: sdist wheel
