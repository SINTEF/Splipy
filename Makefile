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


# Test targets

.PHONY: pytest
pytest:
	uv run pytest --benchmark-skip

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
test: pytest examples


# Build targets (used from CI)

.PHONY: sdist
sdist:
	uv build --all-packages --sdist

.PHONY: wheel
wheel:
	uv build --all-packages --wheel

.PHONY: build
build: sdist wheel
