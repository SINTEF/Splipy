.PHONY: install mypy lint pytest test fmt fmtcheck doc

install:
	poetry install --with=dev

mypy:
	poetry run mypy splipy

lint:
	poetry run ruff splipy

pytest:
	poetry run pytest --benchmark-skip

bench:
	poetry run pytest --benchmark-only

fmt:
	poetry run black splipy
	poetry run isort splipy

fmtcheck:
	poetry run black splipy --check
	poetry run isort splipy --check

test: pytest

doc:
	$(MAKE) -C doc html


# Build targets (used from CI)

sdist:
	poetry build -f sdist

wheel:
	poetry build -f wheel

build: sdist wheel
