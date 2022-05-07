PYTHON = python3
FLAKE8 = ${PYTHON} -m flake8
export PYTHONPATH = ${ROOT}/lib
TWINE = ${PYTHON} -m twine
VULTURE = vulture

PYTEST_FLAGS = -s -vv --tb=native
PYTEST = ${PYTHON} -m pytest ${PYTEST_FLAGS}
