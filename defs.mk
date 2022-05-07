MAKEFLAGS = --no-builtin-rules
SHELL = /bin/bash

devenv = ${ROOT}/dev-venv
testenv = ${ROOT}/test-venv

SYS_PYTHON = python3
PYTHON = source ${devenv}/bin/activate && ${devenv}/bin/python3
PIP = ${PYTHON} -m pip
FLAKE8 = ${PYTHON} -m flake8
export PYTHONPATH = ${ROOT}/lib
TWINE = ${PYTHON} -m twine
VULTURE = vulture

PYTEST_FLAGS = -s -vv --tb=native
PYTEST = ${PYTHON} -m pytest ${PYTEST_FLAGS}
