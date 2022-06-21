MAKEFLAGS = --no-builtin-rules
SHELL = /bin/bash -e

pycbio = ${HOME}/compbio/code/pycbio

dev_venv = ${ROOT}/dev-venv
dev_venv_acct = source ${dev_venv}/bin/activate

SYS_PYTHON = python3
PYTHON = ${dev_venv_acct} && python3
PIP = ${PYTHON} -m pip
FLAKE8 = ${PYTHON} -m flake8
export PYTHONPATH = ${ROOT}/lib
TWINE = ${PYTHON} -m twine
VULTURE = ${PYTHON} -m vulture

PYTEST_FLAGS = -s -vv --tb=native
PYTEST = ${PYTHON} -m pytest ${PYTEST_FLAGS}

primers_juju = ${dev_venv_acct} && ${ROOT}/bin/primers-juju
