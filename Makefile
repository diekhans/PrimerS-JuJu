ROOT = .
include ${ROOT}/defs.mk

pyprogs = $(shell file -F $$'\t' bin/* | awk '/Python script/{print $$1}')
pytests =  tests/libtests/*.py
pypi_url = https://upload.pypi.org/simple/
testpypi_url = https://test.pypi.org/simple/

version = $(shell PYTHONPATH=lib ${PYTHON} -c "import primersjuju; print(primersjuju.__version__)")

# mdl is an uncommon program to verify markdown
have_mdl = $(shell which mdl >&/dev/null && echo yes || echo no)
have_mdlinkcheck = $(shell which markdown-link-check >&/dev/null && echo yes || echo no)


help:
	@echo "clean - remove all build, test, coverage and Python artifacts"
	@echo "devenv - setup local venv"
	@echo "doc - build various pieces of the doc"
	@echo "lint - check style with flake8"
	@echo "lint-doc - check documentation"
	@echo "lint-all - lint plus lint-doc"
	@echo "lint-pages - lint github pages, must have been pushed first"
	@echo "vulture - find unused code"
	@echo "test - run tests quickly with the default Python"
	@echo "install - install the package to the active Python's site-packages"
	@echo "dist - package"
	@echo "test-${PIP} - test install the package using pip"
	@echo "release-testpypi - test upload to testpypi"
	@echo "test-release-testpypi - install from testpypi"
	@echo "release - package and upload a release"
	@echo "test-release - test final release"


lint:
	${FLAKE8} ${pyprogs} ${pytests} lib

.PHONY: devenv
devenv:
	test -x ${devenv} || ${SYS_PYTHON} -m venv ${devenv}
	${PIP} install --upgrade pip
	${PIP} install --upgrade wheel
	${PIP} install -r requirements-dev.txt
	${PIP} install -e .

# this gets a lot of false-positive, just use for code cleanup rather than making it
# standard
vulture:
	${VULTURE} ${pyprogs} ${pytests} lib

test:
	cd tests && ${MAKE} test

clean: test_clean
	rm -rf build/ dist/ ${testenv}/ lib/*.egg-info/ lib/primersjuju/__pycache__/ dev-venv/

test_clean:
	cd tests && ${MAKE} clean

define envsetup
	@rm -rf ${testenv}/
	mkdir -p ${testenv}
	${PYTHON} -m virtualenv ${testenv}
endef
envact = source ${testenv}/bin/activate

dist_tar = dist/primersjuju-${version}.tar.gz
dist_whl = dist/primersjuju-${version}-py3-none-any.whl
pkgver_spec = primersjuju-tools==${version}

dist: clean
	${PYTHON} setup.py sdist
	${PYTHON} setup.py bdist_wheel
	@ls -l ${dist_tar}
	@ls -l ${dist_whl}

# test install locally
test-pip: dist
	${envsetup}
	${envact} && cd ${testenv} && ${PIP} install --no-cache-dir $(realpath ${dist_tar})
	${envact} && cd tests ${MAKE} test

# test release to testpypi
release-testpypi: dist
	${TWINE} upload --repository=testpypi ${dist_whl} ${dist_tar}

# test release install from testpypi
test-release-testpypi:
	${envsetup}
	${envact} && cd ${testenv} && ${PIP} install --no-cache-dir --index-url=${testpypi_url} --extra-index-url=https://pypi.org/simple ${pkgver_spec}
	${envact} && cd tests && ${MAKE} test

release: dist
	${TWINE} upload --repository=pypi ${dist_whl} ${dist_tar}

test-release:
	${envsetup}
	${envact} && cd ${testenv} && ${PIP} install --no-cache-dir ${pkgver_spec}
	${envact} && cd tests && ${MAKE} test

