# limit test, specify:
# make lib-tests tests="test_input.py::test_get_track_annot" 

ROOT = ../..
include ${ROOT}/defs.mk


test: lib-tests

test-fast: lib-tests-fast

lib-tests:
	${PYTEST} ${tests}

lib-tests-fast:
	${PYTEST} --ignore=test_design_primers.py ${tests}

lint:
	cd ${ROOT} && ${MAKE} lint

clean:
	rm -rf output
