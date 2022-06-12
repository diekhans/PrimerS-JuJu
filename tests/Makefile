# limit test, specify: tests="libtests/test_input.py::test_get_track_annot" 

ROOT = ..
include ${ROOT}/defs.mk


test: lib-tests

fast-test: lib-tests-fast

lib-tests:
	${PYTEST} ${tests}

lib-tests-fast:
	${PYTEST} --ignore=libtests/test_design_primers.py ${tests}

lint:
	cd .. && ${MAKE} lint

clean:
	rm -rf output
