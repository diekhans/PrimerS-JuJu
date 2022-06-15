# limit test, specify: tests="libtests/test_input.py::test_get_track_annot" 

ROOT = ..
include ${ROOT}/defs.mk


test: lib-tests prog-tests

fast-test: lib-tests-fast prog-tests-fast

lib-tests:
	${PYTEST} ${tests}

lib-tests-fast:
	${PYTEST} --ignore=libtests/test_design_primers.py ${tests}

prog-tests:
	cd progtests && ${MAKE} prog-tests

prog-tests-fast:
	cd progtests && ${MAKE} prog-tests-fast


lint:
	cd .. && ${MAKE} lint

clean:
	rm -rf output
	cd progtests && ${MAKE} clean
