# limit test, specify:
# make lib-tests tests="libtests/test_input.py::test_get_track_annot" 

ROOT = ..
include ${ROOT}/defs.mk


test: lib-tests prog-tests

test-fast: lib-tests-fast prog-tests-fast

lib-tests:
	cd libtests && ${MAKE} $@

lib-tests-fast:
	cd libtests && ${MAKE} $@

prog-tests:
	cd progtests && ${MAKE} $@

prog-tests-fast:
	cd progtests && ${MAKE} $@


lint:
	cd ${ROOT} && ${MAKE} lint

clean:
	cd libtests && ${MAKE} $@
	cd progtests && ${MAKE} $@

update:
	cd libtests/output && (find . -type f | cpio -pdum ../expected)
	cd progtests/output && (find . -type f | cpio -pdum ../expected)
