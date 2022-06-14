ROOT = ../..
include ${ROOT}/defs.mk

DEBUG = --debug

prog-tests: prog-tests-fast prog-tests-slow

prog-tests-fast: test-check-good

prog-tests-slow: test-design-good3 test-design-subset

test-check-good:
	${primers_juju} ${DEBUG} --conf=input/config1.py check ../input/LRGASP_WTC11_Target_set1.tsv

test-design-good3: mkdirs
	rm -rf output/$@
	${primers_juju} ${DEBUG} --conf=input/config1.py design input/LRGASP_WTC11_Target_subset3.tsv output/$@

test-design-subset: mkdirs
	rm -rf output/$@
	${primers_juju} ${DEBUG} --target=BBC3+1 --conf=input/config1.py design input/LRGASP_WTC11_Target_subset3.tsv output/$@


mkdirs:
	@mkdir -p output
