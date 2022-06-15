ROOT = ../..
include ${ROOT}/defs.mk

DEBUG = --debug

hg38_hub_url_args = --huburl=http://conesalab.org/LRGASP/LRGASP_hub/hub.txt --huburl=https://hgwdev.gi.ucsc.edu/~markd/lrgasp/juju-hub/hub.txt

prog-tests: prog-tests-fast prog-tests-slow

prog-tests-fast: test-check-good

prog-tests-slow: test-design-good3 test-design-subset

test-check-good:
	${primers_juju} ${DEBUG} --conf=input/config1.py check ../input/LRGASP_WTC11_Target_set1.tsv

test-design-good3: mkdirs
	rm -rf output/$@
	${primers_juju} ${DEBUG} ${hg38_hub_url_args} --conf=input/config1.py design input/LRGASP_WTC11_Target_subset3.tsv output/$@

test-design-subset: mkdirs
	rm -rf output/$@
	${primers_juju} ${DEBUG} ${hg38_hub_url_args} --target=BBC3+1 --conf=input/config1.py design input/LRGASP_WTC11_Target_subset3.tsv output/$@


mkdirs:
	@mkdir -p output

clean:
	rm -rf output
