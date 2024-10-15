ROOT = ../..
include ${ROOT}/defs.mk

DEBUG = --debug

hg38_hub_url_args = --huburl=http://conesalab.org/LRGASP/LRGASP_hub/hub.txt --huburl=https://hgwdev.gi.ucsc.edu/~markd/lrgasp/juju-hub/hub.txt
hg38_twobit = https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit

suffixes = .debug.txt designs.tsv primers.bed target.bed uniqueness.bed

define diff_results =
@for f in $$(cd output && find $@ -type f) ; do diff -u expected/$$f output/$$f ; done
@tawk 'NF != 18{print "Expect 18 columns, got " NF ":" NR " in "FILENAME; exit 1}' output/$@/*/*.designs.tsv
endef

test: prog-tests
test-fast: prog-tests-fast

prog-tests: prog-tests-fast prog-tests-slow

prog-tests-fast: test-check-good test-design-subset test-design-skip-error test-disabled

prog-tests-slow: test-design-good3

test-check-good:
	${primers_juju} ${DEBUG} --conf=input/config1.py check ../input/LRGASP_WTC11_Target_set1.tsv

# also check if bigBed can be built
test-design-good3: mkdirs
	rm -rf output/$@
	${primers_juju} ${DEBUG} --conf=input/config1.py design ${hg38_hub_url_args} input/LRGASP_WTC11_Target_subset3.tsv output/$@
	${diff_results}
	bedToBigBed -tab -type=bed12+ -as=${ROOT}/etc/bigPrimerBed.as -sizesIs2Bit -extraIndex=name \
	   output/test-design-good3/ZBTB45+1/ZBTB45+1.primers.bed ${hg38_twobit} output/$@/ZBTB45+1.bb
	rm output/$@/ZBTB45+1.bb

test-design-subset: mkdirs
	rm -rf output/$@
	${primers_juju} ${DEBUG} --conf=input/config1.py design --target=BBC3+1 ${hg38_hub_url_args} input/LRGASP_WTC11_Target_subset3.tsv output/$@
	${diff_results}

test-design-skip-existing: mkdirs
	rm -rf output/$@
	${primers_juju} ${DEBUG} --conf=input/config1.py design --target=BBC3+1 ${hg38_hub_url_args} input/LRGASP_WTC11_Target_subset3.tsv output/$@ >& output/$@.design.out
	${primers_juju} ${DEBUG} --conf=input/config1.py design --target=BBC3+1 ${hg38_hub_url_args} input/LRGASP_WTC11_Target_subset3.tsv output/$@ >& output/$@.skip.out
	diff expected/$@.design.out output/$@.design.out
	diff expected/$@.skip.out output/$@.skip.out

test-design-skip-error: mkdirs
	rm -rf output/$@
	${primers_juju}  --conf=input/config1.py design ${hg38_hub_url_args} input/skip_error.tsv output/$@ >& output/$@.design.out  && (echo "Error: BBC3+2 should have failed" 1>&2; exit 1) || true
	@test -e output/$@/BBC3+2/BBC3+2.designs.tsv && (echo "Error: BBC3+2 show not have produced output/$@/BBC3+2/BBC3+2.designs.tsv file" 1>&2; exit 1) || true
	diff expected/$@.design.out output/$@.design.out
	${diff_results}

test-disabled: mkdirs
	rm -rf output/$@
	${primers_juju} ${DEBUG} --conf=input/config1.py design ${hg38_hub_url_args} input/disabled.tsv output/$@
	${diff_results}


mkdirs:
	@mkdir -p output

clean:
	rm -rf output
