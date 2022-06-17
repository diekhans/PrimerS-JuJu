
## Test data:

The directory data should contain:

* hg38.2bit https://hgdownload.soe.ucsc.edu/gbdb/hg38/hg38.2bit
* GCF_000001405.39_GRCh38.p13_assembly_report.txt https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt
* hg38 assembly report https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt
* gencodeV39.bb https://hgdownload.soe.ucsc.edu/gbdb/hg38/gencode/gencodeV39.bb
* hg38KgSeqV39.2bit https://hgdownload.soe.ucsc.edu/gbdb/hg38/targetDb/hg38KgSeqV39.2bit
* WTC11_consolidated.bigBed http://conesalab.org/LRGASP/LRGASP_hub/hg38/Human_samples/WTC11_consolidated.bigBed


## Notes

### JSON input for fast test.

Some tests have serialized the slow to produce isPcr results for later tests.  These
are in input/xxx.serial.json.  If these are removed, the will be regenerated when 
the tests run.  

See:
* libtests/test_design_primers_guts.py
