
## Test data:

The directory data/ should contain:

* hg38.2bit https://hgdownload.soe.ucsc.edu/gbdb/hg38/hg38.2bit
* GCF_000001405.40_GRCh38.p14_assembly_report.txt https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt
* gencodeV41.bb https://hgdownload.soe.ucsc.edu/gbdb/hg38/gencode/gencodeV41.bb
* hg38KgSeqV41.2bit https://hgdownload.soe.ucsc.edu/gbdb/hg38/targetDb/hg38KgSeqV41.2bit
* WTC11_consolidated.bigBed http://conesalab.org/LRGASP/LRGASP_hub/hg38/Human_samples/WTC11_consolidated.bigBed


## Notes

    gfServer start hgwdev.gi.ucsc.edu 12201 -stepSize=5 -log=lrgasp-hg38-pcr.log /hive/users/markd/gencode/projs/lrgasp/primers/primer-design/data/hg38/hg38_transcriptome.2bit &
