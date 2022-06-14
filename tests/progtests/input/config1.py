# configuration file with hg38

from os import path as osp
from primersjuju.config import PrimersJuJuConfig, GenomeConfig
from primersjuju.genome_data import GenomeData
from primersjuju.uniqueness_query import IsPcrServerSpec

data_dir = osp.join(osp.dirname(configPyFile), "../../data")

def data_path(fname):
    return osp.join(data_dir, fname)

hg38gd = GenomeData("hg38",
                    data_path("hg38.2bit"),
                    assembly_report=data_path("GCF_000001405.39_GRCh38.p13_assembly_report.txt"))
hg38gd.add_track("gencodeV39",
                 data_path("gencodeV39.bb"),
                 "https://hgdownload.soe.ucsc.edu/gbdb/hg38/gencode/gencodeV39.bb")
hg38gd.add_track("WTC11_consolidated",
                 data_path("WTC11_consolidated.bigBed"),
                 "http://conesalab.org/LRGASP/LRGASP_hub/hg38/Human_samples/WTC11_consolidated.bigBed")

hg38config = GenomeConfig(hg38gd,
                          genome_ispcr_spec=IsPcrServerSpec("blat1d.soe.ucsc.edu", 17903, data_dir),
                          transcriptome_ispcr_spec=IsPcrServerSpec("blat1a.soe.ucsc.edu", 17907, data_dir,  # hg38KgSeqV39
                                                                   trans_bigbed=data_path("gencodeV39.bb")))
config = PrimersJuJuConfig()
config.add_genome(hg38config)
