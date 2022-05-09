import re
import os.path as osp
from primersjuju.genome_data import GenomeData

mydir = osp.join(osp.dirname(__file__))

def _test_data_file(fname):
    return osp.join(mydir, "../data", fname)

_test_genome_data = None  # lazy load

def get_test_data():
    global _test_genome_data
    if _test_genome_data is None:
        _test_genome_data = GenomeData("hg38",
                                       _test_data_file("hg38.2bit"),
                                       "https://hgdownload.soe.ucsc.edu/gbdb/hg38/hg38.2bit")
        _test_genome_data.add_track("gencodeV39",
                                    _test_data_file("gencodeV39.bb"),
                                    "https://hgdownload.soe.ucsc.edu/gbdb/hg38/gencode/gencodeV39.bb")
        _test_genome_data.add_track("WTC11_consolidated",
                                    _test_data_file("WTC11_consolidated.bigBed"),
                                    "http://conesalab.org/LRGASP/LRGASP_hub/hg38/Human_samples/WTC11_consolidated.bigBed")

    return _test_genome_data


def assert_except_msg(exinfo, cls, msgre):
    "check causes"
    ex = exinfo.value
    while ex is not None:
        if isinstance(ex, cls) and re.match(msgre, str(ex)):
            return
        ex = ex.__cause__
    raise Exception(f"Caught exception can't match: {str(cls)} {msgre}") from exinfo.value
