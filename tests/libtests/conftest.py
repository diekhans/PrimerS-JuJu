import pytest
import os.path as osp
from primersjuju.genome_data import GenomeData
from . import mydir
from primersjuju.primer_target_spec import primer_targets_specs_read

def _test_data_file(fname):
    return osp.join(mydir, "../data", fname)

@pytest.fixture(scope="session")
def gdata():
    test_gdata = GenomeData("hg38",
                            _test_data_file("hg38.2bit"),
                            "https://hgdownload.soe.ucsc.edu/gbdb/hg38/hg38.2bit")
    test_gdata.add_track("gencodeV39",
                         _test_data_file("gencodeV39.bb"),
                         "https://hgdownload.soe.ucsc.edu/gbdb/hg38/gencode/gencodeV39.bb")
    test_gdata.add_track("WTC11_consolidated",
                         _test_data_file("WTC11_consolidated.bigBed"),
                         "http://conesalab.org/LRGASP/LRGASP_hub/hg38/Human_samples/WTC11_consolidated.bigBed")

    return test_gdata

@pytest.fixture(scope="session")
def example_targets_specs():
    return primer_targets_specs_read(osp.join(mydir, "../../docs/primer-targets-example.tsv"))


@pytest.fixture(scope="session")
def example_wtc11_targets_specs_set1():
    return primer_targets_specs_read(osp.join(mydir, "../input/LRGASP_WTC11_Target_set1.tsv"))
