"""
Query for uniqueness in genome and transcriptome
"""
from dataclasses import dataclass
import pipettor
from pycbio.hgdata.psl import PslReader
from .genome_data import big_bed_read_by_names

@dataclass
class IsPcrServerSpec:
    """Specification of an isPCR server, static or dynamic;
    for genome or transcriptome"""
    host: str
    port: str
    target_seq_dir: str  # contains 2bit matching one returned by server
    dyn_name: str = None  # target or transcriptome name for dynamic server
    dyn_data_dir: str = None   # dynamic blat data dir
    item_bigbed: str = None  # big bed file/URL for transcriptome

class UniquenessQuery:
    """Interface to UCSC isPCR server to query for uniqueness."""
    def __init__(self, genome_spec, transcriptome_spec):
        self.genome_spec = genome_spec
        self.transcriptome_spec = transcriptome_spec

    def _gfPcr(self, spec, name, left_primer, right_primer, max_size):
        "returns PSL records"
        cmd = ["gfPcr", f"-maxSize={max_size}", "-out=psl", f"-name={name}"]
        if spec.dyn_name is not None:
            cmd.append(f"-genome={spec.dyn_name}")
        if spec.dyn_data_dir is not None:
            cmd.extend(f"-genomeDataDir={spec.dyn_data_dir}")
        cmd.extend([spec.host, spec.port, spec.target_seq_dir, left_primer, right_primer, "/dev/stdout"])
        with pipettor.Popen(cmd) as fh:
            return [p for p in PslReader(fh)]

    def query_genome(self, name, left_primer, right_primer, max_size):
        genome_psls = self._gfPcr(self.genome_spec, name, left_primer, right_primer, max_size)
        return genome_psls

    def query_transcriptome(self, name, left_primer, right_primer, max_size):
        trans_psls = self._gfPcr(self.transcriptome_spec, name, left_primer, right_primer, max_size)
        return trans_psls
