"""
Query for uniqueness in genome and transcriptome
"""
from dataclasses import dataclass
import pipettor
from pycbio.hgdata.bed import BedReader

@dataclass
class IsPcrServerSpec:
    """Specification of an isPCR server, static or dynamic;
    for genome or transcriptome"""
    host: str
    port: str
    target_2bit: str
    target_name: str = None  # target or transcriptome name for dynamic server
    target_data_dir: str = None   # dynamic blat data dir
    item_bigbed: str = None  # big bed file/URL for transcriptome

class UniquenessQuery:
    """Interface to UCSC isPCR server to query for uniqueness."""
    def __init__(self, genome_spec, transcriptome_spec):
        self.genome_spec = genome_spec
        self.transcriptome_spec = transcriptome_spec

    def _gfPcr(self, spec, name, left_primer, right_primer, max_size):
        cmd = ["gfPcr", f"-maxSize={max_size}", "-out=bed", f"-name={name}"]
        if spec.target_name is not None:
            cmd.append(f"-genome={spec.target_name}")
        if spec.target_data_dir is not None:
            cmd.extend(f"-genomeDataDir={spec.target_data_dir}")
        cmd.extend([spec.host, spec.port, spec.genome, left_primer, right_primer,
                    "/dev/stdout"])
        with pipettor.Popen(cmd) as fh:
            return [b for b in BedReader(fh)]

    def query_genome(self, name, left_primer, right_primer, max_size):
        return self._gfPcr(self.genome_spec, name, left_primer, right_primer, max_size)

    def query_transcriptome(self, name, left_primer, right_primer, max_size):
        beds = self._gfPcr(self.transcriptome_spec, name, left_primer, right_primer, max_size)
