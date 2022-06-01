"""
Primer selection for a target.
"""
from typing import Sequence
from dataclasses import dataclass
from . import PrimersJuJuDataError
from .primer3_interface import Primer3Result, primer3_design
from primersjuju.target_transcripts import TargetTranscripts, TargetTranscript


@dataclass
class PrimerDesign:
    """information collected on one primer pair"""
    primer3_result: Primer3Result

@dataclass
class PrimerDesigns:
    """information collected on all primers for a target"""
    target_id: str
    target_transcripts: TargetTranscripts
    target_transcript: TargetTranscript
    designs: Sequence[PrimerDesign]

def _build_primer_designs(target_transcripts, target_transcript, primer3_results):
    return PrimerDesigns(target_transcripts.target_id,
                         target_transcripts, target_transcript,
                         [PrimerDesign(r) for r in primer3_results])

def design_primers(genome_data, target_transcripts, uniqueness_query=None):
    """design transcripts """
    target_transcript = target_transcripts.transcripts[0]
    primer3_results = primer3_design(target_transcript)

    return _build_primer_designs(target_transcripts, target_transcript,
                                 primer3_results)
