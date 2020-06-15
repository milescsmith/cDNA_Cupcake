__author__ = "lachesis"

from .coordinate_mapper import (
    consistute_genome_seq_from_exons,
    get_base_to_base_mapping_from_sam,
    get_exon_coordinates,
    iter_cigar_string,
    make_exons_from_base_mapping,
)
from .MPileUpVariantCaller import MPileUPVariant
from .MummerSNPReader import SNPReader, SNPRecord, write_snp_to_vcf
from .SAMMPileUpReader import MPileUpReader, MPileUpRecord
from .VariantPhaseCleaner import (
    calc_hap_diff,
    error_correct_haplotypes,
    get_hap_model,
    infer_haplotypes_via_exhaustive_diploid_only,
    infer_haplotypes_via_min_diff,
    make_haplotype_counts,
)
from .VariantPhaser import Haplotypes, VariantPhaser, phase_isoforms, type_fa_or_fq
