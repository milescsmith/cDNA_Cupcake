__author__ = "lachesis"

from .coordinate_mapper import (
    iter_cigar_string,
    make_exons_from_base_mapping,
    get_base_to_base_mapping_from_sam,
    get_exon_coordinates,
    consistute_genome_seq_from_exons,
)

from .MPileUpVariantCaller import MPileUPVariant

from .MummerSNPReader import (
    SNPRecord,
    SNPReader,
    write_snp_to_vcf
)

from .SAMMPileUpReader import (
    MPileUpRecord,
    MPileUpReader
)

from .VariantPhaser import (
    type_fa_or_fq,
    VariantPhaser,
    phase_isoforms,
    Haplotypes,
)

from .VariantPhaseCleaner import (
    make_haplotype_counts,
    calc_hap_diff,
    get_hap_model,
    infer_haplotypes_via_exhaustive_diploid_only,
    infer_haplotypes_via_min_diff,
    error_correct_haplotypes,
)