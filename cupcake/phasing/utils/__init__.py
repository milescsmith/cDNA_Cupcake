from .evaluate_phase_switch import eval_isophase_phaseswitch, main_eval, read_config
from .evaluate_snp_with_genome import (
    eval_isophase,
    get_positions_to_recover,
    main_brangus,
    main_maize,
    read_fake_mapping,
)
from .select_loci_to_phase import (
    getargs,
    make_fake_genome,
    read_flnc_fastq,
    read_GFF,
    read_read_stat,
    select_loci_to_phase,
)
