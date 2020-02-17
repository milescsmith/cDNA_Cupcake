from .evaluate_phase_switch import (
    read_config,
    eval_isophase_phaseswitch,
    main_eval,
)

from .evaluate_snp_with_genome import (
    read_fake_mapping,
    get_positions_to_recover,
    eval_isophase,
    main_brangus,
    main_maize,
)

from .select_loci_to_phase import (
    read_flnc_fastq,
    read_read_stat,
    read_GFF,
    make_fake_genome,
    select_loci_to_phase,
    getargs,
)