# -*- coding: utf-8 -*-
from setuptools import setup

from build import *

package_dir = {"": "src"}

packages = [
    "cupcake",
    "cupcake.SequelQC",
    "cupcake.annotation",
    "cupcake.bacteria",
    "cupcake.beta",
    "cupcake.phasing",
    "cupcake.phasing.io",
    "cupcake.phasing.utils",
    "cupcake.post_isoseq_cluster",
    "cupcake.sequence",
    "cupcake.simulate",
    "cupcake.singlecell",
    "cupcake.targeted",
    "cupcake.tofu",
    "cupcake.tofu.branch",
    "cupcake.tofu.counting",
]

package_data = {"": ["*"], "cupcake.annotation": ["test_data/*"]}

install_requires = [
    "Cython>=0.29.21,<0.30.0",
    "bcbio-gff>=0.6.6,<0.7.0",
    "biopython>=1.78,<2.0",
    "bx-python>=0.8.9,<0.9.0",
    "numpy>=1.20.1,<2.0.0",
    "pandas>=1.2.2,<2.0.0",
    "pbcommand @ git+https://github.com/PacificBiosciences/pbcommand@master",
    "pbcore @ git+https://github.com/PacificBiosciences/pbcore@master",
    "pbcoretools @ git+https://github.com/PacificBiosciences/pbcoretools@master",
    "pysam>=0.16.0,<0.17.0",
    "pyvcf @ git+https://github.com/milescsmith/PyVCF@master",
    "scikit-learn>=0.24.1,<0.25.0",
]

entry_points = {
    "console_scripts": [
        "chain_fusion_samples = " "cupcake.tofu.counting.chain_fusion_samples:main",
        "chain_samples = cupcake.tofu.counting.chain_samples:main",
        "collapse_isoforms_by_sam = " "cupcake.tofu.collapse_isoforms_by_sam:main",
        "color_bed12_post_sqanti = " "cupcake.tofu.color_bed12_post_sqanti:main",
        "evaluate_alignment_sam = " "cupcake.annotation.alignment_stats_from_sam:main",
        "filter_away_subset = " "cupcake.tofu.filter_away_subset:main",
        "filter_by_count = cupcake.tofu.filter_by_count:main",
        "fusion_collate_info = " "cupcake.tofu.fusion_collate_info:main",
        "fusion_finder = cupcake.tofu.fusion_finder:main",
        "get_abundance_post_collapse = "
        "cupcake.tofu.get_abundance_post_collapse:main",
        "make_fake_genome = " "cupcake.phasing.create_fake_genome:main",
        "make_file_for_subsample = "
        "cupcake.annotation.make_file_for_subsampling_from_collapsed:main",
        "parse_matchAnnot = " "cupcake.annotation.parse_matchAnnot:main",
        "run_phaser = cupcake.phasing.run_phaser:main",
        "scrub_sample_GFFs = " "cupcake.tofu.counting.scrub_sample_GFF_junctions:main",
        "simple_stats_post_collapse = " "cupcake.tofu.simple_stats_post_collapse:main",
        "subsample = cupcake.annotation.subsample:main",
        "subsample_with_category = " "cupcake.annotation.subsample_with_category:main",
        "summarize_junctions = "
        "cupcake.tofu.counting.summarize_sample_GFF_junctions:main",
    ]
}

setup_kwargs = {
    "name": "cupcake",
    "version": "19.0.2",
    "description": "Miscellaneous collection of Python and R scripts for processing Iso-Seq data",
    "long_description": None,
    "author": "Elizabeth Tseng",
    "author_email": "etseng@pacb.com",
    "maintainer": None,
    "maintainer_email": None,
    "url": None,
    "package_dir": package_dir,
    "packages": packages,
    "package_data": package_data,
    "install_requires": install_requires,
    "entry_points": entry_points,
    "python_requires": ">=3.7.9,<4.0.0",
}

build(setup_kwargs)

setup(**setup_kwargs)
