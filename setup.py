from pathlib import Path

import numpy as np
from setuptools import Extension, setup, find_packages

from Cython.Build import cythonize

__author__ = "etseng@pacb.com"
version = "12.1.2"

ext_modules = [
    Extension(
        " = cdna_cupcake.cupcake.tofu.branch.intersection_unique",
        ["cdna_cupcake/cupcake/tofu/branch/intersection_unique.pyx"],
    ),
    Extension(
        "cdna_cupcake.cupcake.tofu.branch.c_branch",
        ["cdna_cupcake/cupcake/tofu/branch/c_branch.pyx"],
    ),
    Extension(
        "cdna_cupcake.cupcake.ice.find_ECE", ["cdna_cupcake/cupcake/ice/find_ECE.pyx"]
    ),
]

EXT_MODULES = [
    Extension(
        " = cdna_cupcake.cupcake.tofu.branch.intersection_unique",
        ["cdna_cupcake/cupcake/tofu/branch/intersection_unique.pyx"],
    ),
    Extension(
        "cdna_cupcake.cupcake.tofu.branch.c_branch",
        ["cdna_cupcake/cupcake/tofu/branch/c_branch.pyx"],
    ),
    Extension(
        "cdna_cupcake.cupcake.ice.find_ECE", ["cdna_cupcake/cupcake/ice/find_ECE.pyx"]
    ),
]

setup(
    name="cupcake",
    version=version,
    author="Elizabeth Tseng",
    author_email="etseng@pacb.com",
    description="Miscellaneous collection of Python and R scripts for processing Iso-Seq data",
    long_description=Path("README.rst").read_text("utf-8"),
    url="https://github.com/Magdoll/cDNA_Cupcake",
    ext_modules=cythonize(EXT_MODULES),
    include_dirs=[np.get_include()],
    zip_safe=False,
    python_requires="~=3.7",
    license="BSD3",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords=["isoseq", "rnaseq", "pacbio", "long reads"],
    packages=find_packages(include=["cdna_cupcake", "cdna_cupcake.*"]),
    package_dir={"cdna_cupcake": "cdna_cupcake"},
    package_data={"": ["cupcake/test_data/*.*",]},
    include_package_data=True,
    setup_requires=["numpy", "cython"],
    install_requires=[
        l.strip() for l in Path("requirements.txt").read_text("utf-8").splitlines()
    ],
    # this might be better to reformat into a series of subcommands
    entry_points={
        "console_scripts": [
            "evaluate_alignment_sam = annotation.alignment_stats_from_sam.evaluate_alignment_sam:evaluate_alignment_sam",
            "make_file_for_subsample = cdna_cupcake.annotation.make_file_for_subsampling_from_collapsed:make_file_for_subsample",
            "parse_matchAnnot = cdna_cupcake.annotation.parse_matchAnnot:parse_matchAnnot",
            "subsample_with_category = cdna_cupcake.annotation.subsample_with_category:subsample",
            "subsample = cdna_cupcake.annotation.subsample:subsample",
            "collapse_isoforms_by_sam = cdna_cupcake.cupcake.tofu.collapse_isoforms_by_sam:main",
            "get_abundance_post_collapse = cdna_cupcake.cupcake.tofu.get_abundance_post_collapse:get_abundance_post_collapse",
            "filter_by_count = cdna_cupcake.cupcake.tofu.filter_by_count:filter_by_count",
            "filter_away_subset = cdna_cupcake.cupcake.tofu.filter_away_subset:main",
            "fusion_finder = cdna_cupcake.cupcake.tofu.fusion_finder:fusion_main",
            "chain_samples = cdna_cupcake.cupcake.tofu.counting.chain_samples:chain_samples_multithread",
            "chain_fusion_samples = cdna_cupcake.cupcake.tofu.counting.chain_fusion_samples:chain_fusion_samples",
            "summarize_junctions = cdna_cupcake.cupcake.tofu.counting.summarize_sample_GFF_junctions:summarize_junctions",
            "scrub_sample_GFFs = cdna_cupcake.cupcake.tofu.counting.scrub_sample_GFF_junctions:scrub_sample_GFFs",
            "run_Consensus = cdna_cupcake.cupcake2.tofu2.ice_pbdagcon2:runConsensus",
            "run_preCluster = cdna_cupcake.cupcake2.tofu2.run_preCluster:main",
            "run_IceInit2 = cdna_cupcake.cupcake2.tofu2.run_IceInit2:run_IceInit2",
            "run_IceIterative2 = cdna_cupcake.cupcake2.tofu2.run_IceIterative2:run_IceIterative2",
            "run_IcePartial2 = cdna_cupcake.cupcake2.tofu2.run_IcePartial2:main",
            "run_IceArrow2 = cdna_cupcake.cupcake2.tofu2.run_IceArrow2:main",
            "SeqSplitter = cdna_cupcake.cupcake2.io.SeqSplitter:main",
            "picking_up_ice2 = cdna_cupcake.cupcake2.tofu2.picking_up_ice2:main",
            "make_fake_genome = cdna_cupcake.phasing.create_fake_genome:main",
            # " = cdna_cupcake.phasing.run_phaser",
            # " = cdna_cupcake.sequence.sam_to_bam",
            # " = cdna_cupcake.sequence.err_correct_w_genome",
            # " = cdna_cupcake.sequence.sam_to_gff3",
            # " = cdna_cupcake.sequence.STAR",
            # " = cdna_cupcake.sequence.BED",
            # " = cdna_cupcake.sequence.coordinate_mapper",
        ],
    },
)
