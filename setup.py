from pathlib import Path

import numpy as np
from setuptools import Extension, setup, find_packages

from Cython.Build import cythonize

__author__ = "etseng@pacb.com"

# why is this defined twice?  commenting out until I figure out why
# ext_modules = [
#     Extension(
#         " = cupcake.cupcake.tofu.branch.intersection_unique",
#         ["cdna_cupcake/cupcake/tofu/branch/intersection_unique.pyx"],
#     ),
#     Extension(
#         "cupcake.cupcake.tofu.branch.c_branch",
#         ["cdna_cupcake/cupcake/tofu/branch/c_branch.pyx"],
#     ),
#     Extension(
#         "cupcake.cupcake.ice.find_ECE", ["cdna_cupcake/cupcake/ice/find_ECE.pyx"]
#     ),
# ]

EXT_MODULES = [
    Extension(
        "cupcake.cupcake.tofu.branch.intersection_unique",
        ["cupcake/cupcake/tofu/branch/intersection_unique.pyx"],
    ),
    Extension(
        "cupcake.cupcake.tofu.branch.c_branch",
        ["cupcake/cupcake/tofu/branch/c_branch.pyx"],
    ),
    Extension("cupcake.cupcake.ice.find_ECE", ["cupcake/cupcake/ice/find_ECE.pyx"]),
]

setup(
    name="cupcake",
    version="12.2.1",
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
    packages=find_packages(),
    package_dir={"cupcake": "cupcake"},
    package_data={"": ["cupcake/cupcake/test_data/*.*"]},
    include_package_data=True,
    setup_requires=["numpy", "cython"],
    install_requires=[
        l.strip() for l in Path("requirements.txt").read_text("utf-8").splitlines()
    ],
    # this might be better to reformat into a series of subcommands
    entry_points={
        "console_scripts": [
            # annotation submodule
            "evaluate_alignment_sam = cupcake.annotation.alignment_stats_from_sam.evaluate_alignment_sam:main",
            "make_file_for_subsample = cupcake.annotation.make_file_for_subsampling_from_collapsed:main",
            "parse_matchAnnot = cupcake.annotation.parse_matchAnnot:main",
            "subsample_with_category = cupcake.annotation.subsample_with_category:main",
            "subsample = cupcake.annotation.subsample:main",
            # cupcake submodule
            "collapse_isoforms_by_sam = cupcake.cupcake.tofu.collapse_isoforms_by_sam:main",
            "get_abundance_post_collapse = cupcake.cupcake.tofu.get_abundance_post_collapse:main",
            "filter_by_count = cupcake.cupcake.tofu.filter_by_count:main",
            "filter_away_subset = cupcake.cupcake.tofu.filter_away_subset:main",
            "fusion_finder = cupcake.cupcake.tofu.fusion_finder:main",
            "chain_samples = cupcake.cupcake.tofu.counting.chain_samples:main",
            "chain_fusion_samples = cupcake.cupcake.tofu.counting.chain_fusion_samples:main",
            "summarize_junctions = cupcake.cupcake.tofu.counting.summarize_sample_GFF_junctions:main",
            "scrub_sample_GFFs = cupcake.cupcake.tofu.counting.scrub_sample_GFF_junctions:main",
            # cupcake2 submodule
            "run_Consensus = cupcake.cupcake2.tofu2.ice_pbdagcon2:runConsensus",
            "run_preCluster = cupcake.cupcake2.tofu2.run_preCluster:main",
            "run_IceInit2 = cupcake.cupcake2.tofu2.run_IceInit2:main",
            "run_IceIterative2 = cupcake.cupcake2.tofu2.run_IceIterative2:main",
            "run_IcePartial2 = cupcake.cupcake2.tofu2.run_IcePartial2:main",
            "run_IceArrow2 = cupcake.cupcake2.tofu2.run_IceArrow2:main",
            "SeqSplitter = cupcake.cupcake2.io.SeqSplitter:main",
            "picking_up_ice2 = cupcake.cupcake2.tofu2.picking_up_ice2:main",
            # uh, phasing! are we still doing phasing?
            "make_fake_genome = cupcake.phasing.create_fake_genome:main",
            # " = cupcake.phasing.run_phaser",
            # " = cupcake.sequence.sam_to_bam",
            # " = cupcake.sequence.err_correct_w_genome",
            # " = cupcake.sequence.sam_to_gff3",
            # " = cupcake.sequence.STAR",
            # " = cupcake.sequence.BED",
            # " = cupcake.sequence.coordinate_mapper",
        ]
    },
)
