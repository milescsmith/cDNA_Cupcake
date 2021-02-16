from glob import glob
from pathlib import Path

import numpy as np
from Cython.Build import cythonize
from setuptools import Extension, find_packages, setup

try:
    exec(Path("src", "cupcake", "__about__.py").read_text())
except FileNotFoundError:
    __author__ = "Elizabeth Tseng"
    __email__ = "etseng@pacb.com"
    __version__ = "19.0.0"


EXT_MODULES = [
    Extension(
        "cupcake.tofu.branch.intersection_unique",
        ["src/cupcake/tofu/branch/intersection_unique.pyx"],
    ),
    Extension(
        "cupcake.tofu.branch.c_branch", ["src/cupcake/tofu/branch/c_branch.pyx"],
    ),
]

setup(
    name="cupcake",
    version=__version__,
    author=__author__,
    author_email=__email__,
    description="Miscellaneous collection of Python and R scripts for processing Iso-Seq data",
    long_description=Path("README.md").read_text("utf-8"),
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
    packages=find_packages("src"),
    package_data={"": ["tests/test_data/*.*"]},
    package_dir={"cupcake": "src/cupcake"},
    # py_modules=[Path(path).stem for path in glob("src/cupcake/*.py")],
    include_package_data=True,
    setup_requires=["numpy", "cython", "setuptools_scm"],
    # setup_requires=["cython", "setuptools_scm"],
    install_requires=[
        line.strip()
        for line in Path("requirements.txt").read_text("utf-8").splitlines()
    ],
    # this might be better to reformat into a series of subcommands
    entry_points={
        "console_scripts": [
            # annotation submodule
            "evaluate_alignment_sam      = cupcake.annotation.alignment_stats_from_sam:main",
            "make_file_for_subsample     = cupcake.annotation.make_file_for_subsampling_from_collapsed:main",
            "parse_matchAnnot            = cupcake.annotation.parse_matchAnnot:main",
            "subsample_with_category     = cupcake.annotation.subsample_with_category:main",
            "subsample                   = cupcake.annotation.subsample:main",
            # tofu submodule
            "collapse_isoforms_by_sam    = cupcake.tofu.collapse_isoforms_by_sam:main",
            "get_abundance_post_collapse = cupcake.tofu.get_abundance_post_collapse:main",
            "filter_by_count             = cupcake.tofu.filter_by_count:main",
            "filter_away_subset          = cupcake.tofu.filter_away_subset:main",
            "fusion_finder               = cupcake.tofu.fusion_finder:main",
            "chain_samples               = cupcake.tofu.counting.chain_samples:main",
            "chain_fusion_samples        = cupcake.tofu.counting.chain_fusion_samples:main",
            "summarize_junctions         = cupcake.tofu.counting.summarize_sample_GFF_junctions:main",
            "scrub_sample_GFFs           = cupcake.tofu.counting.scrub_sample_GFF_junctions:main",
            # uh, phasing! are we still doing phasing?
            "make_fake_genome            = cupcake.phasing.create_fake_genome:main",
            "simple_stats_post_collapse  = cupcake.tofu.simple_stats_post_collapse:main",
            "fusion_collate_info         = cupcake.tofu.fusion_collate_info:main",
            "color_bed12_post_sqanti     = cupcake.tofu.color_bed12_post_sqanti:main",
            "run_phaser                  = cupcake.phasing.run_phaser:main",
        ]
    },
)
