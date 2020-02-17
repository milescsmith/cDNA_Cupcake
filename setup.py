from pathlib import Path
import numpy as np
from Cython.Build import cythonize
from setuptools import setup, Extension

__author__ = "etseng@pacb.com"
VERSION = "10.0.4"

EXT_MODULES = [
    Extension(
        "cupcake.tofu.branch.intersection_unique",
        ["cupcake/tofu/branch/intersection_unique.pyx"],
    ),
    Extension("cupcake.tofu.branch.c_branch", ["cupcake/tofu/branch/c_branch.pyx"]),
    Extension("cupcake.ice.find_ECE", ["cupcake/ice/find_ECE.pyx"]),
]

setup(
    name="cupcake",
    version=VERSION,
    author="Elizabeth Tseng",
    author_email="etseng@pacb.com",
    ext_modules=cythonize(EXT_MODULES),
    include_dirs=[np.get_include()],
    zip_safe=False,
    packages=[
        "annotation",
        "cupcake",
        "cupcake.io",
        "cupcake.ice",
        "cupcake.tofu",
        "cupcake.tofu.branch",
        "cupcake.tofu.counting",
        "cupcake2",
        "cupcake2.io",
        "cupcake2.ice2",
        "cupcake2.tofu2",
        "phasing",
        "phasing.io",
        "sequence",
    ],
    setup_requires=["numpy", "cython"],
    install_requires=[
        l.strip() for l in Path("requirements.txt").read_text("utf-8").splitlines()
    ],
    scripts=[
        "annotation/alignment_stats_from_sam.py",
        "annotation/make_file_for_subsampling_from_collapsed.py",
        "annotation/parse_matchAnnot.py",
        "annotation/subsample_with_category.py",
        "annotation/subsample.py",
        "cupcake/tofu/collapse_isoforms_by_sam.py",
        "cupcake/tofu/get_abundance_post_collapse.py",
        "cupcake/tofu/filter_by_count.py",
        "cupcake/tofu/filter_away_subset.py",
        "cupcake/tofu/fusion_finder.py",
        "cupcake/tofu/counting/chain_samples.py",
        "cupcake/tofu/counting/chain_fusion_samples.py",
        "cupcake/tofu/counting/summarize_sample_GFF_junctions.py",
        "cupcake/tofu/counting/scrub_sample_GFF_junctions.py",
        "cupcake2/tofu2/ice_pbdagcon2.py",
        "cupcake2/tofu2/run_preCluster.py",
        "cupcake2/tofu2/run_IceInit2.py",
        "cupcake2/tofu2/run_IceIterative2.py",
        "cupcake2/tofu2/run_IcePartial2.py",
        "cupcake2/tofu2/run_IceArrow2.py",
        "cupcake2/io/SeqSplitter.py",
        "cupcake2/tofu2/picking_up_ice2.py",
        "phasing/create_fake_genome.py",
        "phasing/run_phaser.py",
    ],
)
