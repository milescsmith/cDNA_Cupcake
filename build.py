# lifted from https://github.com/sdispater/pendulum
import shutil
from distutils.command.build_ext import build_ext
from distutils.core import Distribution, Extension
from distutils.errors import CCompilerError, DistutilsExecError, DistutilsPlatformError
from pathlib import Path

import numpy as np
from Cython.Build import cythonize

# C Extensions

extensions = [
    Extension(
        "cupcake.tofu.branch.intersection_unique",
        include_dirs=[np.get_include()],
        sources=["src/cupcake/tofu/branch/intersection_unique.pyx"],
    ),
    Extension(
        "cupcake.tofu.branch.c_branch",
        include_dirs=[np.get_include()],
        sources=["src/cupcake/tofu/branch/c_branch.pyx"],
    ),
]


class BuildFailed(Exception):
    pass


class ExtBuilder(build_ext):
    # This class allows C extension building to fail.

    built_extensions = []

    def run(self):
        try:
            build_ext.run(self)
        except (DistutilsPlatformError, FileNotFoundError):
            print(
                "  Unable to build the C extensions, "
                "cDNA_Cupcake will use the pure python code instead."
            )

    def build_extension(self, ext):
        try:
            build_ext.build_extension(self, ext)
        except (CCompilerError, DistutilsExecError, DistutilsPlatformError, ValueError):
            print(
                f"  Unable to build the '{ext.name}' C extension, cDNA_Cupcake "
                f"will use the pure python version of the extension."
            )


def build(setup_kwargs):
    """
    This function is mandatory in order to build the extensions.
    """
    distribution = Distribution(
        {
            "name": "src/cupcake",
            "ext_modules": cythonize(
                extensions, compiler_directives={"language_level": "3"}
            ),
        }
    )
    distribution.package_dir = "cupcake"

    cmd = ExtBuilder(distribution)
    cmd.ensure_finalized()
    cmd.run()

    # Copy built extensions back to the project
    for output in cmd.get_outputs():
        output = Path(output)
        relative_extension = Path("src").joinpath(output.relative_to(cmd.build_lib))
        if not output.exists():
            continue

        shutil.copyfile(output, relative_extension)
        mode = relative_extension.stat().st_mode
        mode |= (mode & 0o444) >> 2
        relative_extension.chmod(mode)

    return setup_kwargs


if __name__ == "__main__":
    build({})
