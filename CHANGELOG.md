# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [21.05.6] - 2021-04-15

### Changed
- Removed pbcore, pbcoretools, and pbcommand from the list of dependencies as
  they were not being used.
  
## [21.0.5] - 2021-04-07

### Fixed
- Hot fix for importlib.metadata.metadata usage

## [21.0.4] - 2021-04-07

### Changed
- Replaced all instances of Pathobj.open() with open(Pathobj) (to catch
  instances where a string was being passed by an outside-cupcake module)
- All instances of `__version__`, `__author__`, and `__email__` now use the
  entries in `cupcake.__about__` as their single source, and `cupcake.__about__`
  uses `importlib.metadata` (or `import_metadata`) to extract that information
- Update dependencies

### Fixed
- minor typo fixes

## [21.0.3] - 2021-03-31

### Changed
- Replaced the inactive `PyVCF` module with the actively developed `vcfpy`

## [21.0.2] - 2021-03-26

### Fixed
- Fixed bugs in cupcake.sequence.GFF and cupcake.tofu.branch.branch_simple2 that
  were causing collapse_isoforms_by_sam to fail

## [21.0.1] - 2021-03-25

### Fixed
- Changed the "scripts" in pyproject.toml to call `app` instead of `main`
  since no scripts with typer appear to work otherwise
- Cleanup of `typer.Argument` objects that were given alternative parameter names
- Stray bugs picked up by flake8

## [21.0.0] - 2021-03-24

### Changed
- Completely replaced argparse usage with typer
- Completely replaced printing to stderr with a logger
- Partial replacement of using the os module with pathlib
## [19.1.2] - 2021-03-22

### Changed
- submodules in cupcake.tofu.counting, cupcake.tofu.branch, cupcake.singlecell, cupcake.annotation:
  - Function args and return have typing values
  - context managers and pathlib.Path for file handling
  - linting


## [19.1.1] - 2021-03-01

### Changed
- Most changes in cupcake.tofu.counting.chain_samples and 
  cupcake.tofu.counting.chain_fusion_samples
  - Replaced argparse with typer
  - os.path with pathlib.Path
  - format strings with f-strings
  - Typing

## [19.1.0] - 2021-02-26

### Added
- CHANGELOG.md

### Changed
- Started total replacement of `argparse` in favor of `typer`
  - So far, only completed the root `cupcake.tofu` module
- Further replacements of older .format style strings and string concatination in favor of f-strings
- Attempted to mostly replace instances of `chr` and `type` that conflicting with Python builtins with
  the more correct (according to the GFF/GTF definitions) `seqname` and `feature`, respectively
- Replaced instances of `print(..., file=sys.stderr)` with the use of a logger

### Fixed
- Miscellaneous linting problems
- removed unnecessary list comprehensions
- replaced calls to `dict` with dictionary comprehensions
- convert strings for regex to raw strings

## [19.0.2] - 2021-02-20

### Changed
- Reformmated to follow PEP517 - that is, replaced setup.py with pyproject.toml and build.py
- Now using Poetry for dependency managment


[21.0.5]: https://github.com/milescsmith/cDNA_Cupcake/compare/21.0.4...21.0.5
[21.0.4]: https://github.com/milescsmith/cDNA_Cupcake/compare/21.0.3...21.0.4
[21.0.3]: https://github.com/milescsmith/cDNA_Cupcake/compare/21.0.2...21.0.3
[21.0.2]: https://github.com/milescsmith/cDNA_Cupcake/compare/21.0.1...21.0.2
[21.0.1]: https://github.com/milescsmith/cDNA_Cupcake/compare/21.0.0...21.0.1
[21.0.0]: https://github.com/milescsmith/cDNA_Cupcake/compare/19.1.2...21.0.0
[19.1.2]: https://github.com/milescsmith/cDNA_Cupcake/compare/19.1.1...19.1.2
[19.1.1]: https://github.com/milescsmith/cDNA_Cupcake/compare/19.1.0...19.1.1
[19.1.0]: https://github.com/milescsmith/cDNA_Cupcake/compare/19.0.2...19.1.0
[19.0.2]: https://github.com/milescsmith/cDNA_Cupcake/compare/19.0.1...19.0.2
