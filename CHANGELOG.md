# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),

## [21.7.0] - 2021-07-14

### Fixed
- Fix the final reporting of the output location by
  `cupcake.annotation.make_file_for_subsampling_from_collapsed.make_file_for_subsample()`
- Logging now works as originally intended
- Fixed several instances of missing or extra parentheses being present in the
  arguments for `main()`

### Added
- All `main()` functions have a `--version` argument, all of which call a
  common `version_callback()` in `__init__.py`

### Changed
- Renamed "logging.py" submodule "logger.py"

## [21.6.0] - 2021-07-13

### Fixed
- Actually create the `good_ids` variable before it is referenced in
  `cupcake.annotation.make_file_for_subsampling_from_collapsed.make_file_for_subsample()`

## [21.5.0] - 2021-06-29

Version bump for Docker container

## [21.4.0] - 2021-06-24

### Changed
- `cupcake.tofu.filter_away_subset.sanity_check_collapse_input()` now just checks
  that the files are where the user says they are but still verifies their
  records all match
- `cupcake.tofu.filter_away_subset.main()` now requires the path
  to the needed files be explicitly passed instead of inferring by use of a prefix
  (needed for Nextflow pipeline)

### Fixed
- FASTA files written by - `cupcake.tofu.filter_away_subset.main()` now have
  the correct ".fa" suffix

### Added
- pre-commit hooks
- type checking/linting with pyre


## [21.3.0] - 2021-06-23

### Changed
- `cupcake.tofu.get_abundance_post_collapse.get_abundance_post_collapse()` now takes the path to the group_file created by `collapse_isoforms_by_sam()` instead of a prefix (makes it better compatible with pipelines)
- Removed `count_header` from the return values of
  `cupcake.tofu.counting.chain_samples.read_count_info()` as it was not being used
- Update dependencies

## [21.2.8] - 2021-06-11

### Changed
- Stopped err_correct_w_genome() from printing every sequence it finds to the log

## [21.2.7] - ?

## [21.2.6] - 2021-05-18

### Fixed
- Fixed Python 2 `print` stragglers in cython code

## [21.2.5] - 2021-05-18

### Fixed
- Updated Dockerfile to match the working one in [nf-isoseq](milescsmith/nf-isoseq)

## [21.2.4] - 2021-05-17

### Changed
- Revised logging and now incorporating the coloredlogs module

### Fixed
- Removed most instances of `Path.write_text()`, replacing it with a standard
  open file and write using a context manager.

## [21.2.3] - 2021-05-10

### Changed
- Removed all dev dependencies because I cannot figure out how to get a module that
  relies on cupcake to ignore them on installation (ipython was killing my
  SQANTI3 docker build times)

## [21.2.2] - 2021-05-10

### Fixed
- Fixed a problem with opening gzipped FASTA/FASTQ files in tofu/collapse_isoforms_by_sam

### Fixed
- Reverted a `[]` to `list()` change in branch_simple2.py

## [21.2.1] - 2021-05-10

### Fixed
- Reverted a `[]` to `list()` change in branch_simple2.py

## [21.2.1] - 2021-05-07

### Fixed
- Cython did not like some of the c_branch.pyx code.

## [21.2.0] - 2021-05-07

### Added
- Added `OpenFile` class, an enhanced version of the standard `open` function
  that (should) transparently handle (b)gzipped files.  Currently only used
  when reading FASTA/FASTQ files.

## [21.1.1] - 2021-04-15

### Fixed
- Hotfix to allow `fa2fq` and `fq2fa` to correctly accept files with an
  `.fa` and `.fq` extension


## [21.1.0] - 2021-04-15

### Added
- added `fa2fq`, `fq2fa`, `get_seq_stats`, `rev_comp`, `get_seqs_from_list`,
  `err_correct_w_genome`, and `group_ORF_sequences` to the list of scripts

## [21.0.6] - 2021-04-15

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

[21.2.8]: https://github.com/milescsmith/cDNA_Cupcake/compare/21.6.0...21.7.0
[21.2.8]: https://github.com/milescsmith/cDNA_Cupcake/compare/21.5.0...21.6.0
[21.2.8]: https://github.com/milescsmith/cDNA_Cupcake/compare/21.4.0...21.5.0
[21.2.8]: https://github.com/milescsmith/cDNA_Cupcake/compare/21.3.0...21.4.0
[21.2.8]: https://github.com/milescsmith/cDNA_Cupcake/compare/21.2.8...21.3.0
[21.2.8]: https://github.com/milescsmith/cDNA_Cupcake/compare/21.2.7...21.2.8
[21.2.7]: https://github.com/milescsmith/cDNA_Cupcake/compare/21.2.6...21.2.7
[21.2.6]: https://github.com/milescsmith/cDNA_Cupcake/compare/21.2.5...21.2.6
[21.2.5]: https://github.com/milescsmith/cDNA_Cupcake/compare/21.2.4...21.2.5
[21.2.4]: https://github.com/milescsmith/cDNA_Cupcake/compare/21.2.3...21.2.4
[21.2.3]: https://github.com/milescsmith/cDNA_Cupcake/compare/21.2.2...21.2.3
[21.2.2]: https://github.com/milescsmith/cDNA_Cupcake/compare/21.2.1...21.2.2
[21.2.1]: https://github.com/milescsmith/cDNA_Cupcake/compare/21.2.0...21.2.1
[21.2.0]: https://github.com/milescsmith/cDNA_Cupcake/compare/21.1.1...21.2.0
[21.1.1]: https://github.com/milescsmith/cDNA_Cupcake/compare/21.1.0...21.1.1
[21.1.0]: https://github.com/milescsmith/cDNA_Cupcake/compare/21.0.5...21.1.0
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
