# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://github.com/olivierlacan/keep-a-changelog/compare/19.0.2...19.1.0
[Unreleased]: https://github.com/olivierlacan/keep-a-changelog/compare/19.0.1...19.0.2