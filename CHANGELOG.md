# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- added --geno option to set filtering criteria of SNPs in calculation of PFB file,
  equivalent of the PLINK --geno flag.
### Fixed

## [0.0.2] - 2022-10-12

### Added 

- Now PFB calculation is approx the same as PENN-CNV perl scripts,
  but with(hopefully) more performant polars code. 

## [0.0.1] - 2022-10-10

### Added

- Initial nextflow dsl2 pipeline to generate CNV calls
- Takes Illumina Genome Studio report files with Gtype, 
  B allele frequency, and LRR cols
- Implements PENN CNV pipline until generation of raw CNV calls,
  but with improved paralellisation.
  
