# methylCRF 1.1

## 1.1 update: 
* Replace olapBed with mapBed in bedtools
* Fix a compiling error of crfasgd on mac 

## about this tool:
* Estimate genome wide methylation level at each CpG sites using MeDIP-seq and MRE-seq data. Please read [methylcrf website](methylcrf.wustl.edu) for more information.
* Publication: [Stevens_GenomeResearch2013](http://genome.cshlp.org/content/23/9/1541.abstract)

## How to install:

* step1: open Makefile with your favorite text editor, and replace "/usr/local/bin/mapBed" in line 29 with the path to mapBed on your machine.
* step2: make

Enjoy!
