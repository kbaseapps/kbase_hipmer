# Hipmer release notes
=========================================

2.3.0
-----
* reactivated App
* updated slurm submit specs for Perlmutter (from Cori)
* added an authorized user list which can exceed the input size limit of 500 Gbp

2.2.4
-----
* hide App

2.2.3
-----
* remove deprecated ContigSet output type
* tidy up App Docs

2.2.2
-----
* use release AssemblyUtil

2.2.1
-----
* updated icon
* added arg to limit input size by Gbp

2.2.0
-----
* added mer size defaults of 21,33,55,77,99
* added scaff mer sizes with default of 99,33

2.1.0
-----
* Update to MetaHipMer v2.1.0.1-357
* dropping the 127 kmer arg
* use PrgEnv-cray/6.0.10
* use upcxx/2022.3.0

2.0.0
-----
* Update to MetaHipMer v2.0.2a

1.2.1
-----
* Split out metagenome and plant assembly into separate apps

1.2.0
-----
* Switch to new HPC model for execution
* Update to 1.2.2-7

1.1.1
-----
* Made input mandatory

1.1.0
-----
* Added parameter to filter short length contigs from assembly

1.0.7
-----
* Update to 1.2.1.48
* Support for new runner

1.0.6
-----
* Removed single-end from allowed types

0.1.0
-----
* Update to new SDK and Python3
* Add debug option
