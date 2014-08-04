KineticWavelets
===============

Prerequisites
==============

Bioconductor packages
---------------------

Biostrings and ShortRead packages from Bioconductor.

If you are unfamiliar with bioconductor merely run the following
commands at your R prompt.

`source("http://bioconductor.org/biocLite.R")`
`biocLite(c('shortRead','Biostrings'))`

h5r
---

Install hdf5 library for your OS, for me using OSX and the homebrew package manager
the relevant command was `brew install hdf5` 

Run the command.

`wget http://cran.r-project.org/src/contrib/Archive/h5r/h5r_1.4.7.tar.gz`

Then install from within R using install.packages setting repos to NULL and
the type to 'source'.

`install.packages('<path to h5r.tar.gz>', repos = NULL, type="source")`

Any packages available from CRAN will be installed by default. 

R-pbh5
------

Run the command

`git clone https://github.com/PacificBiosciences/R-pbh5`

then

`R CMD INSTALL R-pbh5`


Installation
============

In a R prompt type

`install.packages('devtools')`

`library(devtools)`

`install_github('KineticWavelets',username='sterlo')`

to use the package type.

`library(KineticWavelets)`

Tutorial
========

Download a reference fasta file and some CMPh5 data. 

Create a new KineticWavelets object

`kin = KineticWavelets(h5="h5 path",reff="reference path")1

Example wave128 window run.

`wav = wave128Window(kin,DNAPattern="GGCGGC")`
`plotSmoothAverage(wav)`
`plotDetailWave(wav)`

Example baseCorrelation run.

`bas = baseCorrelation(kin,DNAPattern="GGCGGC)`
`plotPatternCorrelation(bas)`
`plotBaseCorrelation(bas)`


