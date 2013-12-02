KineticWavelets
===============

Prerequisites
==============

Any packages available from CRAN will be installed by default. 

Other packages that are needed are the R-pbh5 library from Pacific Biosciences.

https://github.com/PacificBiosciences/R-pbh5

And the Biostrings and ShortRead packages from Bioconductor.

Installation
============

In a R prompt type

install.packages('devtools')

library(devtools)

install_github('KineticWavelets',username='sterlo')

to use the package type.

library(KineticWavelets)

Tutorial
========

Download a reference fasta file and some CMPh5 data. 

Create a new KineticWavelets object


kin = KineticWavelets(h5="h5 path",reff="reference path")

Example wave128 window run.

wav = wave128Window(kin,DNAPattern="GGCGGC")
plotSmoothAverage(wav)
plotDetailWave(wav)

Example baseCorrelation run.

bas = baseCorrelation(kin,DNAPattern="GGCGGC)
plotPatternCorrelation(bas)
plotBaseCorrelation(bas)


