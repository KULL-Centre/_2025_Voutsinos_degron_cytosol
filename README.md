# _2025_Voutsinos_degron_cytosol
Scripts and output from "A complete map of human cytosolic degrons and characterization of their exposure and relevance for disease"
by Vasileios Voutsinos, Kristoffer E. Johansson, Fia B. Larsen, Martin Grønbæk-Thygesen, Nicolas Jonsson, Giulio Tesei, Amelie Stein, Douglas M. Fowler, Kresten Lindorff-Larsen and Rasmus Hartmann-Petersen.

Contents
--------
- **library**: Script and input files for makeing the DNA libraries
- **counts**: Scripts for processing FASTQ files. Output counts and raw FASTQ files are available on ERDA
- **score**: Scripts for calculating scores. FACS data files are available on ERDA
- **models**: Scripts for training models
- **pathogenic**: Scripts for analysing pathogenic variants from ClinVar
- **proteome**: Scripts for building the "cytosolic proteome" and structural analysis 
- **plots**: Scripts for plotting 

Sequencing counts, FASTQ files and FACS data are available on [ERDA](https://sid.erda.dk/sharelink/FUolgBZn9W)

Peptide abundance predictor (PAP)
---------------------------------

Predictive models described in the paper are available in the script models/PAP.py. The neural network model requires a weight file compressed in models/pap_weights.tgz and depends on Keras2 available in tensorflow version 2.14.

PAP webserver
-------------

The file PapLab.ipynb is made to run as a webservice at Google colab available at
https://colab.research.google.com/github/KULL-Centre/_2025_Voutsinos_degron_cytosol/blob/main/PapLab.ipynb

This requires a login for Googles services, e.g. a gmail.
