# Start R


# first install superFreq if running first time
# load devtools that allows installation from github (may need to install devtools first with install.packages("devtools"))
library(devtools)

#there are sometimes conflicts between the github install and bioconductor in different version
#so safer to manually install bioconductor dependencies.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("GenomeInfoDb")
BiocManager::install("GenomicFeatures")
BiocManager::install("VariantAnnotation")

#then install superFreq
install_github('ChristofferFlensburg/superFreq')

#load and test superFreq
library(superFreq)
?superFreq

library(superFreq)

#maximum number of threads. Limited speed up above ~5 cpus for exomes and RNA-Seq and ~10-20 for genomes.
#Better to parallelise across individuals for cohorts, see the cohort section in the github README.
cpus=1

#this is the meta data input. See ?superFreq for how to set it up.
metaDataFile = '/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/Superfreq/superfreq.tsv'

#this is the capture region for bed file input. 
captureRegions = '/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/reference/fixed.bed'

#This directory with (links to) the reference normals needs to be created and set up. See ?superFreq
normalDirectory = '/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/Superfreq/normal'

#The reference fasta and name. Only hg19, hg38 and mm10 available atm.
reference = '/bgfs/alee/LO_LAB/Personal/Daisong/Reference_genome/Mouse/GRCm38.p4.genome/GRCm38.p4.genome.fasta'
genome = 'mm10'

#The directory where the log file and saved .Rdata is stored. Will be created.
Rdirectory = '/ihome/alee/dal413/R/Superfreq'
#The directory where all the plots and tables from the analysis go. Will be created.
plotDirectory = '/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/Superfreq/plot'

#The mode. Default 'exome' is for exomes, while 'RNA' has some minor changes when running on RNA.
#There is also a "genome" mode for genomes: ~24h for cancer-normal at 10 cpus, 200GB memory.
mode = 'exome'

#this performs the actual analysis. output goes to Rdirectory and plotDirectory.
#runtime is typically less than 6 hours at 4 cpus for a cancer-normal exome, but can vary significantly depending on input.
#For a typical cancer-normal exome, 5-10GB of memory is used per cpus, but again, can vary significantly depending on input.
#later runs typically a bit faster as the setup and part of the analysis on the reference normals can be reused.
data =
  superFreq(metaDataFile, normalDirectory=normalDirectory,
            Rdirectory=Rdirectory, plotDirectory=plotDirectory, reference=reference, genome=genome,
            cpus=cpus, mode=mode)