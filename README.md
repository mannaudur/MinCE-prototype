# MinCE
## Minhash with Cluster Extensions

MinCE is a fast, memory efficient tool to identify species and strains of unknown sequencing data. It accepts fasta and fastq files and even large, unprocessed fastq files with several unknown genomes. It is based on Mash (https://github.com/marbl/Mash), which offers quick comparison between genomes, and extends its resolution to the point of being able to distinguish between strains within species. This is achieved by pre-processing the dataset on which it operates and as such, MinCE is dataset-specific.

A pre-processed dataset of 285.339 genomes (254.025 Eubacteria and 4313 Archaea) is available at https://www.dropbox.com/s/uzxit07gt9bnhjy/MinCE_GTDB202_285K.zip?dl=0, though be advised that this is already present in the current Docker image. The dataset is based on the Genome Taxonomy Database 202 Release, but a few files are absent. For a description of the pre-processed dataset's content and our compilation of it, turn to https://www.dropbox.com/s/1rndrglvglmp10m/MinCE%20databank%20info.zip?dl=0. 

A dataset is necessary to run MinCE and as of this moment, this is the only pre-processed dataset available. In the setup/ folder are basic instructions on how to create other datasets and you're welcome to contact us with any questions in that regard.

The dataset of 285K genomes is roughly 3 Gb and, when queried with an unkown genome, will return the GenBank (GCA_) og RefSeq (GCF_) assembly accession string for the given genome, if it exists in the dataset. This identifier translates to NCBI searches. If the exact genome is not present in the dataset, MinCE will return its closest relative within the set. By default, the top 5 results are given. 

For single-genome files, MinCE gives results within 5 or 10 minutes, depending on your hardware. The easiest way to set up MinCE is to use Docker. 