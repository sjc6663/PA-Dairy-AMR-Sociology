# Analysis of Whole Genome Sequencing (Shotgun Metagenomics) Data

This data is from Pennsylvania Dairy Farms and is related to the antimicrobial resistance genes found on those farms and thier relationship to sociological factors associated with the farm and farm operator/manager. The results are presented in a manuscript submitted to Journal of Rural Studies:

Ransom, E., Bierly, S., Ganda, E., Exploring the Relevance of Gender in Agriculture for Antimicrobial Resistance: A Case Study of PA Dairy Farms. Journal of Rural Studies. Submitted for Review January 2024. 

This repository includes scripts related to analyzing that data. Analysis was done using AMR++ for gene alignment and downstream analysis was completed with R in RStudio. Below is a breakdown of what each file is relevant for related to analysis. 


## Data

### CIA_genes.csv
This file lists all of the clinically signficant genes selected based on the WHO list of clinically significant antimicrobial classes. 

### PADairy-AMR-Metadata.xlsx
This is the metadata file with anonymous survey results for the study. 

### countmatrix-cleanedall.txt
This is the count table for gene counts for all samples. This is used to fill the "otu_table" slot in the phyloseq object.

### geneinfo-all.txt
This is the gene table for all gene information for all samples. This is used to fill the "tax_table" slot in the phyloseq object. 

### decontam-ps.rds
This is the decontaminated phyloseq object. 

### sig-decontam-ps.rds
This is the decontaminated phyloseq object only containing information for the clinically significant genes. 

## Scripts

### bash scripts folder
These bash scripts are used on the raw data to run the samples through AMRPlusPlus in Penn State's Roar Collab Super Computer. Slight modification can be performed to fit the job to any linux based system. 

### merge-outputs.R
This script is used to generate the batches of output data from AMRPlusPlus into one matrix gene file and one matrix count file to import into phyloseq.

### generate-phyloseq.R
This script is used to take the count matrix, gene info matrix, and metadata to generate a single phyloseq object that can be used to manipulate data for downstream analysis. 

### decontam.R
This script utilizes the decontam R package to remove contaminant genes based on the Positive and Negative Control samples sequenced and put through the AMRPlusPlus pipeline.

### sig-ps-creation.R
This script is used to combine the Clinically Significant Genes csv file (CIA_genes.csv) with the phyloseq object to only select genes that are deemed clinically significant for downstream analysis. 

### figure2.R
This script is used for the generation of bar plots for cows and calves separating the samples by farm and Broadclass of AMR genes. 

### figure3.R
This script is used for the generation of a bar plot, alpha diversity statistical test, violin plot, beta diversity statistical test, and PCA of data based on gender of primary farm operator. 

### figure4.R 
This script is used for the generation of a bar plot, alpha diversity statistical test, violin plot, beta diversity statistical test, and PCA of data based on cow vs. calf identification of samples. 

### figure5.R
This script is used for the generation of a bar plot, alpha diversity statistical test, violin plot, and beta diversity statistical test of data based on type of labor at the farm (family vs. non-family). 

### figure6.R
This script is used for the generation of a bar plot, alpha diversity statistical test, violin plot, and beta diversity statistical test of data based on whether or not farms experience a language barrier between primary operator and employees. 

### figure7.R
This script is used for the generation of a bar plot, alpha diversity statistical test, violin plot, beta diversity statistical test, and PCA of data based on gender of primary farm operator for clincally significant genes. 

