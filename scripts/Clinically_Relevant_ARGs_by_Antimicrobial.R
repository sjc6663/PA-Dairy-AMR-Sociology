#MEGARes v3.0 Extracting Unique Genes Based on WHO 2018 Clinically Important Antimicrobials
#Samantha Seibel 20230620

#Load packages
library(dplyr)
library(stringr)
library(tidylog)
library(tidyr)

#Read MEGARes v3 only file
MEGAResv3only <- read.table("C:\\Users/saman/OneDrive - The Pennsylvania State University/Ganda Lab/Spring 2023/Thesis Research/rhAMR/ARGs/Clinically_Relevant_ARGs/Data/20230601_annotations/megv3_only")

#Fix the hyphens to underscores
MEGAResv3only$Class <- gsub('Multi-', 'Multi_', MEGAResv3only$Class)

#Remove Unconfirmed SNPs
#MEGAResv3only_confirmed <- MEGAResv3only[is.na(MEGAResv3only$SNPConfirmation), ]

#Filter by drug class according to WHO Clinically Relevant Antimicrobials

Aminoglycosides <- MEGAResv3only %>%
  filter(Class == "Aminoglycosides")

BetaLactams <- MEGAResv3only %>%
  filter(Class == "betalactams")

Fluoroquinolones <- MEGAResv3only %>%
  filter(Class == "Fluoroquinolones")

Fosfomycins <- MEGAResv3only %>%
  filter(Class == "Fosfomycin")

Glycopeptides <- MEGAResv3only %>%
  filter(Class == "Glycopeptides")

Lipopeptides <- MEGAResv3only%>%
  filter(Class == "Lipopeptides")

Macrolide_Lincosamide_Streptogramins <- MEGAResv3only %>%
  filter(Class == "MLS")

Rifampin <-  MEGAResv3only %>%
  filter(Class == "Rifampin")

Oxazolidinones <-  MEGAResv3only %>%
  filter(Class == "Oxazolidinone")

Multi_drug_resistance <- MEGAResv3only %>%
  filter(Class == "Multi_drug_resistance")

Drug_and_biocide_resistance <- MEGAResv3only %>%
  filter(Class == "Drug_and_biocide_resistance")

Multi_biocide_resistance <- MEGAResv3only %>%
  filter(Class == "Multi_biocide_resistance")


#Remove duplicated rows based on Genes
Aminoglycosides_unique <- Aminoglycosides %>% 
  distinct(Gene, .keep_all = TRUE)

BetaLactams_unique <- BetaLactams %>% 
  distinct(Gene, .keep_all = TRUE)

Fluoroquinolones_unique <- Fluoroquinolones %>% 
  distinct(Gene, .keep_all = TRUE)

Fosfomycins_unique <- Fosfomycins %>% 
  distinct(Gene, .keep_all = TRUE)

Glycopeptides_unique <- Glycopeptides %>% 
  distinct(Gene, .keep_all = TRUE)

Lipopeptides_unique <- Lipopeptides %>% 
  distinct(Gene, .keep_all = TRUE)

MLS_unique <- Macrolide_Lincosamide_Streptogramins %>% 
  distinct(Gene, .keep_all = TRUE)

Rifampin_unique <- Rifampin %>% 
  distinct(Gene, .keep_all = TRUE)

Oxazolidinones_unique <- Oxazolidinones %>%
  distinct(Gene, .keep_all = TRUE)

Multi_drug_resistance_unique <- Multi_drug_resistance %>%
  distinct(Gene, .keep_all = TRUE)

Drug_and_biocide_resistance_unique <- Drug_and_biocide_resistance %>%
  distinct(Gene, .keep_all = TRUE)

Multi_biocide_resistance_unique <- Multi_biocide_resistance %>%
  distinct(Gene, .keep_all = TRUE)

#Combine into one file of all unique genes associated with clinically important Antimicrobials
all_imp_genes <- rbind(Aminoglycosides_unique, BetaLactams_unique, Fluoroquinolones_unique, Fosfomycins_unique, Glycopeptides_unique, Lipopeptides_unique, MLS_unique, Oxazolidinones_unique, Rifampin_unique, Multi_biocide_resistance_unique, Drug_and_biocide_resistance_unique, Multi_drug_resistance_unique)

View(all_imp_genes)

#Create a text file for all the clinically important antimicrobial associated genes
write.table(all_imp_genes, file = "C:\\Users/saman/OneDrive - The Pennsylvania State University/Ganda Lab/Spring 2023/Thesis Research/rhAMR/ARGs/Clinically_Relevant_ARGs/Results/CIA_genes", sep = "\t") 
