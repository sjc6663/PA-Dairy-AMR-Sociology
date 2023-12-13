# Create Phyloseq Object of Taxa
# SAB 12/13/2023

# load packages
library(phyloseq)
library(microViz)
library(stringr)
library(readxl)

#read in biom file
ps <- import_biom("taxonomy/bracken.biom")

#fix sample names
sample_names(ps) <- str_remove(sample_names(ps), ".non.host") #drop non.host extension

#fix tax names
colnames(tax_table(ps)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

meta <- read_excel("/Users/stephanieclouser/OneDrive - The Pennsylvania State University/Shared-Projects/Ransom-AMR/3.Sample-Collection/Ransom-AMR-Metadata.xlsx")

samp <- phyloseq::sample_data(meta)
rownames(samp) <- meta$`Sample-ID`

#fix meta to match AMR ++
sample_data(ps) # has the file extenstion stuff so fix
ps <- ps %>%
  ps_mutate(Sample.ID = rownames(sample_data(ps))) #add column to join
ps <- ps %>%
  ps_select(-c(Id)) #drop "Id" column with file extension stuff

sample_data(ps) <- samp

# in the tax table all the taxa have a letter __ in the names, let's remove that
tax_table(ps)

# find and substitute
tax_table(ps) <- gsub(tax_table(ps), pattern = "k__", replacement = "") %>% 
  gsub(tax_table(ps), pattern = "p__", replacement = "") %>% 
  gsub(tax_table(ps), pattern = "c__", replacement = "") %>% 
  gsub(tax_table(ps), pattern = "o__", replacement = "") %>% 
  gsub(tax_table(ps), pattern = "f__", replacement = "") %>% 
  gsub(tax_table(ps), pattern = "g__", replacement = "") %>% 
  gsub(tax_table(ps), pattern = "s__", replacement = "")

tax_table(ps)

save(ps, file = "taxonomy/R/ps_tax.RData") #save
