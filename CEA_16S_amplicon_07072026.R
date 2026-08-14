#SAVING THE ENVIRONMENT IN R
#https://www.techcoil.com/blog/how-to-save-and-load-environment-objects-in-r/
save.image(file='/local1/workdir1/tw488/CEA_RData_RScript/CEA_check_07072026.RData')
#VIEW THE SAVED FILE IN A THE DIRECTORY
#dir()
#TERMINATE R SESSION
#quit(save='no')
#Loading back the entire list environment of objects
load('/local1/workdir1/tw488/CEA_RData_RScript/CEA_check_07072026.RData')

#Load packages
#https://benjjneb.github.io/dada2/tutorial.html
#https://benjjneb.github.io/dada2/ITS_workflow.html
#if (!requireNamespace("BiocManager", quietly = TRUE))y
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Maaslin2")
BiocManager::install("dada2")
install.packages('theseus')
library("theseus")
library(dada2); packageVersion("dada2")
library(tidyverse)
library(ggplot2)
library(ShortRead)
library(Maaslin2)
library(phyloseq)
library(tools)
library(vegan)
library(reshape2)

#This workflow is from  https://benjjneb.github.io/dada2/ITS_workflow.html
#Map the path to the sequencing files on my drive
path <- "/local1/workdir1/tw488/CEA_16S_amplicon/For_Analysis_v2/"
list.files(path)

#Generate matched lists of the forward and reverse read files
fnFs <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq.gz", full.names = TRUE))

list(fnFs)
list(fnRs)

#Identify the primers used to make these amplicons
FWD <- "GTGCCAGCMGCCGCGGTAA"  ## 
REV <- "GGACTACHVGGGTWTCTAAT"  ##
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

#remove poly Ns from reads, create output filtN
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

#Heads up- this step takes a long time
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

#counting primers in F and R reads
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#the above step took awhile too
#What are the results? Are all F reads in F position, and all R reads in R position?

#remove primers with cutadapt- 
cutadapt <-"/workdir1/tw488/CIDA_Sequencing_Runs/cutadapt-venv/bin/cutadapt"
system2(cutadapt, args = "--version")

#create output filenames for the cutadapt-ed files,
path.cut <- file.path(path, "cutadapt_CEAcheck_07072026")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads 
                             "-m", 1, "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

#double check for primers, output should all be 0,

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))


# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq.gz", full.names = TRUE))
list(cutFs)

q<-basename(file_path_sans_ext(cutFs))
#this only got rid of .gz, need to get rid of .fastq as well
sample.names <-(file_path_sans_ext(q))
#.fastq.gz should be removed
head(sample.names)
list(sample.names)
type(sample.names)



#inspect quality reads
plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])

#Start following this workflow: https://benjjneb.github.io/dada2/tutorial.html
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path.cut, "filtered_CEAcheck07072026", basename(cutFs))
filtRs <- file.path(path.cut, "filtered_CEAcheck07072026", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs,  maxN = 0, maxEE = c(2, 2),truncLen=c(230,230), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)

list(out)
plotQualityProfile(filtFs[1:2])
plotQualityProfile(filtRs[1:2])
class(out)

list(filtFs)
class(filtFs)
list(filtRs)
class(filtRs)
list(sample.names)

#Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


#run DADA
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)


#Merge paired reads
#Most of your reads should successfully merge. If that is not the case upstream parameters may need to be revisited: Did you trim away the overlap between your reads?
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)


#Track reads through the pipeline- tracks how many reads were filtered/discarded at what step
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
list(track)


#Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/local1/workdir1/tw488/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


## 1. ASV Table 
# Prep the ASV table! 
samples_out <- rownames(seqtab.nochim)
head(samples_out)
# Pull out sample names from the fastq file name 
sample_names_reformatted <-gsub("_R1.fastq.gz","",samples_out)
head(sample_names_reformatted)


# Replace the names in our seqtable 
#Check the names and sort
seqtab.nochimdf<-as.data.frame(seqtab.nochim)
sample_names_reformatted.df<-as.data.frame(sample_names_reformatted)
#sort
#replace
rownames(seqtab.nochim) <- sample_names_reformatted
seqtab.nochimdfv2<-as.data.frame(seqtab.nochim)
head(seqtab.nochim)

#Do the numbers match up with the new names? 
#If so then continue; if not then figure out where the sorting went wrong
### intuition check 
h<-stopifnot(rownames(seqtab.nochim) == sample_names_reformatted)
h

# Modify the ASV names and then save a fasta file! 
# Give headers more manageable names
# First pull the ASV sequences
asv_seqs <- colnames(seqtab.nochim)
# make headers for our ASV seq fasta file, which will be our ASV names
asv_headers <- vector(dim(seqtab.nochim)[2], mode = "character")
# loop through vector and fill it in with ASV names 
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")
}
# intuition check
asv_headers

# Rename ASVs in table then write out our ASV fasta file! 
#Look at the original files and the long ASV names (consisting of the DNA)
View(seqtab.nochim)
asv_tab <- t(seqtab.nochim)
View(asv_tab)
## Rename our asvs! 
row.names(asv_tab) <- sub(">", "", asv_headers)
View(asv_tab)

## 2. Taxonomy Table 
#prepare-tax-table
View(taxa)
# Add the ASV sequences from the rownames to a column 
new_tax_tab <- taxa %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASVseqs") 
head(new_tax_tab)

# intuition check 
h1<-stopifnot(new_tax_tab$ASVseqs == colnames(seqtab.nochim))
h1
list(asv_tab)
type(asv_tab)


# Now let's add the ASV names 
rownames(new_tax_tab) <- rownames(asv_tab)
View(new_tax_tab)
### Final prep of tax table. Add new column with ASV names 
asv_tax <- 
  new_tax_tab %>%
  # add rownames from count table for phyloseq handoff
  mutate(ASV = rownames(asv_tab)) %>%
  # Resort the columns with select
  dplyr::select(Kingdom, Phylum, Class, Order, Family, Genus, ASV, ASVseqs)
View(asv_tax)
# Intution check
j<-stopifnot(asv_tax$ASV == rownames(asv_tax), rownames(asv_tax) == rownames(asv_tab))
j


#Import the sample dataframe with sample metadata
samdf <-read.csv("/local1/workdir1/tw488/CEA_16S_amplicon/CEA_sampledata_forimport_07092025.csv", header=TRUE,row.names = 1)

#We now construct a phyloseq object directly from the dada2 outputs.
asv_tab<-as.matrix(asv_tab)
asv_tax<-as.matrix(asv_tax)
ps <- phyloseq(otu_table(asv_tab, taxa_are_rows=TRUE), 
               sample_data(samdf), 
               tax_table(as.matrix(asv_tax)))

sample_variables(ps)

ps

#code to remove ASVs==1, NC samples and chloroplast
ps.1 <- ps %>% subset_taxa( Family!= "mitochondria" | is.na(Family) & Class!="Chloroplast" | is.na(Class) )

ps.2 <- subset_samples(ps.1, SampleType!= "C")

#make ps object of just NC samples
ps.NC <- subset_samples(ps.1, SampleType== "C")

# remove taxa with 1 read count, from https://deneflab.github.io/Diversity_Productivity/analysis/OTU_Removal_Analysis.html
ps.3<-prune_taxa(taxa_sums(ps.2) > 1, ps.2) 

# check to see if above worked by summing rows and columns, as well as 
smdta3<-as.data.frame(ps.3@sam_data)
otutbl3<-as.data.frame(ps.3@otu_table)
txtbl3<-as.data.frame(ps.3@tax_table)
ps.3otu_colsums<-as.data.frame(colSums(otutbl3))
ps.3otu_rowsums<-as.data.frame(rowSums(otutbl3))

#Rarefy samples

ps.noncontam4 = rarefy_even_depth(ps.3, rngseed=1, sample.size=min(sample_sums(ps.3)), replace=TRUE, trimOTUs = TRUE, verbose = TRUE)
ps.noncontam4otu<-ps.noncontam4@otu_table
ps.noncontam4samdata<-as.data.frame(ps.noncontam4@sam_data)
ps.noncontam4otu_colsums<-as.data.frame(colSums(ps.noncontam4otu))
ps.noncontam4otu_rowsums<-as.data.frame(rowSums(ps.noncontam4otu))

class(ps.noncontam4otu) <- "matrix" 

## rarefaction curve
ps.noncontam4otu_flipped<-as.data.frame(t(ps.noncontam4otu))
ps.noncontam4otu.rarecurve = rarecurve(ps.noncontam4otu_flipped, step = 10, col = "blue", cex = 0.6, label = T)

#check distribution of OTU rowsums
otu_tbl_rowsums<-tibble::rownames_to_column(ps.noncontam4otu_rowsums, "ASV")
colnames(otu_tbl_rowsums)[2] <-"rowsums"
ggplot(otu_tbl_rowsums, aes(x=rowsums)) + geom_histogram()


#### Calculate relative abundance ####

#count_seqs function for use in calculating relative abundance
count_seqs <- function(pool){
  count_dat <- sample_sums(pool)
  names <- names(count_dat)
  values <- unname(count_dat)
  counts_1 <- cbind(names, values) %>% as.data.frame()
  colnames(counts_1) = c("sample", "counts") 
  
  counts_1 <- counts_1 %>%
    arrange(as.numeric(counts)) %>%
    mutate(counts = as.numeric(counts))
  counts_1$n <- seq.int(nrow(counts_1))
  
  plot <- ggplot(counts_1, aes(n, counts)) +
    geom_line()+
    theme_classic()
  
  return(list("counts"=counts_1, "plot"=plot))
  
}
ps.noncontam4.counts <- count_seqs(ps.noncontam4@otu_table)

abund.ps.noncontam4 <- ps.noncontam4@otu_table %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("ASV_") %>%
  tidyr::pivot_longer((-ASV_), names_to = "sample") %>%
  merge(ps.noncontam4.counts$counts, by="sample") %>%
  dplyr::select(-n) %>%
  mutate(newcol=sample)%>%
  tidyr::separate(newcol, c("Pond_Type","Pond_Number","Location"), "_",remove = TRUE)%>%
  group_by(sample) %>%
  mutate(rel_abund = value/counts) %>% #relative abundance by treatment/day/rep
  ungroup() 

#check that relative abundance is calculated correctly- everything should add up to '1'
abund.ps.noncontam4 %>% group_by(sample) %>% summarise(sum = sum(rel_abund))


#get taxonomy
ps.noncontam4tax <- ps.noncontam4@tax_table %>%
as.data.frame() %>%
 tibble::rownames_to_column("ASV_") 
#### NEW: ADJUST TAXONOMY LOGIC FOR DADA2 MATRIX IN PHYLOSEQ ####
# Row-by-row function to find the lowest valid higher rank if Genus is missing
#backtrack_taxa <- function(row) {
#ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family")

# genus_val <- row["Genus"]

# 1. If Genus is valid (not NA, not blank, and doesn't contain "unclassified"), return it as-is
#  if (!is.na(genus_val) && 
# stringr::str_trim(genus_val) != "" && 
#  !stringr::str_detect(stringr::str_to_lower(genus_val), "unclassified")) {
# return(genus_val)
#}

# 2. Backtrack from Family down to Phylum/Kingdom if Genus is missing/NA
#  for (r in rev(ranks)) {
#   current_val <- row[r]
#  
# if (!is.na(current_val) && 
#    stringr::str_trim(current_val) != "" && 
#   stringr::str_to_lower(current_val) != "incertae_sedis" && 
#  !stringr::str_detect(stringr::str_to_lower(current_val), "unclassified")) {

# We combine capitalized "Unclassified" with the original, un-mutated string value
#   return(paste("Unclassified", current_val))
#  }
# }

#  return("Unclassified Bacteria")
#}

# Apply the backtracking adjustment directly to overwrite your 'Genus' column
#ps.noncontam4tax$Genus <- apply(ps.noncontam4tax, 1, backtrack_taxa)

#head(ps.noncontam4tax)

#### CONTINUE ORIGINAL WORKFLOW ####

#merge taxonomy with rarefied relative abundance
count_tax.ps.noncontam4 <- merge(abund.ps.noncontam4,ps.noncontam4tax, by="ASV_")
head(count_tax.ps.noncontam4)

# Get summary abundance of each taxa
taxa.meta.sumps.noncontam4 <- count_tax.ps.noncontam4%>%
  dplyr::select(-c(value, counts)) %>%
  tidyr::pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "ASV_"),
                      names_to = "level", 
                      values_to = "taxon")

#format phylum names for Genus
genus_rel_abund_ps.noncontam4 <- taxa.meta.sumps.noncontam4 %>%
  filter(level == "Genus") %>%
  group_by(taxon,sample,Pond_Type,Pond_Number,Location) %>% 
  summarise(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  mutate(taxon = stringr::str_replace(taxon,
                                      "^unclassified (.*)", "Unclassified *\\1*"),
         taxon = stringr::str_replace(taxon, "^(\\S*)$", "\\1")) %>%
  ungroup()
genus_rel_abund_ps.noncontam4$taxon <-genus_rel_abund_ps.noncontam4$taxon%>%tidyr::replace_na('Unclassified')

#check how many taxa there are total
genus_rel_abund_ps.noncontam4_alltaxa <- genus_rel_abund_ps.noncontam4%>% filter(rel_abund!="0") %>%distinct(taxon)

#pool taxon <5% relative abundance into "<5% Relative Abundance" category
taxon_pool_ps.noncontam4 <- genus_rel_abund_ps.noncontam4 %>%
  group_by(sample, taxon, Pond_Type,Pond_Number,Location) %>%
  summarise(mean = mean(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarise(pool = max(mean) <5,
            mean=mean(mean), 
            .groups="drop") 

#arrange data frame to group samples by relative abundance within each taxon; 
#this will aid visual comparison between the barplots
join_taxon_rel_abund_ps.noncontam4 <- inner_join(genus_rel_abund_ps.noncontam4, taxon_pool_ps.noncontam4, by="taxon") %>%
  mutate(taxon = ifelse(pool, "<5% Relative Abundance", taxon)) %>%
  group_by(sample, taxon,Pond_Type,Pond_Number,Location) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  mutate(taxon = factor(taxon),
         taxon = forcats::fct_reorder(taxon, rel_abund, .desc=TRUE),
         taxon = forcats::fct_shift(taxon, n=1)) %>%
  ungroup() %>%
  unique() %>%
  group_by(taxon) %>%
  arrange(desc(rel_abund), .by_group=TRUE) 


#Make Tables of genera for organic and conventional samples in genus_rel_abund_ps.noncontam4

genera_org<-genus_rel_abund_ps.noncontam4%>%filter(Pond_Type=="Organic")%>%
  group_by(taxon)%>%
  mutate(Avg = mean(rel_abund))%>%distinct(taxon,Avg)%>%
  arrange(desc(Avg))%>%filter(Avg!="0")
colnames(genera_org)[1] <- "Organic Sample Genera"
colnames(genera_org)[2] <- "Average Relative Abundance (%)"
write.csv(genera_org,"/local1/workdir1/tw488/CEA_RData_RScript/manuscript_data_analysis/Organic_genera_relab182026", row.names = TRUE)

genera_convent<-genus_rel_abund_ps.noncontam4%>%
  filter(Pond_Type=="Conventional")%>%group_by(taxon)%>%
  mutate(Avg = mean(rel_abund))%>%distinct(taxon,Avg)%>%
  arrange(desc(Avg))%>%filter(Avg!="0")
colnames(genera_convent)[1] <- "Conventional Sample Genera"
colnames(genera_convent)[2] <- "Average Relative Abundance (%)"
write.csv(genera_convent,"/local1/workdir1/tw488/CEA_RData_RScript/manuscript_data_analysis/Conventional_genera_relab182026", row.names = TRUE)

#make taxa and sample data tables from phyloseq object
ps.noncontam4_sampledata<-data.frame(ps.noncontam4@sam_data)

#check to make sure all are actually data.frames
class(ps.noncontam4_sampledata)

#DCAST TO CONVERT THE DATA FRAME
#HERE
genus_rel_abund_ps.noncontam4_4maaslin<-dcast(genus_rel_abund_ps.noncontam4 ,sample~taxon,value.var="rel_abund")
genus_rel_abund_ps.noncontam4_4maaslin<-genus_rel_abund_ps.noncontam4_4maaslin%>% column_to_rownames(var="sample")

#run maaslin2 with relative abundance at the genus level

ps.noncontam4_genus_diffab = Maaslin2(input_data = genus_rel_abund_ps.noncontam4_4maaslin, 
                                      input_metadata = ps.noncontam4_sampledata,
                                      analysis_method = "LM",
                                      normalization  = "NONE", 
                                      transform="NONE",
                                      min_prevalence = 0,
                                      min_abundance = 0,
                                      save_models="TRUE",
                                      save_scatter="TRUE",
                                      random_effects = c("Pond_Type_Number"),
                                      fixed_effects  = "Pond_Type",
                                      reference      = c("Conventional"),
                                      output         = "/local1/workdir1/tw488/CEA_RData_RScript/Maaslin2/ps.noncontam4_genus_diffab_recheck07072026/",
                                      cores=2) 

#Import results in and filter by a <0.05 FDR
ps.noncontam4_genus_diffab_results<-read.table("/local1/workdir1/tw488/CEA_RData_RScript/Maaslin2/ps.noncontam4_genus_diffab_recheck07072026/all_results.tsv",header=TRUE)
ps.noncontam4_genus_diffab_results0.05<-ps.noncontam4_genus_diffab_results%>%filter(qval<0.05)%>%
  select("feature", "coef","qval")%>% mutate(log2_fold_change = log2(exp(coef)))%>%select("feature","qval","log2_fold_change")


genus_rel_abund_ps.noncontam4_TAXAmeanbyOvC<-genus_rel_abund_ps.noncontam4%>% group_by(taxon,Pond_Type) %>% summarise(AvgbyDay= mean(rel_abund))
genus_rel_abund_ps.noncontam4_TAXAmeanbyOvC1<-dcast(genus_rel_abund_ps.noncontam4_TAXAmeanbyOvC,taxon~Pond_Type,value.var="AvgbyDay")

colnames(ps.noncontam4_genus_diffab_results0.05)[1] <- "Genus"
colnames(genus_rel_abund_ps.noncontam4_TAXAmeanbyOvC1)[1] <- "Genus"
Genus_diffab_relab_table<-merge(genus_rel_abund_ps.noncontam4_TAXAmeanbyOvC1,ps.noncontam4_genus_diffab_results0.05,by="Genus")
colnames(Genus_diffab_relab_table)[2] <- "Conventional Sample Average Relative Abundance (%)"
colnames(Genus_diffab_relab_table)[3] <- "Organic Sample Average Relative Abundance (%)"
colnames(Genus_diffab_relab_table)[4] <- "Significance Level"
colnames(Genus_diffab_relab_table)[5] <- "Log2 Fold Change"
write.csv(Genus_diffab_relab_table,"/local1/workdir1/tw488/CEA_RData_RScript/manuscript_data_analysis/Genus_diffab_relab_table", row.names = TRUE)


#format phylum names for FAMILY________________________________________________________________
fam_rel_abund_ps.noncontam4 <- taxa.meta.sumps.noncontam4 %>%
  filter(level == "Family") %>%
  group_by(taxon,sample,Pond_Type,Pond_Number,Location) %>% 
  summarise(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  mutate(taxon = stringr::str_replace(taxon,
                                      "^unclassified (.*)", "Unclassified *\\1*"),
         taxon = stringr::str_replace(taxon, "^(\\S*)$", "\\1")) %>%
  ungroup()
fam_rel_abund_ps.noncontam4$taxon <-fam_rel_abund_ps.noncontam4$taxon%>%tidyr::replace_na('Unclassified')


#pool taxon <5% relative abundance into "<5% relative abundance" category
famtaxon_pool_ps.noncontam4 <- fam_rel_abund_ps.noncontam4 %>%
  group_by(sample, taxon, Pond_Type,Pond_Number,Location) %>%
  summarise(mean = mean(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarise(pool = max(mean) <5,
            mean=mean(mean), 
            .groups="drop") 

#arrange data frame to group samples by relative abundance within each taxon; 
#this will aid visual comparison between the barplots
join_famtaxon_rel_abund_ps.noncontam4 <- inner_join(fam_rel_abund_ps.noncontam4, famtaxon_pool_ps.noncontam4, by="taxon") %>%
  mutate(taxon = ifelse(pool, "<5% Relative Abundance", taxon)) %>%
  group_by(sample, taxon,Pond_Type,Pond_Number,Location) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  mutate(taxon = factor(taxon),
         taxon = forcats::fct_reorder(taxon, rel_abund, .desc=TRUE),
         taxon = forcats::fct_shift(taxon, n=1)) %>%
  ungroup() %>%
  unique() %>%
  group_by(taxon) %>%
  arrange(desc(rel_abund), .by_group=TRUE) 


color_map <- c( "<5% Relative Abundance"="black", "Chitinophagaceae" = "#665191",
                "Flavobacteriaceae"="#8DA0CB", "Micrococcaceae"= "#8F9E6F", "Devosiaceae"="#E9C46A", 
                "Pleomorphomonadaceae"= "#4D4D4D", "Lachnospiraceae"= "#E76F51","Mycobacteriaceae"="#8C564B",
                "Spirosomaceae"= "#3A86A8","env.OPS 17"="#AA4499","Sphingomonadaceae"= "#F95D6A",
                "Rhodobacteraceae"= "#009E73","Rhizobiaceae"="#E7298A","Unclassified"="#CDBE70",
                "Xanthobacteraceae" = "gray")

#### Relative abundance plot ####
levels(join_famtaxon_rel_abund_ps.noncontam4$taxon)
join_famtaxon_rel_abund_ps.noncontam4 %>%
  mutate(taxon = factor(taxon, levels = c(
    "<5% Relative Abundance","Chitinophagaceae","Devosiaceae","env.OPS 17",
    "Flavobacteriaceae","Lachnospiraceae",
    "Micrococcaceae","Mycobacteriaceae","Pleomorphomonadaceae",
    "Rhizobiaceae","Rhodobacteraceae",
    "Sphingomonadaceae","Spirosomaceae",
    "Xanthobacteraceae","Unclassified"))) %>%
  ggplot() +
  geom_col(mapping = aes(x = sample, y = rel_abund, fill = taxon), 
           color = "black", position = "fill", show.legend = TRUE, width = .4) +
  facet_grid(cols = vars(Pond_Type), scales = "free_x", space = "free_x") +
  ylab("Proportion of Community") +
  xlab("CEA Pond Water Samples") +
  scale_x_discrete(labels = c(
    "Conventional_1_H"="Conventional Pond 1 H", 
    "Conventional_1_M"="Conventional Pond 1 M", 
    "Conventional_1_S"="Conventional Pond 1 S",
    "Conventional_2_H"="Conventional Pond 2 H",
    "Conventional_2_M"="Conventional Pond 2 M", 
    "Conventional_2_S"= "Conventional Pond 2 S",
    "Organic_1_H"="Organic Pond 1 H", 
    "Organic_1_M"="Organic Pond 1 M", 
    "Organic_1_S"="Organic Pond 1 S",
    "Organic_2_H"="Organic Pond 2 H", 
    "Organic_2_M"="Organic Pond 2 M",
    "Organic_2_S"="Organic Pond 2 S" )) +
  theme_linedraw() +
  scale_fill_manual(
    values = color_map,
    labels = c(
      "<5% Relative Abundance" = "<5% Relative Abundance",
      "Unclassified" = "Unclassified",
      "env.OPS 17" = expression(italic("env.OPS 17")),
      "Lachnospiraceae"=expression(italic("Lachnospiraceae")),
      "Chitinophagaceae"=expression(italic("Chitinophagaceae")),
      "Devosiaceae"=expression(italic("Devosiaceae")),
      "Flavobacteriaceae"=expression(italic("Flavobacteriaceae")),
      "Micrococcaceae"=expression(italic("Micrococcaceae")),
      "Mycobacteriaceae"=expression(italic("Mycobacteriaceae")),
      "Pleomorphomonadaceae"=expression(italic("Pleomorphomonadaceae")),
      "Rhizobiaceae"=expression(italic("Rhizobiaceae")),
      "Rhodobacteraceae"=expression(italic("Rhodobacteraceae")),
      "Sphingomonadaceae"=expression(italic("Sphingomonadaceae")),
      "Spirosomaceae"=expression(italic("Spirosomaceae")),
      "Xanthobacteraceae"=expression(italic("Xanthobacteraceae")) )) +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.text.y = element_text(size = 30, color = "black"),
    strip.text = element_text(size = 50),
    axis.title.x = element_text(size = 50, color = "black", margin = margin(t = 30, r = 0, b = 30, l = 0)),
    axis.title.y = element_text(size = 50, color = "black", margin = margin(t = 0, r = 30, b = 0, l = 0)),
    axis.text.x = element_text(size = 35, angle = 90, color = "black"),   
    legend.text = element_text(size = 43), # Removed face = "italic"
    legend.position = "bottom",
    legend.spacing.x = unit(0.1, 'mm'),
    legend.spacing.y = unit(0.05, 'mm'),
    legend.location = "plot",
    legend.title = element_text(face = "bold", size = 45)
  ) +
  guides(fill = guide_legend(ncol = 6, byrow = TRUE)) + 
  labs(fill = 'Bacterial Groups')

ggsave(
  filename = "/local1/workdir1/tw488/CEA_RData_RScript/manuscript_data_analysis/samples_Fam_relabbarplot772026.tiff", 
  plot = last_plot(),       
  device = "tiff",          
  units = "in", 
  width = 45, 
  height = 25, 
  dpi = 700, 
  compression = "lzw",
  type = "cairo")


#list of genera for organic and conventional samples in fam_rel_abund_ps.noncontam4- do not want under 5% grouped together

fam_org<-fam_rel_abund_ps.noncontam4%>%filter(Pond_Type=="Organic")%>%
  group_by(taxon)%>%
  mutate(Avg = mean(rel_abund))%>%distinct(taxon,Avg)%>%
  arrange(desc(Avg))%>%filter(Avg!="0")
colnames(fam_org)[1] <- "Organic Sample Families"
colnames(fam_org)[2] <- "Average Relative Abundance (%)"
write.csv(fam_org,"/local1/workdir1/tw488/CEA_RData_RScript/manuscript_data_analysis/organic_fam_relab", row.names = TRUE)


fam_convent<-fam_rel_abund_ps.noncontam4%>%
  filter(Pond_Type=="Conventional")%>%group_by(taxon)%>%
  mutate(Avg = mean(rel_abund))%>%distinct(taxon,Avg)%>%
  arrange(desc(Avg))%>%filter(Avg!="0")
colnames(fam_convent)[1] <- "Conventional Sample Families"
colnames(fam_convent)[2] <- "Average Relative Abundance (%)"
write.csv(fam_convent,"/local1/workdir1/tw488/CEA_RData_RScript/manuscript_data_analysis/Conventional_fam_relab", row.names = TRUE)


#ALPHA DIVERSITY INDICES AND TFSOS
alphadiv<-estimate_richness(ps.noncontam4, split = TRUE, measures = NULL)
alphadiv$sample <- row.names(alphadiv)  
alphadiv <- alphadiv %>% separate(sample, c("Pond_Type", "Pond_Number", "Sample_Location"), "_")
alphadiv$sample <- row.names(alphadiv) 
alphadiv$Pond_Type_Number = paste(alphadiv$Pond_Type, alphadiv$Pond_Number, sep="_")

#calculate Pielou's Eveness
alphadiv <- alphadiv %>% mutate (Pielou=(alphadiv$Shannon)/log(alphadiv$Observed))

#ALPHA DIV
r<-plot_richness(ps.noncontam4, color="Pond_Type", measures=c("Observed","Eveness"))+ 
  geom_boxplot()+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))+
  theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))
r
#ggsave("alphadiv.tiff", units="in", width=12, height=6, dpi=300, compression = 'lzw')


#Pond Type: Conventional vs Organic
#-	No treatment or Permute pond_type_number together
#Sample Location: Harvest, Middle, Seeding
#Need to account for organic/ conventional and which pond it came form, or just organic/ conventional


#PERMANOVA for pond type

ps.noncontotu<-pstoveg_otu(ps.noncontam4)
ps.noncontsample<-pstoveg_sample(ps.noncontam4)
ps.noncontsd<-pstoveg_sd(ps.noncontam4)

#make in to distance object
ps.noncontotu_vegdist<-vegdist(ps.noncontotu, method="bray")

#combine otu and sd
ps.noncontsd_otu<-cbind(ps.noncontsd,ps.noncontotu)

#how function from permute package
how_pondtype <- with(ps.noncontsd_otu, how(plots=Plots(strata=Pond_Type_Number, type="free")))
how_pondtype #within=Within(type="free"),

#ensure replicates are permuted together
testsd<-ps.noncontsd_otu
testsd <- as_tibble(testsd, rownames = "SampleID")
#the number below should be the number of samples you have (the number of rows in testsd)
testsd$permute1<-shuffle(12,control=how_pondtype)
testsd$origorder<-1:12
testsd<-as.data.frame(testsd)

#visualize permutations, change 'group' and 'color' to whatever variable that you want to group samples by
testsd%>%pivot_longer(permute1:origorder, names_to="column", values_to="position", cols_vary = "slowest") %>%
  ggplot(aes(x=column, y=position, group=SampleID, color=Pond_Type_Number))+geom_line()


#check to see even # of samples
table(ps.noncontsd$Pond_Type_Number)
table(ps.noncontsd$Sample_Location)


#plug in here
ps.noncontotu_perm_pondtype <- adonis2(ps.noncontotu_vegdist~Pond_Type*Sample_Location, data=ps.noncontsd, permutations = how_pondtype, method = "bray", 
                                       sqrt.dist = FALSE, add = FALSE, by = "terms",
                                       parallel = getOption("mc.cores"), na.action = na.fail)
ps.noncontotu_perm_pondtype


#NC sample analysis, for unrarefied NC samples:ps.noncontamNC 
ps.NCsampledata<-as.data.frame(ps.NC@sam_data)
ps.NCtbl<-as.data.frame(ps.NC@otu_table)
ps.NCtxtbl<-as.data.frame(ps.NC@tax_table)
ps.NCotu_colsums<-as.data.frame(colSums(ps.NCtbl))
ps.NCotu_rowsums<-as.data.frame(rowSums(ps.NCtbl))

ps.NC.counts <- count_seqs(ps.NC@otu_table)

abund.ps.noncontamNC <- ps.NC@otu_table %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("ASV_") %>%
  tidyr::pivot_longer((-ASV_), names_to = "sample") %>%
  merge(ps.NC.counts$counts, by="sample") %>%
  dplyr::select(-n) %>%
  mutate(newcol=sample)%>%
  tidyr::separate(newcol, c("Pond_Type","Pond_Number","Location"), "_",remove = TRUE)%>%
  group_by(sample) %>%
  mutate(rel_abund = value/counts) %>% #relative abundance by treatment/day/rep
  ungroup() 

#check that relative abundance is calculated correctly- everything should add up to one
abund.ps.noncontamNC %>% group_by(sample) %>% summarise(sum = sum(rel_abund))

#get taxonomy
ps.noncontamNCtax <- ps.NC@tax_table %>%
  as.data.frame() %>%
  rownames_to_column("ASV_") 
head(ps.noncontamNCtax )

#merge taxonomy with rarefied relative abundance
count_tax.ps.noncontamNC <- merge(abund.ps.noncontamNC,ps.noncontamNCtax, by="ASV_")
head(count_tax.ps.noncontamNC)

# Get summary abundance of each taxa
taxa.meta.sumps.noncontamNC <- count_tax.ps.noncontamNC%>%
  dplyr::select(-c(value, counts)) %>%
  tidyr::pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "ASV_"),
                      names_to = "level", 
                      values_to = "taxon")

#format phylum names for GENUS
genus_rel_abund_ps.noncontamNC <- taxa.meta.sumps.noncontamNC %>%
  filter(level == "Genus") %>%
  group_by(taxon,sample,Pond_Type) %>% 
  summarise(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  mutate(taxon = stringr::str_replace(taxon,
                                      "^unclassified (.*)", "Unclassified *\\1*"),
         taxon = stringr::str_replace(taxon, "^(\\S*)$", "\\1")) %>%
  ungroup()
genus_rel_abund_ps.noncontamNC$taxon <-genus_rel_abund_ps.noncontamNC$taxon%>%tidyr::replace_na('Unclassified')


#Filter out genera that have 0% relative abundance
genus_rel_abund_ps.noncontamNC_no0relab <- genus_rel_abund_ps.noncontamNC%>% filter(rel_abund!="0") 

#rel ab of genera for org and convent samples in join_taxon_rel_abund_ps.noncontam4

#genera_org<-genus_rel_abund_ps.noncontamNC_no0relab

genus_rel_abund_ps.noncontamNC_no0relab_4export<-dcast(genus_rel_abund_ps.noncontamNC_no0relab,sample~taxon,value.var="rel_abund")
genus_rel_abund_ps.noncontamNC_no0relab_4export<-genus_rel_abund_ps.noncontamNC_no0relab_4export%>%replace(is.na(.), 0)
write.csv(genus_rel_abund_ps.noncontamNC_no0relab_4export,"/local1/workdir1/tw488/CEA_RData_RScript/manuscript_data_analysis/NC_samples_relab", row.names = TRUE)

#Family analysis for NC samples
#NC sample analysis, for unrarefied NC samples:ps.noncontamNC 

#format phylum names for family
fam_rel_abund_ps.noncontamNC <- taxa.meta.sumps.noncontamNC %>%
  filter(level == "Family") %>%
  group_by(taxon,sample,Pond_Type) %>% 
  summarise(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  mutate(taxon = stringr::str_replace(taxon,
                                      "^unclassified (.*)", "Unclassified *\\1*"),
         taxon = stringr::str_replace(taxon, "^(\\S*)$", "\\1")) %>%
  ungroup()
fam_rel_abund_ps.noncontamNC$taxon <-fam_rel_abund_ps.noncontamNC$taxon%>%tidyr::replace_na('Unclassified')


#pool taxon <5% relative abundance into "other" category
fam_rel_abund_ps.noncontamNC_pool <- fam_rel_abund_ps.noncontamNC %>%
  group_by(sample, taxon, Pond_Type) %>%
  summarise(mean = mean(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarise(pool = max(mean) <5,
            mean=mean(mean), 
            .groups="drop")

join_fam_rel_abund_ps.noncontamNC <- inner_join(fam_rel_abund_ps.noncontamNC,fam_rel_abund_ps.noncontamNC_pool, by="taxon") %>%
  mutate(taxon = ifelse(pool, "<5% Relative Abundance", taxon)) %>%
  group_by(sample, taxon,Pond_Type) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  mutate(taxon = factor(taxon),
         taxon = forcats::fct_reorder(taxon, rel_abund, .desc=TRUE),
         taxon = forcats::fct_shift(taxon, n=1)) %>%
  ungroup() %>%
  unique() %>%
  group_by(taxon) %>%
  arrange(desc(rel_abund), .by_group=TRUE) 




#rel ab genera for org and convent samples in join_taxon_rel_abund_ps.noncontam4
#filter out families that have 0% rel ab
fam_rel_abund_ps.noncontamNC_no0relab <- fam_rel_abund_ps.noncontamNC%>% filter(rel_abund!="0") 

fam_rel_abund_ps.noncontamNC_no0relab_4export<-dcast(fam_rel_abund_ps.noncontamNC_no0relab,sample~taxon,value.var="rel_abund")
fam_rel_abund_ps.noncontamNC_no0relab_4export<-fam_rel_abund_ps.noncontamNC_no0relab_4export%>%replace(is.na(.), 0)
write.csv(fam_rel_abund_ps.noncontamNC_no0relab_4export,"/local1/workdir1/tw488/CEA_RData_RScript/manuscript_data_analysis/NC_samples_relab_fam", row.names = TRUE)
