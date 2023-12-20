# Downloaded breast cancer Dataset ('Breast Invasive Carcinoma (TCGA, PanCancer Atlas)') from cbioportal
#1.	Download the dataset on: https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018
# Set own directory path 
path = "S:\Sathya\UCD\Autumn Sememster - 2023\ANAT40040-Bio Principles & Cellular Org-2\Assignment\Assignment 2"

# set file name 
file_location = "brca_tcga_pan_can_atlas_2018.tar.gz"

# 2.	Untar the folder and extract the files / unzip / extract ZIP file and pointing to the extraxted location
untar(file_location)

new_file_location = paste(getwd(), "brca_tcga_pan_can_atlas_2018", sep = "/")

setwd(new_file_location)

patient_data  <- read.delim("data_clinical_patient.txt")

Patient_age <- which(colnames(patient_data) =="Diagnosis.Age")

# we will use following files:
# data_clinical_patient.txt, data_mrna_seq_v2_rsem.txt, data_mutations.txt and data_cna.txt


#3.	Read the RNASeq file: data_mrna_seq_v2_rsem.txt
rnaseq = read.delim("data_mrna_seq_v2_rsem.txt")

# we will delete the genes for which there's more than one Hugo Symbol
# These are typically genes with no Hugo Symbol ("" as an entry) or pseudogenes.

keep = !duplicated(rnaseq[,1])

rnaseq = rnaseq[keep,]

# set rownames of rnaseq to hugo symbols

rownames(rnaseq)  = rnaseq[,1]

#4.	Read the Patient Data file: data_clinical_patient.txt
clinical = read.delim("data_clinical_patient.txt")

#5.	Read the Copy Number Aberrations Data: data_cna.txt
cna = read.delim("data_cna.txt")

# 7.	Create metadata using the CNA level of ERBB2+ (greater than 0 means amplified).
erbb2_indx = which(cna[,1] == 'ERBB2')

# Plot histogram to visualize explore the data.

hist(as.numeric(cna[erbb2_indx,-c(1,2)]))

#6.	Match the RNASeq patient ids with the CNA ids and the Patient Data ids.
rna_cna_id = which(is.element(colnames(rnaseq[,-c(1,2)]), colnames(cna[,-c(1,2)])))

# select only the rna cases which have cna data.
rna_cna_sub = rnaseq[,2+rna_cna_id]

# check all patients in rna_can_sub are in cna
no_pats_in_rna_cna_sub_and_cna = sum(is.element(colnames(rnaseq[,2+rna_cna_id]), colnames(cna[,-c(1,2)]))) 

# sanity check.This will print an error if the result is not the same.

sanity_check = no_pats_in_rna_cna_sub_and_cna == dim(rna_cna_sub)[2]

# Pre-allocate memory for ERBB2
meta_erbb2 = matrix(0,length(rna_cna_id),1)

for (i in 1:length(rna_cna_id)){
  # access the colnames of i
  col_i = colnames(rna_cna_sub)[i]
  # get the index in cna for the same patient
  col_cna = which(colnames(cna)==col_i)
  # store if they're amplified.
  meta_erbb2[i,] = 1*(cna[erbb2_indx,col_cna]>0)
  
}

# This are some checks you can do to make sure your code worked.
# There's some more systematic checks you can do. See unit testing.

# simple checks to make sure. 

col_i = colnames(rna_cna_sub)[1]
col_cna = which(colnames(cna)==col_i)

# sanity check
(cna[erbb2_indx,col_cna]>0) == meta_erbb2[1,1]

# see now if a positive meta_erbb2 is amplified.
pos_example = which(meta_erbb2==1)[1]
col_i = colnames(rna_cna_sub)[pos_example]
col_cna = which(colnames(cna)==col_i)

# sanity check
(cna[erbb2_indx,col_cna]>0) == meta_erbb2[pos_example,1]

# batch checks should print true.

# We will add a title to the metadata.
colnames(meta_erbb2) = 'ERBB2Amp'

# transform into integers
rna_cna_sub = round(rna_cna_sub)

#8. Normalize data using DESeq2.
#The purpose of normalization is to eliminate systematic effects that are not associated with the biological differences of interest.
# Install DESeq2.
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install DeSeq2

BiocManager::install("DESeq2")

library(DESeq2)
#rna_cna_id
#Refer - https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
#Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = meta, 
                              design = ~ sampletype)

#we nested it within the View() function so that rather than getting printed in the console we can see it in the script editor
View(counts(dds))
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

#9. Obtain Differentially Expressed Genes.

dds <- DESeq(dds)
res <- result(dds)

#Change DESeq Object to R Object (data.frame)
res <- as.data.frame(res)
class(res)

head(res)

sigs <- na.omit(res)
sigs <- sigs[sigs$padj < 0.05 & sigs$baseMean > 50,]


# explore result
summary(res)

#10. Perform a Pathway Enrichment Analysis

# Install Cluster profiler
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Clusterprofiler
BiocManager::install("ClusterProfiler")

# Install AnnotationProfiler
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install AnnotationProfiler
BiocManager::install("AnnotationProfiler")

# Install Homosapiens Database
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install org.Hs.eg.db
BiocManager::install("org.Hs.eg.db") 

library(clusterProfiler)
library(AnnotationProfiler)
library(org.Hs.eg.db)

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "Id")
as.data.frame(GO_results)
fit <- plot(barplot(GO_results, showCategory = 15))

png("out.png", res = 250, width = 1400, height = 1800)
print(fit)
dev.off()

fit


#11. Get the variance stabilised transformed expression values.
vsd <- vst(dds.blind=FALSE)

#12. With the vst values obtain a PCA plot.
#Refer: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#data-quality-assessment-by-sample-clustering-and-visualization
# search with Heatmap of the count matrix
library("pheatmap")
plotPCA(vsd,intgroup=c("sequencing","Treatment"))

#Top 10 genes
top_10 <- res[order(res$padj),][1:10,]
top_10 <- row.names(top_10)
top_10

rld <- rlog(dds,blind=FALSE)

pheatmap(assay(rld)[top_10,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)

#with Clustering
pheatmap(assay(rld)[top_10,])

annot_info <- as.data.frame(colData(dds),c('S', 'T'))
pheatmap(assay(rld)[top_10,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=annot_info)



#Optional for Additional Marks

#13 Cluster the data and show in PCA



#14 With the vst values of the DE genes generate an overall survival model.




#15 Use lasso cross validation to find the set of genes which predict survival the best.

df.t <- structure(list(hsa_miR_105_5p = c(3.58497328179801, 5.73145238130165, 
1.19037294682376, -1.28586123284671, 1.27004401721869, 0.958088884635556
), hsa_miR_17_3p = c(1.21345556145455, 4.71642723353062, 5.87616915208789, 
0.776249937585565, 4.86437477300888, 1.71876771352689), hsa_miR_3916 = c(6.74863569372315, 
3.23155618956527, -0.105259761381448, -1.28586123284671, 4.60953338597123, 
2.95060221832751), hsa_miR_1295a = c(-1.35668910756094, 0.147551018264645, 
2.44220202218853, -1.28586123284671, 5.47367734142336, -0.135507425889107
)), row.names = c("86", "175", "217", "394", "444", "618"), class = "data.frame")

Time <- structure(c(1796, 1644.04166666667, 606.041666666667, 1327.04166666667, 
665, 2461), class = "difftime", units = "days")

Status <- c(0L, 0L, 1L, 0L, 1L, 0L)


cox.out <- capture.output(for(i  in colnames(df.t)){
  print(summary(coxph(as.formula(paste0("Surv(Time, Status)~", i )),  data=as.data.frame(df.t))))
})









