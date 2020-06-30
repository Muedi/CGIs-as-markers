
##### RnBeads for differentially methylated HCCs and blood ####
## This script exctracts information from .idat raw data from Illumina and from GEO_series data and computes CGI-methylation.
## it returns said methylation for the sites and islands as .csv-files.

library(RnBeads)
library(RnBeads.hg19)
library(GEOquery)
library(stringr)
library(minfi)


#####Ghostscript and ZIP environment, see RnBeads documentation/FAQ
Sys.setenv("R_GSCMD" = "C:/Program Files/gs/gs9.50/bin/gswin64c.exe")
Sys.setenv("R_ZIPCMD" = "C:/RTools/bin/zip.exe")
Sys.setenv("R_UNZIPCMD" = "C:/RTools/bin/unzip.exe")

##### parallel
CORES <- 6

parallel.setup(CORES)

#####set your R options####
rnb.options(identifiers.column="sample_Name",
            # change this to F if your computer has much memory
            disk.dump.big.matrices = T,
            normalization.method="none",
            differential.comparison.columns = "sample_Group",
            differential.enrichment.lola = F,
            covariate.adjustment.columns = c("sample_Name", "sample_Group"),
            region.types = "cpgislands"
)

#####setting up folders #####
data.dir <- "C:/Users/User/project_folder"

idat.dir <- file.path("C:/Users/User/project_folder/idat")
analysis.dir <- file.path(data.dir, "Analysis")
save.dir <- file.path(data.dir, "saved_sets")

# initialize Sampla annotations and data sources
# the sample sheets need to be produced by you for the respective data you want to read in.
# RnBeads provides an example on how to do that in the Documantation
# PBMC
sample.annotation.PBMC  <- file.path(data.dir, "Sample_sheet_PBMC.csv")
# whole blood and cfDNA control
sample.annotation.whole_blood  <- file.path(data.dir, "Sample_sheet_wb.csv")
sample.annotation.cfDNA_ctrl <- file.path(data.dir, "Sample_sheet_cf_moss.csv")
# solid tumor
sample.annotation.HCC_2 <- file.path(data.dir, "Sample_sheet_HCC_2.csv")

# collects idats from idat-folder only needed for data which is present in idat format
data.source.PBMC <- c(idat.dir, sample.annotation.PBMC)
data.source.whole_blood <- c(idat.dir, sample.annotation.whole_blood)
data.source.cfDNA_crtl <- c(idat.dir, sample.annotation.cfDNA_ctrl)
data.source.HCC_2 <- c(idat.dir, sample.annotation.HCC_2)

# report directories, one for every datasource (RnBeads needs to initialize a folder to drop its QC and preprocessing data.)
report.dir <- file.path(analysis.dir, "reports_details")
report.dir.2 <- file.path(analysis.dir, "reports_details_2")
report.dir.3 <- file.path(analysis.dir, "reports_details_3")
report.dir.4 <- file.path(analysis.dir, "reports_details_4")
report.dir.5 <- file.path(analysis.dir, "reports_details_5")
report.dir.6 <- file.path(analysis.dir, "reports_details_6")
report.dir.7 <- file.path(analysis.dir, "reports_details_7")
report.dir.8 <- file.path(analysis.dir, "reports_details_8")

###### load data from idat: whole blood (healthy control) data from GEO ######
rnb.initialize.reports(report.dir)
logger.start(fname=NA)

result <- rnb.run.import(data.source=data.source.whole_blood,
                         data.type="infinium.idat.dir",
                         dir.reports=report.dir)

rnb.set.whole_blood <- result$rnb.set
remove(result)
# qc + preprocessing
rnb.run.qc(rnb.set.whole_blood, report.dir)
tmp <- rnb.run.preprocessing(rnb.set.whole_blood, dir.reports=report.dir)
rnb.set.whole_blood <- tmp$rnb.set
remove(tmp)

####load data from idat: peripheral blood (leucocyte DNA, older patients) data from GEO ######
rnb.initialize.reports(report.dir.2)
logger.start(fname=NA)

result2 <- rnb.run.import(data.source=data.source.PBMC,
                          data.type="infinium.idat.dir",
                          dir.reports=report.dir.2)

rnb.set.PBMC <- result2$rnb.set
# sample_group of pbmc is tagged as blood in the samplesheets -> PBMC
rnb.set.PBMC@pheno$sample_Group <- "PBMC"

# qc + preprocessing
rnb.run.qc(rnb.set.PBMC, report.dir.2)
tmp <- rnb.run.preprocessing(rnb.set.PBMC, dir.reports=report.dir.2)
rnb.set.PBMC <- tmp$rnb.set
remove(tmp)


###### load data from idat: cfDNA Moss 4 healthy pools, 4 healthy solo, data from GEO ######
rnb.initialize.reports(report.dir.3)
logger.start(fname=NA)

result3 <- rnb.run.import(data.source = data.source.cfDNA_crtl,
                          data.type="infinium.idat.dir",
                          dir.reports = report.dir.3)

rnb.set.cf_moss <- result3$rnb.set
remove(result3)
# qc + preprocessing
rnb.run.qc(rnb.set.cf_moss, report.dir.3)
tmp <- rnb.run.preprocessing(rnb.set.cf_moss, dir.reports=report.dir.3)
rnb.set.cf_moss <- tmp$rnb.set
remove(tmp)


#### Load data from GEO matrix series: cfDNA pooled GSE110185
rnb.set.cf_pooled <- rnb.read.geo("GSE110185_series_matrix.txt.gz")
samples_vs_title <- data.frame(samples(rnb.set.cf_pooled), rnb.set.cf_pooled@pheno$title)
to_remove <- which(!(samples_vs_title$rnb.set.cf_pooled.pheno.title %in% str_subset(samples_vs_title$rnb.set.cf_pooled.pheno.title, "NCF")))
rnb.set.cf_pooled <- remove.samples(rnb.set.cf_pooled, to_remove)


rnb.set.cf_pooled@pheno$sample_ID <- rnb.set.cf_pooled@pheno$geo_accession
rnb.set.cf_pooled@pheno$sample_Group <- "Blood"
rnb.set.cf_pooled@pheno$sample_Name <- rnb.set.cf_pooled@pheno$title
rnb.set.cf_pooled@pheno$sub_Group <- "cf_DNA_pool"

# init the reportfolder
rnb.initialize.reports(report.dir.4)
logger.start(fname=NA)

# qc + preprocessing
rnb.run.qc(rnb.set.cf_pooled, report.dir.4)
tmp <- rnb.run.preprocessing(rnb.set.cf_pooled, dir.reports=report.dir.4)
rnb.set.cf_pooled <- tmp$rnb.set
remove(tmp)

#### Load data from GEO matrix series: cfDNA HCC GSE129374
rnb.set.cf_HCC <- rnb.read.geo("GSE129374_series_matrix.txt.gz")
samples_vs_title <- data.frame(samples(rnb.set.cf_HCC), rnb.set.cf_HCC@pheno$title)
to_remove <- which(!(samples_vs_title$rnb.set.cf_HCC.pheno.title %in% str_subset(samples_vs_title$rnb.set.cf_HCC.pheno.title, "HCC")))
rnb.set.cf_HCC <- remove.samples(rnb.set.cf_HCC, to_remove)

rnb.set.cf_HCC@pheno$sample_ID <- rnb.set.cf_HCC@pheno$geo_accession
rnb.set.cf_HCC@pheno$sample_Group <- "cf_HCC"
rnb.set.cf_HCC@pheno$sample_Name <- rnb.set.cf_HCC@pheno$title

# init the reportfolder
rnb.initialize.reports(report.dir.5)
logger.start(fname=NA)

# qc + preprocessing
rnb.run.qc(rnb.set.cf_HCC, report.dir.5)
tmp <- rnb.run.preprocessing(rnb.set.cf_HCC, dir.reports=report.dir.5)
rnb.set.cf_HCC <- tmp$rnb.set
remove(tmp)
save.rnb.set(rnb.set.cf_HCC, file.path(save.dir, "rnb.set.cf_HCC"), archive = TRUE)

#### Load data from GEO matrix series: solid tumor HCC GSE99036
rnb.set.HCC_solid <- rnb.read.geo("GSE99036_series_matrix.txt.gz")
samples_vs_title <- data.frame(samples(rnb.set.HCC_solid), rnb.set.HCC_solid@pheno$title)
to_remove <- which(!(samples_vs_title$rnb.set.HCC_solid.pheno.title %in% str_subset(samples_vs_title$rnb.set.HCC_solid.pheno.title, "carcinoma")))
rnb.set.HCC_solid <- remove.samples(rnb.set.HCC_solid, to_remove)

rnb.set.HCC_solid@pheno$sample_ID <- rnb.set.HCC_solid@pheno$geo_accession
rnb.set.HCC_solid@pheno$sample_Group <- "HCC_solid"
rnb.set.HCC_solid@pheno$sample_Name <- rnb.set.HCC_solid@pheno$title

# init the reportfolder
rnb.initialize.reports(report.dir.6)
logger.start(fname=NA)

# qc + preprocessing
rnb.run.qc(rnb.set.HCC_solid, report.dir.6)
tmp <- rnb.run.preprocessing(rnb.set.HCC_solid, dir.reports=report.dir.6)
rnb.set.HCC_solid <- tmp$rnb.set
remove(tmp)
save.rnb.set(rnb.set.HCC_solid, file.path(save.dir, "rnb.set.HCC_solid"), archive = TRUE)

## Load data directly from idat: solid tumor HCC GSE77269
# the data was not online with the complete barcode and sentrix ID which is needed by RnBeads
# therefore i used minfi to load the data and convert the rgset to a rnb.set

hcc_files <- list.files(idat.dir, pattern = "GSM20470")
basenames <- grep("_Grn.idat", hcc_files, value = T)
basenames <- str_replace_all(basenames, "_Grn.idat", "")
minfi.set.HCC <- read.metharray(file.path(idat.dir, basenames), force = T)
rnb.set.HCC_solid_2 <- as(minfi.set.HCC, "RnBeadRawSet")
remove(minfi.set.HCC)
sample.IDs <- list()
for (i in 1:20) {
  sample.IDs[i] <- paste("HCC_2_", i, sep = "", collapse = NULL)
}

rnb.set.HCC_solid_2@pheno$sample_ID <- sample.IDs
rnb.set.HCC_solid_2@pheno$sample_Group <- "HCC_solid"
rnb.set.HCC_solid_2@pheno$sample_Name <- sample.IDs

# init the reportfolder
rnb.initialize.reports(report.dir.7)
logger.start(fname=NA)
# qc + preprocessing
rnb.run.qc(rnb.set.HCC_solid_2, report.dir.7)
tmp <- rnb.run.preprocessing(rnb.set.HCC_solid_2, dir.reports=report.dir.7)
rnb.set.HCC_solid_2 <- tmp$rnb.set
remove(tmp)
save.rnb.set(rnb.set.HCC_solid_2, file.path(save.dir, "rnb.set.HCC_solid_2"), archive = TRUE)



#### Load data directly from GEO matrix series: whole blood 2 GSE40279
wb_GSE40279 <- rnb.read.geo("GSE40279_series_matrix.txt.gz")
# extract the patients underyounger than 49
samples_vs_age <- data.frame(samples(wb_GSE40279), wb_GSE40279@pheno$`age (y)`)
to_remove <- which(!(samples_vs_age$wb_GSE40279.pheno..age..y.. < 49))
wb_GSE40279 <- remove.samples(wb_GSE40279, to_remove)

wb_GSE40279@pheno$sample_ID <- wb_GSE40279@pheno$geo_accession
wb_GSE40279@pheno$sample_Group <- "Blood"
wb_GSE40279@pheno$sample_Name <- wb_GSE40279@pheno$geo_accession

rnb.initialize.reports(report.dir.8)
logger.start(fname=NA)

# qc + preprocessing
rnb.run.qc(wb_GSE40279, report.dir.8)

tmp <- rnb.run.preprocessing(wb_GSE40279, dir.reports=report.dir.8)
wb_GSE40279 <- tmp$rnb.set
remove(tmp)
# save set and remove to keep light profile
save.rnb.set(wb_GSE40279, file.path(save.dir, "wb_GSE40279"), archive=TRUE)
# remove(wb_GSE40279)

####### after the first run all sets are saved and can be loaded again to save time, if any problems occured or data will be added. ########
# load sets
# rnb.set.PBMC <- load.rnb.set(file.path(save.dir, "rnb.set.PBMC.zip"))
# rnb.set.whole_blood <- load.rnb.set(file.path(save.dir, "rnb.set.whole_blood.zip"))
# rnb.set.cf_moss <- load.rnb.set(file.path(save.dir, "rnb.set.cf_moss.zip"))
# rnb.set.cf_pooled <- load.rnb.set(file.path(save.dir, "cf_GSE110185.zip"))
# wb_GSE40279 <- load.rnb.set(file.path(save.dir, "wb_GSE40279.zip"))

# rnb.set.cf_HCC <- load.rnb.set(file.path(save.dir, "rnb.set.cf_HCC.zip"))
# rnb.set.HCC_solid <- load.rnb.set(file.path(save.dir, "rnb.set.HCC_solid.zip"))
# rnb.set.HCC_solid_2 <- load.rnb.set(file.path(save.dir, "rnb.set.HCC_solid_2.zip"))


# list of RnBeadRawSet
rnb.sets <-list()
rnb.sets[['PBMC']]<- rnb.set.PBMC
rnb.sets[['whole_blood']]<- rnb.set.whole_blood
rnb.sets[['whole_blood_2']]<- wb_GSE40279
rnb.sets[['cf_moss']]<- rnb.set.cf_moss
rnb.sets[['cf_pooled']]<- rnb.set.cf_pooled

rnb.sets[['HCC_solid']]<- rnb.set.HCC_solid
rnb.sets[['HCC_solid_2']]<- rnb.set.HCC_solid_2
rnb.sets[['cf_HCC']]<- rnb.set.cf_HCC

# complete samples as main input for python
rnb.set.complete <- rnb.sets[['PBMC']]

for (i in 2:8) {
  rnb.set.complete <- rnb.combine.arrays(rnb.set.complete, rnb.sets[[i]], type="common")
}

# adding information about the data to the new rnb.set
rnb.set.complete<-addPheno(
  rnb.set.complete, unlist(sapply(names(rnb.sets[1:9]), function(nn) rep(nn, length(samples(rnb.sets[[nn]]))))),
  "Platform")

rnb.set.complete<-addPheno(rnb.set.complete, paste(pheno(rnb.set.complete)[["sample_Name"]], pheno(rnb.set.complete)[["type"]], sep=" "), "SampleID")


# execute normalization over the respective combinations to increase comparability
rnb.set.complete <- rnb.execute.normalization(rnb.set.complete, method = "bmiq")

# saving
# save.dir <- file.path(data.dir, "saved_sets")
# save.rnb.set(rnb.set.complete, file.path(save.dir, "rnb.set.complete"), archive = TRUE)

# same as above
rnb.set.complete <- load.rnb.set(file.path(save.dir, "rnb.set.complete.zip"))
# rnb.set.blood_vs_tumor <- load.rnb.set(file.path(save.dir, "rnb.set.blood_vs_tumor.zip"))
# rnb.set.tumor_vs_PBMC <- load.rnb.set(file.path(save.dir, "rnb.set.tumor_vs_PBMC.zip"))
# rnb.set.blood_vs_PBMC <- load.rnb.set(file.path(save.dir, "rnb.set.blood_vs_PBMC.zip"))

# writes out the cpg sites with correspondiong index of the beta values in the second table.

# methylation level
meth.vals <- meth(rnb.set.complete, row.names=TRUE)
dimnames(meth.vals)[[2]] <- rnb.set.complete@pheno$sample_Name 
write.table(meth.vals, file = "C:/Users/User/project_folder/Output_RnBeads/complete_HCC_sites_meth_new2.csv", quote = FALSE, sep = ",")

# Annotation of islands, saving as mean meth of islands
# RnBeads does not return the meth values with CGIs as index,
# but the indices of the respective dataframes are corresponding to the same CGI, so they can be fused together in python without problem
# index = annotation
# saving anno
island_anno <- annotation(rnb.set.complete, type="cpgislands")
write.table(island_anno, file = "C:/Users/User/project_folder/Output_RnBeads/complete_HCC_islands_anno_new2.csv", quote = FALSE, sep = ",")
# saving meth-values
write.table(as.ffdf(rnb.set.complete@meth.regions$cpgislands), file ="C:/Users/User/project_folder/Output_RnBeads/complete_HCC_islands_meth_new2.csv", quote = FALSE, sep = ",")
