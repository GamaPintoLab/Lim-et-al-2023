# -- LARC data -- #
# IMPORTANT FEATURES OF THE DATA:
# This script is used to do the analysis of the CTC dataset from the LARC study
# 3 data points for baseline (1base), during treatment (2dur) and after the treatment (3aft)
# 128 features including 16 BLANKS that were removed from original EDS files before the analysis
# 112 features corresponding to 106 miRNAs and 6 endogenous controls (3 endo controls in duplicates)
# 3 batches with the following patterns: first batch - TEV, second batch - VCE, third batch - YTK03, YTK04

#IMPORTANT FEATURES OF THIS SCRIPT:
# Packages to be used: HTqPCR, DataCombine, genefilter, limma
#.........................................................................................................................
#----------------------------------------------------------------- FIRST PART --------------------------------------------------------------------------------------------------------------------------------------

                                # BEGIN HERE FOR THE HTqPCR EXPLORATORY AND DIFFERENTIAL EXPRESSION ANALYSIS 
                                                                      # |
                                                                      # |
                                                                      # V

#.................................................................................................................................................................................................................
# -- Loading of data -- # BEGIN -----------------------------------------------------------------------------------


setwd("/LARC_project/")

TEV_64=readxl::read_xlsx("TEV64_CTC_21_QuantStudio_export.xlsx")
temp=TEV_64[19:3091,]
colnames(temp)=temp[1,]
tev_64=temp[2:3073,]
rm(temp, TEV_64)
head(tev_64)

TEV_65=readxl::read_xlsx("TEV65_CTC_21_QuantStudio_export.xlsx")
temp=TEV_65[19:3091,]
colnames(temp)=temp[1,]
tev_65=temp[2:3073,]
rm(temp, TEV_65)
head(tev_65)

TEV_71=readxl::read_xlsx("TEV71_CTC_21_QuantStudio_export.xlsx")
temp=TEV_71[19:3091,]
colnames(temp)=temp[1,]
tev_71=temp[2:3073,]
rm(temp, TEV_71)
head(tev_71)

VCE_88=readxl::read_xlsx("VCE88_labelled_CTC_19_QuantStudio_export.xlsx")
temp=VCE_88[21:3093,]
colnames(temp)=temp[1,]
vce_88=temp[2:3073,]
rm(temp, VCE_88)
head(vce_88)

VCE_93=readxl::read_xlsx("VCE93_labelled_CTC_19_QuantStudio_export.xlsx")
temp=VCE_93[21:3093,]
colnames(temp)=temp[1,]
vce_93=temp[2:3073,]
rm(temp, VCE_93)
head(vce_93)

VCE_95=readxl::read_xlsx("VCE95_labelled_CTC_19_QuantStudio_export.xlsx")
temp=VCE_95[21:3093,]
colnames(temp)=temp[1,]
vce_95=temp[2:3073,]
rm(temp, VCE_95)
head(vce_95)

YTK02=readxl::read_xlsx("YTK02_QuantStudio_export.xlsx")
temp=YTK02[21:3093,]
colnames(temp)=temp[1,]
ytk02=temp[2:3073,]
rm(temp, YTK02)
head(ytk02)

YTK03=readxl::read_xlsx("YTK03_QuantStudio_export.xlsx")
temp=YTK03[21:3093,]
colnames(temp)=temp[1,]
ytk03=temp[2:3073,]
rm(temp, YTK03)
head(ytk03)

YTK04=readxl::read_xlsx("YTK04_QuantStudio_export.xlsx")
temp=YTK04[21:3093,]
colnames(temp)=temp[1,]
ytk04=temp[2:3073,]
rm(temp, YTK04)
head(ytk04)

YTK05=readxl::read_xlsx("YTK05_QuantStudio_export.xlsx")
temp=YTK05[21:3093,]
colnames(temp)=temp[1,]
ytk05=temp[2:3073,]
rm(temp, YTK05)
head(ytk05)

YTK07=readxl::read_xlsx("YTK07_QuantStudio_export.xlsx")
temp=YTK07[21:3093,]
colnames(temp)=temp[1,]
ytk07=temp[2:3073,]
rm(temp, YTK07)
head(ytk07)


tev_64=as.data.frame(tev_64); tev_65=as.data.frame(tev_65); tev_71=as.data.frame(tev_71)

vce_88=as.data.frame(vce_88); vce_93=as.data.frame(vce_93); vce_95=as.data.frame(vce_95)

ytk02=as.data.frame(ytk02); ytk03=as.data.frame(ytk03); ytk04=as.data.frame(ytk04); ytk05=as.data.frame(ytk05);ytk07=as.data.frame(ytk07)

# A tibble: 6 x 16
# Well  `Well Position` Omit  `Sample Name` `Target Name` Task  Reporter Quencher Crt   `Crt Mean` `Crt SD` `Amp Score` `Cq Conf` CRTAMPLITUDE CRTNOISE
# <chr> <chr>           <chr> <chr>         <chr>         <chr> <chr>    <chr>    <chr> <chr>      <chr>    <chr>       <chr>     <chr>        <chr>
#   1 1     A1a1            false CM1           hsa-let-7a    UNKN~ FAM      NFQ-MGB  21.2~ 21.227075~ NA       1.32420543~ 0.781274~ N            N
# 2 2     A1a2            false CM1           hsa-let-7c    UNKN~ FAM      NFQ-MGB  27.7~ 27.738746~ NA       1.47680113~ 0.957805~ N            N
# 3 3     A1a3            false CM1           hsa-let-7f    UNKN~ FAM      NFQ-MGB  23.5~ 23.539237~ NA       1.39804202~ 0.967745~ N            N
# 4 5     A1a5            false CM1           hsa-miR-15a   UNKN~ FAM      NFQ-MGB  25.3~ 25.395061~ NA       1.44830555~ 0.962559~ N            N
# 5 6     A1a6            false CM1           hsa-miR-15b   UNKN~ FAM      NFQ-MGB  22.5~ 22.563753~ NA       1.54591603~ 0.957571~ N            N
# 6 7     A1a7            false CM1           hsa-miR-16    UNKN~ FAM      NFQ-MGB  19.2~ 19.231309~ NA       1.47410729~ 0.927107~ N            N

# -- Loading of data -- # END ...........................................................................................................................................

# -- Creating a list of patients -- BEGIN -------------------------------------------------------------------------------------------------------------------------------
tev_64[,4] -> patients_tev_64
tev_65[,4] -> patients_tev_65
tev_71[,4] -> patients_tev_71

vce_88[,4] -> patients_vce_88
vce_93[,4] -> patients_vce_93
vce_95[,4] -> patients_vce_95

ytk02[,4] -> patients_ytk02
ytk03[,4] -> patients_ytk03
ytk04[,4] -> patients_ytk04
ytk05[,4] -> patients_ytk05
ytk07[,4] -> patients_ytk07

#removing the duplicates
patient_list_tev_64=unique(patients_tev_64);patient_list_tev_64=as.data.frame(patient_list_tev_64)
patient_list_tev_65=unique(patients_tev_65);patient_list_tev_65=as.data.frame(patient_list_tev_65)
patient_list_tev_71=unique(patients_tev_71);patient_list_tev_71=as.data.frame(patient_list_tev_71)

patient_list_vce_88=unique(patients_vce_88);patient_list_vce_88=as.data.frame(patient_list_vce_88)
patient_list_vce_93=unique(patients_vce_93);patient_list_vce_93=as.data.frame(patient_list_vce_93)
patient_list_vce_95=unique(patients_vce_95);patient_list_vce_95=as.data.frame(patient_list_vce_95)

patient_list_ytk02=unique(patients_ytk02);patient_list_ytk02=as.data.frame(patient_list_ytk02)
patient_list_ytk03=unique(patients_ytk03);patient_list_ytk03=as.data.frame(patient_list_ytk03)
patient_list_ytk04=unique(patients_ytk04);patient_list_ytk04=as.data.frame(patient_list_ytk04)
patient_list_ytk05=unique(patients_ytk05);patient_list_ytk05=as.data.frame(patient_list_ytk05)
patient_list_ytk07=unique(patients_ytk07);patient_list_ytk07=as.data.frame(patient_list_ytk07)

rm(patients_tev_64, patients_tev_65, patients_tev_71, patients_vce_88, patients_vce_93,
   patients_vce_95, patients_ytk02, patients_ytk03, patients_ytk04, patients_ytk05,
   patients_ytk07)


# -- Creating individual files for each patient -- BEGIN ------------------------------------------------------------------

library(DataCombine) #this is to use grepl.sub

setwd("Individual_files/") #go to this directory to write the files

#Write individual files for tev_64
for (i in 1:length(t(patient_list_tev_64))){
  for (j in i) {
    temp=grepl.sub(tev_64, patient_list_tev_64[j,], 4) #This greps the patient in j in the whole dataset
    k=patient_list_tev_64[j,] #creates a vector with sample names instead of numbers
    write.table(temp, paste0(k,"_","tev_64_ctc",".txt"),sep="\t", quote=FALSE)  #this assigns the specific sample name to the dataset
  }
  rm(temp, j, i)
}

#Write individual files for tev_65
for (i in 1:length(t(patient_list_tev_65))){
  for (j in i) {
    temp=grepl.sub(tev_65, patient_list_tev_65[j,], 4) #This greps the patient in j in the whole dataset
    k=patient_list_tev_65[j,] #creates a vector with sample names instead of numbers
    write.table(temp, paste0(k,"_","tev_65_ctc",".txt"),sep="\t", quote=FALSE)  #this assigns the specific sample name to the dataset
  }
  rm(temp, j, i, k)
}

#Write individual files for tev_71
for (i in 1:length(t(patient_list_tev_71))){
  for (j in i) {
    temp=grepl.sub(tev_71, patient_list_tev_71[j,], 4) #This greps the patient in j in the whole dataset
    k=patient_list_tev_71[j,] #creates a vector with sample names instead of numbers
    write.table(temp, paste0(k,"_","tev_71_ctc",".txt"),sep="\t", quote=FALSE)  #this assigns the specific sample name to the dataset
  }
  rm(temp, j, i, k)
}

#Write individual files for vce_88
for (i in 1:length(t(patient_list_vce_88))){
  for (j in i) {
    temp=grepl.sub(vce_88, patient_list_vce_88[j,], 4) #This greps the patient in j in the whole dataset
    k=patient_list_vce_88[j,] #creates a vector with sample names instead of numbers
    write.table(temp, paste0(k,"_","vce_88_ctc",".txt"),sep="\t", quote=FALSE)  #this assigns the specific sample name to the dataset
  }
  rm(temp, j, i, k)
}

#Write individual files for vce_93
for (i in 1:length(t(patient_list_vce_93))){
  for (j in i) {
    temp=grepl.sub(vce_93, patient_list_vce_93[j,], 4) #This greps the patient in j in the whole dataset
    k=patient_list_vce_93[j,] #creates a vector with sample names instead of numbers
    write.table(temp, paste0(k,"_","vce_93_ctc",".txt"),sep="\t", quote=FALSE)  #this assigns the specific sample name to the dataset
  }
  rm(temp, j, i, k)
}

#Write individual files for vce_95
for (i in 1:length(t(patient_list_vce_95))){
  for (j in i) {
    temp=grepl.sub(vce_95, patient_list_vce_95[j,], 4) #This greps the patient in j in the whole dataset
    k=patient_list_vce_95[j,] #creates a vector with sample names instead of numbers
    write.table(temp, paste0(k,"_","vce_95_ctc",".txt"), sep="\t",quote=FALSE)  #this assigns the specific sample name to the dataset
  }
  rm(temp, j, i, k)
}

#Write individual files for ykt02
for (i in 1:length(t(patient_list_ytk02))){
  for (j in i) {
    temp=grepl.sub(ytk02, patient_list_ytk02[j,], 4) #This greps the patient in j in the whole dataset
    k=patient_list_ytk02[j,] #creates a vector with sample names instead of numbers
    write.table(temp, paste0(k,"_","ytk02_ctc",".txt"),sep="\t", quote=FALSE)  #this assigns the specific sample name to the dataset
  }
  rm(temp, j, i, k)
}

#Write individual files for ykt03
for (i in 1:length(t(patient_list_ytk03))){
  for (j in i) {
    temp=grepl.sub(ytk03, patient_list_ytk03[j,], 4) #This greps the patient in j in the whole dataset
    k=patient_list_ytk03[j,] #creates a vector with sample names instead of numbers
    write.table(temp, paste0(k,"_","ytk03_ctc",".txt"), sep="\t",quote = FALSE)  #this assigns the specific sample name to the dataset
  }
  rm(temp, j, i, k)
}

#Write individual files for ykt04
for (i in 1:length(t(patient_list_ytk04))){
  for (j in i) {
    temp=grepl.sub(ytk04, patient_list_ytk04[j,], 4) #This greps the patient in j in the whole dataset
    k=patient_list_ytk04[j,] #creates a vector with sample names instead of numbers
    write.table(temp, paste0(k,"_","ytk04_ctc",".txt"),sep="\t",quote = FALSE)  #this assigns the specific sample name to the dataset
  }
  rm(temp, j, i, k)
}

#Write individual files for ykt05
for (i in 1:length(t(patient_list_ytk05))){
  for (j in i) {
    temp=grepl.sub(ytk05, patient_list_ytk05[j,], 4) #This greps the patient in j in the whole dataset
    k=patient_list_ytk05[j,] #creates a vector with sample names instead of numbers
    write.table(temp, paste0(k,"_","ytk05_ctc",".txt"),sep="\t",quote = FALSE)  #this assigns the specific sample name to the dataset
  }
  rm(temp, j, i, k)
}

#Write individual files for ykt05
for (i in 1:length(t(patient_list_ytk07))){
  for (j in i) {
    temp=grepl.sub(ytk07, patient_list_ytk07[j,], 4) #This greps the patient in j in the whole dataset
    k=patient_list_ytk07[j,] #creates a vector with sample names instead of numbers
    write.table(temp, paste0(k,"_","ytk07_ctc",".txt"),sep="\t",quote = FALSE)  #this assigns the specific sample name to the dataset
  }
  rm(temp, j, i, k)
}

rm(tev_64, tev_65, tev_71, vce_88 , vce_93, vce_95, ytk02, ytk03, ytk04, ytk05, ytk07)
rm(patient_list_tev_64, patient_list_tev_65, patient_list_tev_71,patient_list_vce_88, patient_list_vce_93, 
   patient_list_vce_95, patient_list_ytk02, patient_list_ytk03,patient_list_ytk04, 
   patient_list_ytk05, patient_list_ytk07)

# -- Creating individual files for each patient -- END ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#----------------------------------------------------------------- SECOND PART --------------------------------------------------------------------------------------------------------------------------------------

                                # BEGIN HERE FOR THE HTqPCR EXPLORATORY AND DIFFERENTIAL EXPRESSION ANALYSIS 
                                                                      # |
                                                                      # |
                                                                      # V

#.........................................................................................................................
# HTqPCR PACKAGE -- # BEGIN -------------------------------------------------------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("HTqPCR")

#..................................................................................................................................
# Data levels: All the relevant qPCR objects that are being created below are explained here: 

# raw_data -> raw data with no changes
# data_set -> setting the categories "Unrealiable" and "Undetermined" on the data (does not change the raw data)
# filtered_data -> sets all Unreliable/Undetermined to NA
# exprs -> all the NAs from filtered data are set to 0
# q.norm_raw_data -> raw data was normalised using quantile normalisation
# q.norm_data_filter_out -> filtered data with NA's was normalised using quantile normalisation
# q.norm_data_filter_exprs -> filtered data with NA's set to 0 was normalised with quantile normalisation

# All relevant analyses were done with quantile normalisation because comparing to other normalisation methods
# this one is the best.

#.....................................................................................................................................

library(HTqPCR)
library(genefilter)
library(dplyr)
library(stringr)
path = "Patients_files/" #PATH TO FILES
files = read.delim(file.path(path, "files.txt")) 

# READ QPCR OBJECT ...................................................................................................................
raw_data = readCtData(files=files$File, path=path, n.features = 112, na.value = 40, header = TRUE, 
                      column.info = list(flag=3,feature = 5, type = 6, Ct = 10, position = 2), sep = "\t")

data_set = setCategory(raw_data, Ct.max=38, Ct.min=14,quantile=0.9, flag=TRUE, flag.out="true",
                       groups=files$Treatment)

filtered_data<-filterCategory(data_set) #set all Unreliable/Undetermined to NA
filtered_data_out=filterCtData(data_set, remove.category = c("Unreliable","Undetermined"), n.category = 74)

# Quantile normalization - Chosen normalization for plots ...................................................................................................................


q.norm_raw_data=normalizeCtData(raw_data, norm="quantile")
plotCtDensity(q.norm_raw_data, legend=FALSE, main="Quantile normalisation - Raw data")
data_filter_2=filterCategory(data_set,na.categories = c("Unreliable","Undetermined"))
data_filter_3 <- filterCtData(data_filter_2, remove.category=c("Unreliable","Undetermined"),
                              n.category=100,remove.IQR = 1)
plotCtDensity(data_filter_3, legend=FALSE, main="filtered_data")
qFilt.q.norm_filter=normalizeCtData(data_filter_3, norm="quantile")
plotCtDensity(qFilt.q.norm_filter,legend = FALSE, main="filtered_data_quantile")
exprs(filtered_data)[is.na(exprs(filtered_data))] <- 0
q.norm_data_filter_exprs=normalizeCtData(filtered_data, norm = "quantile")
plotCtDensity(q.norm_data_filter_exprs, legend=FALSE, main="Quantile normalisation - Filtered data=0")



# --> Differential expression between Baseline, During Treatment and After treatment -----------------------------------------------------------------------------------------------

files= read.delim(file.path(path, "files.txt")) 
batch=factor(files$Batches)
cond=factor(files$Treatment, levels=c("baseline","control","during","xafter"))

design_batches=model.matrix(~0+cond_1+batch)

colnames(design_batches) = c("baseline","control","during","xafter","batch2","batch3")

contrasts_batches= makeContrasts(during-baseline, xafter-during, xafter-baseline, 
                                   baseline-control, 
                                   levels=design_batches)

q.norm_data_filter_exprs=normalizeCtData(filtered_data, norm = "quantile")
plotCtDensity(q.norm_data_filter_exprs, legend=FALSE, main="Quantile normalisation - Filtered data=0")

exprs(filtered_data)[is.na(exprs(filtered_data))]=0

qnorm_exprs_diff.exp_3b = limmaCtData(q.norm_data_filter_exprs, design=design_batches,
                                      contrasts=contrasts_batches,sort=FALSE)

qnorm_exprs_diff.exp_3b$`during - baseline`$genes=featureNames(raw_data)
qnorm_exprs_diff.exp_3b$`xafter - during`$genes=featureNames(raw_data)
qnorm_exprs_diff.exp_3b$`xafter - baseline`$genes=featureNames(raw_data)
qnorm_exprs_diff.exp_3b$`baseline - control`$genes=featureNames(raw_data)


# --> Differential expression between good and bad responders  -----------------------------------------------------------------------------------------------

# Baseline - Responders vs Non-responders
path = "LARC_DATA/patient_files_good_bad_responders/" #PATH TO FILES
files = read.delim(file.path(path, "files_responders_baseline.txt")) 

# READ QPCR OBJECT ...................................................................................................................
raw_data = readCtData(files=files$File, path=path, n.features = 112, na.value = 40, header = TRUE, 
                      column.info = list(flag=3,feature = 5, type = 6, Ct = 10, position = 2), sep = "\t")

data_set = setCategory(raw_data, Ct.max=38, Ct.min=14,quantile=0.9, flag=TRUE, flag.out="true",
                       groups=files$Treatment)

file= read.delim(file.path(path, "files_responders_baseline.txt")) 
batch=factor(file$Batches)
condition=factor(file$Treatment, levels=c("GR","BR"))


design_matrix=model.matrix(~0+condition+batch)


colnames(design_matrix) = c("good","bad","batch2","batch3")

contrasts= makeContrasts(good-bad, levels=design_matrix)

par(mfrow=c(1,1))
filtered_data<-filterCategory(data_set) #set all Unreliable/Undetermined to NA
q.norm_data_filter_out=normalizeCtData(raw_data, norm = "quantile")

plotCtDensity(q.norm_data_filter_out, legend=FALSE, main="Quantile normalisation - Filtered data")
qnorm_filter_out_diff_baseline=limmaCtData(raw_data, design=design_matrix,
                                           contrasts = contrasts,
                                           sort = FALSE) 
qnorm_filter_out_diff_baseline$`good-bad`$genes=featureNames(raw_data)


#During treatment - Responders vs Non-responders

path = "f:/PhD/PhD_2019/LARC_DATA/patient_files_good_bad_responders/" #PATH TO FILES
files = read.delim(file.path(path, "files_responders_during.txt")) 

# READ QPCR OBJECT ...................................................................................................................
raw_data = readCtData(files=files$File, path=path, n.features = 112, na.value = 40, header = TRUE, 
                      column.info = list(flag=3,feature = 5, type = 6, Ct = 10, position = 2), sep = "\t")

data_set = setCategory(raw_data, Ct.max=38, Ct.min=14,quantile=0.9, flag=TRUE, flag.out="true",
                       groups=files$Treatment)

file= read.delim(file.path(path, "files_responders_during.txt")) 
batch=factor(file$Batches)
condition=factor(file$Treatment, levels=c("GR","BR"))

design_matrix=model.matrix(~0+condition+batch)


colnames(design_matrix) = c("good","bad","batch2","batch3")

contrasts= makeContrasts(good-bad, levels=design_matrix)

par(mfrow=c(1,1))
filtered_data<-filterCategory(data_set) #set all Unreliable/Undetermined to NA
q.norm_data_filter_out=normalizeCtData(raw_data, norm = "quantile")

plotCtDensity(q.norm_data_filter_out, legend=FALSE, main="Quantile normalisation - Filtered data")
qnorm_filter_out_diff_during=limmaCtData(raw_data, design=design_matrix,
                                         contrasts = contrasts,
                                         sort = FALSE) 
qnorm_filter_out_diff_during$`good-bad`$genes=featureNames(raw_data)
 
#Post-treatment - Responders vs Non-responders
path = "f:/PhD/PhD_2019/LARC_DATA/patient_files_good_bad_responders/" #PATH TO FILES
files = read.delim(file.path(path, "files_responders_after.txt")) 

# READ QPCR OBJECT ...................................................................................................................
raw_data = readCtData(files=files$File, path=path, n.features = 112, na.value = 40, header = TRUE, 
                      column.info = list(flag=3,feature = 5, type = 6, Ct = 10, position = 2), sep = "\t")

data_set = setCategory(raw_data, Ct.max=38, Ct.min=14,quantile=0.9, flag=TRUE, flag.out="true",
                       groups=files$Treatment)

file= read.delim(file.path(path, "files_responders_after.txt")) 
batch=factor(file$Batches)
condition=factor(file$Treatment, levels=c("GR","BR"))

design_matrix=model.matrix(~0+condition+batch)


colnames(design_matrix) = c("good","bad","batch2","batch3")

contrasts= makeContrasts(good-bad, levels=design_matrix)

par(mfrow=c(1,1))
filtered_data<-filterCategory(data_set) #set all Unreliable/Undetermined to NA
q.norm_data_filter_out=normalizeCtData(raw_data, norm = "quantile")

plotCtDensity(q.norm_data_filter_out, legend=FALSE, main="Quantile normalisation - Filtered data")
qnorm_filter_out_diff._after=limmaCtData(raw_data, design=design_matrix,
                                         contrasts = contrasts,
                                         sort = FALSE)
qnorm_filter_out_diff._after$`good-bad`$genes=featureNames(raw_data)


# --> PLOT HEATMAPS --------------------------------------------------------------------------------------------------------------------

sample_names=colnames(raw_data@assayData$featureCategory)
a=gsub("_[0-9]*", "_", sample_names)
names=a %>% stringr::str_remove(c("__ctc"))
names=names%>% stringr::str_remove(c("_ctc"))
names=gsub("tev","TEV", names)
names=gsub("ytk","YTK", names)
names=gsub("vce","VCE", names)

filtered_data<-filterCategory(data_set) 
exprs(filtered_data)[is.na(exprs(filtered_data))]=1

# Unormalized data 
plotCtHeatmap(filtered_data,sample.names = names, dist="euclidean", key=FALSE, offsetRow = 1,
              cexRow = .1, cexCol=1.1, Colv=TRUE, dendrogram="col", mar = c(8,2), gene.names = "",
              lwid = c(0.1,7),
              lhei = c(0.1,0.5))
  
# quantile normalized data
plotCtHeatmap(q.norm_raw_data, dendrogram="col", sample.names = names, key=FALSE, 
              cexRow = .1, cexCol=1.1,lwid = c(0.1,7), gene.names = "",mar = c(8,2),
              lhei = c(0.1,0.5))



