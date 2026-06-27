library(tidyverse)
library(GEOquery)
library(ggplot2)
library(affy)
library()

#getting sup files
gse <- getGEOSuppFiles('GSE148537')


#unzipping .tar file
untar('GSE148537/GSE148537_RAW.tar', exdir = 'data/')

#reading the raw data files
raw_data <- ReadAffy(celfile.path = 'data/')

#normalizing the raw_data
normalized_data <- rma(raw_data)

#get expression estimates
normalized_exprs <- as.data.frame(exprs(normalized_data))

#get gse matrix file
gse <- getGEO('GSE148537', GSEMatrix = TRUE)

#fetch feature data to get IDs - gene mapping symbols
feature_data <- gse$GSE148537_series_matrix.txt.gz@featureData@data

#choosing our columns of interest and making a subset out of it 
subset <- feature_data[,c(1, 11)]

#joining two dataframe into one using their IDs
df <- normalized_exprs%>%
  rownames_to_column(var = 'ID')%>%
  inner_join(., subset, by = 'ID')


