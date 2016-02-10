######################################################################
# methyLiftover - A utility to map whole-genome & reduced representation bisulfite data by genomic position to the Illumina
#     450k Methylation Array. 
# Authors: Alexander J. Titus, E. Andres Houseman, Kevin C. Johnson, Brock C. Christensen
######################################################################

######################################################################
# Install and load dependency packages
######################################################################
setUpWorkSpace <- function() {
  # Load the packages needed for the processing
  if("data.table" %in% rownames(installed.packages()) == FALSE) {install.packages('data.table')}
  if("IlluminaHumanMethylation450kanno.ilmn12.hg19" %in% rownames(installed.packages()) == FALSE) {
    biocLite("IlluminaHumanMethylation450kanno.ilmn12.hg19")}
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(data.table)
}

######################################################################
# liftover450k: Liftover genomic data to the 450k probset CpG set. 
######################################################################
# @param fileDir absolute directory path to the input file, including file name
# @param fileKeys vector containing two column names (class character) in order of c(chromosomeName, genomicPosition) EX/ c('chr', 'pos')
# @param fileSkipLine number of lines to skip before starting input. The default of 1 runs a standard data.table import skipping BED file headers
# @param outputfileDir the directory to save the liftover RData file
# @param outputFileName the name of the resulting output file of the liftover WITHOUT file extension
liftover450k <- function(fileDir, fileKeys, fileSkipLine = 1, outputFileDir, outputFileName){
  
  # Load the Illumina 450k Methylation Array for merging
  illumData <- IlluminaHumanMethylation450kanno.ilmn12.hg19@data
  illumDataFrame <- illumData$Locations
  illumDataFrame <- cbind(illumDataFrame, illumData$Islands.UCSC )
  illumDataFrame$Type <- illumData$Manifest$Type
  illumDataFrame$probe <- rownames(illumDataFrame)
  illumDataTable <- data.table(as.data.frame(illumDataFrame))
  
  inputFile <- fread(fileDir, skip = fileSkipLine)
  WGBSliftover <- merge(illumDataTable, inputFile, by.x=c('chr', 'pos'), by.y=fileKeys)
  rownames(WGBSliftover) <- WGBSliftover$probe
  setwd(outputFileDir)
  fileName <- paste(outputFileName, '.RData', sep='')
  save(WGBSliftover, file = fileName)
}

######################################################################
# liftoverUserFiles: Liftover genomic data and merge two datasets.
######################################################################
# @param fileDir1 absolute directory path to file 1, including file name
# @param fileDir2 absolute directory path to file 1, including file name
# @param file1Keys vector containing two column names (class character) in order of c(chromosomeName, genomicPosition) EX/ c('chr', 'pos')
# @param file2Keys vector containing two column names (class character) in order of c(chromosomeName, genomicPosition) EX/ c('chr', 'pos')
# @param file1SkipLine number of lines to skip before starting file 1 import The default of 1 runs a standard data.table import skipping the BED header row
# @param file2SkipLine number of lines to skip before starting file 2 import The default of 1 runs a standard data.table import skipping the BED header row
# @param outputfileDir the directory to save the liftover RData file
# @param outputFileName the name of the resulting output file of the liftover WITHOUT file extension
liftoverUserFiles <- function(fileDir1, fileDir2, file1Keys, file2Keys, file1SkipLine = 1, file2SkipLine = 1, outputFileDir, outputFileName){
  
  # Read in our files
  file1 <- fread(fileDir1, skip=file1SkipLine)
  file2 <- fread(fileDir2, skip=file2SkipLine)
  
  # Create our merged/mapped file
  liftOver <- merge(file1, file2, by.x=file1Keys, by.y=file2Keys)
  
  # Save the output file
  setwd(outputFileDir)
  fileName <- paste(outputFileName, '.RData', sep='')
  save(liftOver, file = fileName)
}

######################################################################
# Liftover genomic data from files within a directory to the 450k probset CpG set.
######################################################################
# @param fileDir directory where files are stored
# @param fileKeys vector containing two column names (class character) in order of c(chromosomeName, genomicPosition) EX/ c('chr', 'pos')
# @param fileSkipLine number of lines to skip before starting input. The default of 1 runs a standard data.table import skipping the BED header row
# @param outputfileDir the directory to save the liftover RData file
# @param outputFileName the name of the resulting output file of the liftover WITHOUT file extension
liftover450kset <- function(fileDir, fileKeys, fileSkipLine = 1, outputFileDir, outputFileName){
  # Set the directory where the BED files are stored
  setwd(fileDir)
  
  # List all BED files in a directory
  files <- list.files()
  for(i in 1:length(files)){
    outputFileNameUse <- paste(outputFileName, i, sep='_')
    fileDirUse <- paste(fileDir, files[i], sep='/')
    liftover450k(fileDirUse, fileKeys, fileSkipLine, outputFileDir, outputFileNameUse)
  }
}




