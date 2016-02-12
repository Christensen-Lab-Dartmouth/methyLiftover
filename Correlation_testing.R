correlationTest <- function(files, corMeth){
  corrComparison <- data.frame()
  for(i in 1:length(files)){
    list <- c('I', 'II', 'Island', 'Shelf', 'Shore', 'OpenSea')
    load(files[i])
    data <- TCGA_to_WGBS_complete
    rm(TCGA_to_WGBS_complete)
    data[data$Relation_to_island=="N_Shore", ]$Relation_to_island <- 'Shore'
    data[data$Relation_to_island=="S_Shore", ]$Relation_to_island <- 'Shore'
    data[data$Relation_to_island=="N_Shelf", ]$Relation_to_island <- 'Shelf'
    data[data$Relation_to_island=="S_Shelf", ]$Relation_to_island <- 'Shelf'
    data$Beta_value100 <- data$Beta_value*100
    data$percentMeth <- as.numeric(data$percentMeth)
    cols <- c('percentMeth', 'Beta_value100')
    
    for(k in 1:15){
      dataTemp <- data[data$numCTreads >= k,]
      # Calculation correlation
      corrs <- cor(dataTemp[,cols], use = 'pairwise.complete.obs', method = corMeth)[1,2]
      # Build our data.frame
      file <- files[i]
      tcgaCode <- strsplit(file, '_')[[1]][1]
      disease <- strsplit(file, '_')[[1]][2]
      sampleType <- strsplit(tcgaCode, '-')[[1]][4]
      corrFrame <- data.frame(tcgaCode)
      corrFrame$disease <- disease
      corrFrame$sampleType <- sampleType
      corrFrame$correlation <- corrs
      corrFrame$subset <- 'All'
      corrFrame$minReadDepth <- k
      corrFrame$probeType <- 'All'
      corrComparison <- rbind(corrFrame, corrComparison)
      
      for(j in 1:2){
        dataTemp <- data[data$Probe_type==list[j], ]
        dataTemp <- dataTemp[dataTemp$numCTreads >= k,]
        # Calculation correlation
        corrs <- cor(dataTemp[,cols], use = 'pairwise.complete.obs', method = corMeth)[1,2]
        # Build our data.frame
        file <- files[i]
        tcgaCode <- strsplit(file, '_')[[1]][1]
        disease <- strsplit(file, '_')[[1]][2]
        sampleType <- strsplit(tcgaCode, '-')[[1]][4]
        corrFrame <- data.frame(tcgaCode)
        corrFrame$disease <- disease
        corrFrame$sampleType <- sampleType
        corrFrame$correlation <- corrs
        corrFrame$subset <- list[j]
        corrFrame$minReadDepth <- k
        corrFrame$probeType <- list[j]
        corrComparison <- rbind(corrFrame, corrComparison)
        
        for(n in 3:length(list)){
          dataTemp <- data[data$Probe_type==list[j], ]
          dataTemp <- dataTemp[dataTemp$Relation_to_island==list[n], ]
          dataTemp <- dataTemp[dataTemp$numCTreads >= k,]
          # Calculation correlation
          corrs <- cor(dataTemp[,cols], use = 'pairwise.complete.obs', method = corMeth)[1,2]
          # Build our data.frame
          file <- files[i]
          tcgaCode <- strsplit(file, '_')[[1]][1]
          disease <- strsplit(file, '_')[[1]][2]
          sampleType <- strsplit(tcgaCode, '-')[[1]][4]
          corrFrame <- data.frame(tcgaCode)
          corrFrame$disease <- disease
          corrFrame$sampleType <- sampleType
          corrFrame$correlation <- corrs
          corrFrame$subset <- list[n]
          corrFrame$minReadDepth <- k
          corrFrame$probeType <- list[j]
          corrComparison <- rbind(corrFrame, corrComparison)
        }
      }
    }
  }
  return(corrComparison)
}
##########
# Correlation testing
##########
# BLCA correlation testing
setwd('/Users/alexandertitus/Desktop/WGBSliftover_TCGA/Analysis/liftoverSamples/BMIQ_normalized/BLCA')
fileList <- list.files()
corrs1 <- correlationTest(fileList, 'pearson')

# BRCA correlation testing
setwd('/Users/alexandertitus/Desktop/WGBSliftover_TCGA/Analysis/liftoverSamples/BMIQ_normalized/BRCA')
fileList <- list.files()
corrs2 <- correlationTest(fileList, 'pearson')

# UCEC correlation testing
setwd('/Users/alexandertitus/Desktop/WGBSliftover_TCGA/Analysis/liftoverSamples/BMIQ_normalized/UCEC')
fileList <- list.files()
corrs3 <- correlationTest(fileList, 'pearson')

# LUAD correlation testing
setwd('/Users/alexandertitus/Desktop/WGBSliftover_TCGA/Analysis/liftoverSamples/BMIQ_normalized/LUAD')
fileList <- list.files()
corrs4 <- correlationTest(fileList, 'pearson')

# STAD correlation testing
setwd('/Users/alexandertitus/Desktop/WGBSliftover_TCGA/Analysis/liftoverSamples/BMIQ_normalized/STAD')
fileList <- list.files()
corrs5 <- correlationTest(fileList, 'pearson')

correlationMatrix <- rbind(corrs1, corrs2, corrs3, corrs4, corrs5)

setwd('/Users/alexandertitus/Desktop/WGBSliftover_TCGA/Analysis/liftoverSamples/BMIQ_normalized')
write.csv(correlationMatrix, 'TCGA_methyLiftover_correlation_testing_with_probeTypes_pearson_BMIQnorm.csv')

##########
# MAD calculation across samples
##########
MADacrossSets <- function(dir, probeType){
  library('data.table')
  setwd(dir)
  files <- list.files(dir)
  load(files[1])
  absDiffDataFrame <- data.frame(rownames(TCGA_to_WGBS_complete))
  colnames(absDiffDataFrame) <- 'probes'
  rownames(absDiffDataFrame) <- absDiffDataFrame[,1]
  
  for(i in 1:length(files)){ 
    load(files[i])
    
    if(probeType == 'I'){
      TCGA_to_WGBS_complete <- TCGA_to_WGBS_complete[TCGA_to_WGBS_complete$Probe_type=='I',]
    } else if(probeType == 'II') {
      TCGA_to_WGBS_complete <- TCGA_to_WGBS_complete[TCGA_to_WGBS_complete$Probe_type=='II',]
    } else {
      TCGA_to_WGBS_complete <- TCGA_to_WGBS_complete
    }
    
    tempDataFrame <- as.data.frame(TCGA_to_WGBS_complete[,c('percentMeth', 'Beta_value')])
    tempDataFrame$probes <- rownames(tempDataFrame)
    tempDataFrame$percentMeth <- as.numeric(tempDataFrame$percentMeth)
    tempDataFrame$Beta_value <- as.numeric(tempDataFrame$Beta_value)
    tempDataFrame$Beta_valueX100 <- as.numeric(tempDataFrame$Beta_value*100)
    tempDataFrame$Beta_value <- NULL
    label <- paste('absDiff-sample#', substr(files[i], start = 1, stop = 20), sep='-')
    tempDataFrame[[label]] <- abs(tempDataFrame$percentMeth-tempDataFrame$Beta_valueX100)
    
    tempDataFrame$percentMeth <- NULL
    tempDataFrame$Beta_valueX100 <- NULL
    
    absDiffDataFrame <- merge(absDiffDataFrame, tempDataFrame, by='probes')
    rownames(absDiffDataFrame) <- absDiffDataFrame[,1]
  }
  
  MADacross <- data.frame(NA)
  for(j in 2:ncol(absDiffDataFrame)){
    MADacross[j-1,] <- mad(absDiffDataFrame[,j], na.rm = T)
    rownames(MADacross)[j-1] <- colnames(absDiffDataFrame)[j]
  }
  print(mean(MADacross[,1]))
  print(median(MADacross[,1]))
  print(range(MADacross))
  return(MADacross)
}

##########
# MAD calculation within samples
##########
MADwithinSet <- function(dir, probeType) {
  library('data.table')
  setwd(dir)
  files <- list.files(dir)
  MADwithin <- data.frame(NA)
  MADwithin[1,2] <- NA
  colnames(MADwithin) <- c('PercentMethMAD', 'BetaMAD')
  for(i in 1:length(files)){
    load(files[i])
    if(probeType == 'I'){
      TCGA_to_WGBS_complete <- TCGA_to_WGBS_complete[TCGA_to_WGBS_complete$Probe_type=='I',]
    } else if(probeType == 'II') {
      TCGA_to_WGBS_complete <- TCGA_to_WGBS_complete[TCGA_to_WGBS_complete$Probe_type=='II',]
    } else {
      TCGA_to_WGBS_complete <- TCGA_to_WGBS_complete
    }
    
    MADwithin[i,1] <- mad(as.numeric(TCGA_to_WGBS_complete$percentMeth), na.rm = T)
    MADwithin[i,2] <- mad(as.numeric(TCGA_to_WGBS_complete$Beta_value), na.rm = T)
    rownames(MADwithin)[i] <- files[i] 
  }
  return(MADwithin)
}

##########
# MAD testing across all samples
##########
# dir <- '/Users/alexandertitus/Desktop/WGBSliftover_TCGA/Analysis/liftoverSamples'
dir <- '/Users/alexandertitus/Desktop/WGBSliftover_TCGA/Analysis/liftoverSamples/BMIQ_normalized'
probeTypeI <- 'I'
MADacrossI <- MADacrossSets(dir, probeTypeI)
MADwithinI <- MADwithinSet(dir, probeTypeI)

probeTypeII <- 'II'
MADacrossII <- MADacrossSets(dir, probeTypeII)
MADwithinII <- MADwithinSet(dir, probeTypeII)

probeType <- 'All'
MADacross <- MADacrossSets(dir, probeType)
MADwithin <- MADwithinSet(dir, probeType)

combinedReporting <- cbind(MADacross, MADacrossI, MADacrossII)
colnames(combinedReporting) <- c('All', 'Type_I', 'Type_II')
write.csv(combinedReporting, 'MADs_across_sets.csv')

##########
# MAD testing across normal samples
##########
# dir <- '/Users/alexandertitus/Desktop/WGBSliftover_TCGA/Analysis/liftoverSamples/Normal'
dir <- '/Users/alexandertitus/Desktop/WGBSliftover_TCGA/Analysis/liftoverSamples/BMIQ_normalized/Normal'
probeTypeI <- 'I'
MADacrossI <- MADacrossSets(dir, probeTypeI)
MADwithinI <- MADwithinSet(dir, probeTypeI)

probeTypeII <- 'II'
MADacrossII <- MADacrossSets(dir, probeTypeII)
MADwithinII <- MADwithinSet(dir, probeTypeII)

probeType <- 'All'
MADacross <- MADacrossSets(dir, probeType)
MADwithin <- MADwithinSet(dir, probeType)

combinedReporting <- cbind(MADacross, MADacrossI, MADacrossII)
colnames(combinedReporting) <- c('All', 'Type_I', 'Type_II')
write.csv(combinedReporting, 'MADs_across_sets_normal.csv')

##########
# MAD testing across tumor samples
##########
# dir <- '/Users/alexandertitus/Desktop/WGBSliftover_TCGA/Analysis/liftoverSamples/Tumor'
dir <- '/Users/alexandertitus/Desktop/WGBSliftover_TCGA/Analysis/liftoverSamples/BMIQ_normalized/Tumor'
probeTypeI <- 'I'
MADacrossI <- MADacrossSets(dir, probeTypeI)
MADwithinI <- MADwithinSet(dir, probeTypeI)

probeTypeII <- 'II'
MADacrossII <- MADacrossSets(dir, probeTypeII)
MADwithinII <- MADwithinSet(dir, probeTypeII)

probeType <- 'All'
MADacross <- MADacrossSets(dir, probeType)
MADwithin <- MADwithinSet(dir, probeType)

combinedReporting <- cbind(MADacross, MADacrossI, MADacrossII)
colnames(combinedReporting) <- c('All', 'Type_I', 'Type_II')
write.csv(combinedReporting, 'MADs_across_sets_tumor.csv')

##########
# IQR calculations across sets
##########
absDiffAcrossSets <- function(dir, probeType){
  library('data.table')
  setwd(dir)
  files <- list.files(dir)
  load(files[1])
  absDiffDataFrame <- data.frame(rownames(TCGA_to_WGBS_complete))
  colnames(absDiffDataFrame) <- 'probes'
  rownames(absDiffDataFrame) <- absDiffDataFrame[,1]
  
  for(i in 1:length(files)){ 
    load(files[i])
    
    if(probeType == 'I'){
      TCGA_to_WGBS_complete <- TCGA_to_WGBS_complete[TCGA_to_WGBS_complete$Probe_type=='I',]
    } else if(probeType == 'II') {
      TCGA_to_WGBS_complete <- TCGA_to_WGBS_complete[TCGA_to_WGBS_complete$Probe_type=='II',]
    } else {
      TCGA_to_WGBS_complete <- TCGA_to_WGBS_complete
    }
    
    tempDataFrame <- as.data.frame(TCGA_to_WGBS_complete[,c('percentMeth', 'Beta_value')])
    tempDataFrame$probes <- rownames(tempDataFrame)
    tempDataFrame$percentMeth <- as.numeric(tempDataFrame$percentMeth)
    tempDataFrame$Beta_value <- as.numeric(tempDataFrame$Beta_value)
    tempDataFrame$Beta_valueX100 <- as.numeric(tempDataFrame$Beta_value*100)
    tempDataFrame$Beta_value <- NULL
    label <- paste('absDiff-sample#', substr(files[i], start = 1, stop = 20), sep='-')
    tempDataFrame[[label]] <- abs(tempDataFrame$percentMeth-tempDataFrame$Beta_valueX100)
    
    tempDataFrame$percentMeth <- NULL
    tempDataFrame$Beta_valueX100 <- NULL
    
    absDiffDataFrame <- merge(absDiffDataFrame, tempDataFrame, by='probes')
    rownames(absDiffDataFrame) <- absDiffDataFrame[,1]
  }
  
  return(absDiffDataFrame)
}
# IQR calculations
setwd('/Users/alexandertitus/Desktop/WGBSliftover_TCGA/Analysis/liftoverSamples/BMIQ_normalized/All')
absDiff <- absDiffAcrossSets(getwd(), 'All')
high75 <- 0
low75 <- 100
high25 <- 0
low25 <- 100
for(l in 2:ncol(absDiff)){
  if(quantile(absDiff[,l], na.rm=T)[2] > high25) high25 <- quantile(absDiff[,l], na.rm=T)[2]
  if(quantile(absDiff[,l], na.rm=T)[2] < low25) low25 <- quantile(absDiff[,l], na.rm=T)[2]
  if(quantile(absDiff[,l], na.rm=T)[4] > high75) high75 <- quantile(absDiff[,l], na.rm=T)[4]
  if(quantile(absDiff[,l], na.rm=T)[4] < low75) low75 <- quantile(absDiff[,l], na.rm=T)[4]
}

high75-low75
high25-low25
