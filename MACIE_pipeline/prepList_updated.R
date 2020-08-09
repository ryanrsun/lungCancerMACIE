#--------------------------------------------------------#
# Input a list of SNPs with five columns: RS, Chr, BP, Ref, Alt
# Retrive the annotations from CADD database, note which ones cannot
# be retrieved or have different BP/Ref/Alt, format nicely.
# Then run the python script to get the posterior probabilities.
#--------------------------------------------------------#

library(data.table)
library(dplyr)
library(magrittr)

# Change these inputs for your filesystems
inputDir <- "./snpLists"
paramDir <- "./parameters"
caddDir <- "./caddDir"
outputDir <- "./output"
inputFname <- "exList.txt"
outputFnameRoot <- "exList"
# no more changes
#--------------------------------------------------------#

colUsed <- c("GC","CpG","bStatistic","fitCons","cHmmTssA","cHmmTssAFlnk","cHmmTxFlnk","cHmmTx",       
	"cHmmTxWk","cHmmEnhG","cHmmEnh","cHmmZnfRpts","cHmmHet","cHmmTssBiv","cHmmBivFlnk","cHmmEnhBiv",    
	"cHmmReprPC","cHmmReprPCWk","cHmmQuies","EncExp","EncH3K27Ac","EncH3K4Me1","EncH3K4Me3","EncNucleo",     
  "EncOCCombPVal","EncOCDNasePVal","EncOCFairePVal","EncOCpolIIPVal",
	"EncOCctcfPVal","EncOCmycPVal","EncOCDNaseSig","EncOCFaireSig", 
	"EncOCpolIISig","EncOCctcfSig","EncOCmycSig","TFBS",          
	"TFBSPeaks","TFBSPeaksMax","minDistTSS","minDistTSE")
imputeZero <- c("bStatistic","fitCons","EncExp","EncH3K27Ac","EncH3K4Me1","EncH3K4Me3","EncNucleo",
  "EncOCCombPVal","EncOCDNasePVal","EncOCFairePVal","EncOCpolIIPVal",
  "EncOCctcfPVal","EncOCmycPVal","EncOCDNaseSig","EncOCFaireSig",
  "EncOCpolIISig","EncOCctcfSig","EncOCmycSig","TFBS",
  "TFBSPeaks","TFBSPeaksMax")
# needs to be in this order for the rest of the pipeline!
finalEpi <- c("GC", "CpG", "EncExp", "EncH3K27Ac", "EncH3K4Me1", "EncH3K4Me3",   
	"EncOCCombPVal", "EncOCDNasePVal", "EncOCFairePVal",
	"EncOCpolIIPVal","EncOCctcfPVal", "EncOCmycPVal", "EncOCDNaseSig", "EncOCFaireSig","EncOCpolIISig", "EncOCctcfSig","EncOCmycSig",   
	"TFBS", "TFBSPeaks", "TFBSPeaksMax","cHmmTSSOdds","cHmmTxOdds", "cHmmEnhOdds", "cHmmZnfOdds","cHmmReprOdds","bStatistic","minDistTSS", "minDistTSE")  
evoCol <- c("GerpN", "GerpS", "priPhyloP", "mamPhyloP", "verPhyloP", "priPhCons", "mamPhCons", "verPhCons")

# read snpList
setwd(inputDir)
snpList <- fread(inputFname, header=TRUE)

# read preprocess for replacing 0 (will be taking log)
setwd(paramDir)
load("min_nonzero_training.RData")
load("min_cHmm.RData")
load("mean_scores_training_epigenetic.RData")
load("sd_scores_training_epigenetic.RData")

# compile the CADD annotations
uniqueChrs <- unique(snpList$Chr)
annotData <- c()
notFound <- c()
for (chr_it in 1:length(uniqueChrs)) {

	tempChr <- uniqueChrs[chr_it]
	
	# cadd split data dir
	setwd(caddDir)
	tempCaddFile <- fread(paste0("cadd13split", tempChr, ".txt")) %>%
    set_colnames(c("Chrom", "Pos", "Ref", "Anc", "Alt", "Type", "Length", "isTv", "isDerived", "AnnoType", "Consequence", "ConsScore", "ConsDetail", "GC", "CpG", "mapAbility20bp", "mapAbility35bp", "scoreSegDup", "priPhCons", "mamPhCons", "verPhCons", "priPhyloP", "mamPhyloP", "verPhyloP", "GerpN", "GerpS", "GerpRS", "GerpRSpval", "bStatistic", "mutIndex", "dnaHelT", "dnaMGW", "dnaProT", "dnaRoll", "mirSVRScore", "mirSVRE", "mirSVRAln", "targetScan", "fitCons", "cHmmTssA", "cHmmTssAFlnk", "cHmmTxFlnk", "cHmmTx", "cHmmTxWk", "cHmmEnhG", "cHmmEnh", "cHmmZnfRpts", "cHmmHet", "cHmmTssBiv", "cHmmBivFlnk", "cHmmEnhBiv", "cHmmReprPC", "cHmmReprPCWk", "cHmmQuies", "EncExp", "EncH3K27Ac", "EncH3K4Me1", "EncH3K4Me3", "EncNucleo", "EncOCC", "EncOCCombPVal", "EncOCDNasePVal", "EncOCFairePVal", "EncOCpolIIPVal", "EncOCctcfPVal", "EncOCmycPVal", "EncOCDNaseSig", "EncOCFaireSig", "EncOCpolIISig", "EncOCctcfSig", "EncOCmycSig", "Segway", "tOverlapMotifs", "motifDist", "motifECount", "motifEName", "motifEHIPos", "motifEScoreChng", "TFBS", "TFBSPeaks", "TFBSPeaksMax", "isKnownVariant", "ESP_AF", "ESP_AFR", "ESP_EUR", "TG_AF", "TG_ASN", "TG_AMR", "TGAFR", "TG_EUR", "minDistTSS", "minDistTSE", "GeneID", "FeatureID", "CCDS", "GeneName", "cDNApos", "relcDNApos", "CDSpos", "relCDSpos", "protPos", "relProtPos", "Domain", "Dst2Splice", "Dst2SplType", "Exon", "Intron", "oAA", "nAA", "Grantham", "PolyPhenCat", "polyPhenVal", "SIFTcat", "SIFTval", "RawScore", "PHRED"))

	# look for simple match - just the BP
	tempChrList <- snpList %>% filter(Chr == tempChr) %>% distinct()
	bpMatch <- tempCaddFile %>% filter(Pos %in% tempChrList$BP) %>% distinct()

	# look for full match
	fullMatch <- merge(bpMatch, tempChrList, by.x=c("Pos", "Ref", "Alt"), by.y=c("BP", "Ref", "Alt")) %>%
		distinct()

	# merge full and partial matches 
	allMatch <- merge(bpMatch, fullMatch %>% mutate(fullM = 1), by=colnames(bpMatch), all=TRUE) %>%
		mutate(fullM = ifelse(is.na(fullM), 0, fullM)) %>% 		
		distinct()

	# which ones not found
	tempNotFound <- allMatch %>% select(Pos, Ref, Alt, fullM)  %>%
		merge(., tempChrList, by.x=c("Pos", "Ref", "Alt"), by.y=c("BP", "Ref", "Alt"), all=TRUE) %>%
		filter(is.na(fullM) | fullM == 0) %>%
		distinct()

	# note that the the cadd database sometimes has multiple entries for the 
	# same chr/pos/ref/alt (the difference might be consequence classified as "intergenic" vs. 
	# "regulatory" or something equally inane with the rest of the annotations being the same).
	
	# append	
	annotData <- rbind(annotData, allMatch)
	notFound <- rbind(notFound, tempNotFound)

	cat("done", tempChr, '\n')
}

# write output
setwd(outputDir)
write.table(annotData, paste0(outputFnameRoot, "_annot.txt"), append=F, quote=F, row.names=F, col.names=T, sep="\t")
write.table(notFound, paste0(outputFnameRooti "_missing.txt"), apiend=F, quote=F, row.names=F, col.names=T, sep="\t")

# read the annotation file
annotData <- fread(paste0(outputFnameRoot, "_annot.txt"), data.table=F)
noncodingInput <- annotData
noncodingInput <- noncodingInput %>% select(colUsed)

# collapse chromHmm
noncodingInput <- noncodingInput %>% 
	mutate(cHmmTSSclass = cHmmTssA + cHmmTssAFlnk + cHmmTssBiv + cHmmBivFlnk) %>%
	mutate(cHmmTxclass = cHmmTxFlnk + cHmmTx + cHmmTxWk) %>%
	mutate(cHmmEnhclass = cHmmEnhBiv + cHmmEnhG + cHmmEnh) %>%
	mutate(cHmmZnfclass = cHmmZnfRpts) %>%
	mutate(cHmmReprclass = cHmmHet + cHmmReprPC + cHmmReprPCWk) %>%
	mutate(cHmmQuiesclass = cHmmQuies)

# each row should sum to 1
sum(rowSums(noncodingInput[,5:19])) / nrow(noncodingInput)
sum(rowSums(noncodingInput[,41:46])) / nrow(noncodingInput)

# impute missing scores as 0 - fitCons, bStatistic, all EncXXX scores (didn't impute ChromHMM missing)
# don't need to fix GC, CpG, TFBSPeaksMax, minDistTSS, minDistTSE
noncodingInput[,imputeZero][is.na(noncodingInput[,imputeZero])] <- 0
sum(is.na(noncodingInput)) 

# impute the zero values to be min(non-zero value) / 2
for (index in c(1:46)){
  x <- noncodingInput[, index]
  noncodingInput[, index][x==0] <- min_nonzero_training[index]/2

	if (index > 40) {noncodingInput[, index][x==0] <- min_cHmm/2}
}

# calculate odds for the collapsed ChromHmm class
noncodingInput <- noncodingInput %>%
	mutate(cHmmTSSOdds = cHmmTSSclass / cHmmQuiesclass) %>%
	mutate(cHmmTxOdds = cHmmTxclass / cHmmQuiesclass) %>%
	mutate(cHmmEnhOdds = cHmmEnhclass / cHmmQuiesclass) %>%
	mutate(cHmmZnfOdds = cHmmZnfclass / cHmmQuiesclass) %>%
	mutate(cHmmReprOdds = cHmmReprclass / cHmmQuiesclass) %>%
	mutate(EncExp = log(EncExp)) %>%
	mutate(EncH3K27Ac = log(EncH3K27Ac)) %>%
	mutate(EncH3K4Me1 = log(EncH3K4Me1)) %>%
	mutate(EncH3K4Me3 = log(EncH3K4Me3)) %>%
	mutate(EncNucleo = log(EncNucleo)) %>%
	mutate(EncOCDNaseSig = log(EncOCDNaseSig)) %>%
	mutate(EncOCFaireSig = log(EncOCFaireSig)) %>%
	mutate(EncOCpolIISig = log(EncOCpolIISig)) %>%
	mutate(EncOCctcfSig = log(EncOCctcfSig)) %>%
	mutate(EncOCmycSig = log(EncOCmycSig)) %>%
	mutate(TFBS = log(TFBS)) %>%
	mutate(TFBSPeaks = log(TFBSPeaks)) %>%
	mutate(TFBSPeaksMax = log(TFBSPeaksMax)) %>%
	mutate(minDistTSS = -log(minDistTSS)) %>%
	mutate(minDistTSE = -log(minDistTSE)) %>%
	mutate(cHmmTSSOdds = log(cHmmTSSOdds)) %>%
	mutate(cHmmTxOdds = log(cHmmTxOdds)) %>%
	mutate(cHmmEnhOdds = log(cHmmEnhOdds)) %>%
	mutate(cHmmZnfOdds = log(cHmmZnfOdds)) %>% 
	mutate(cHmmReprOdds = log(cHmmReprOdds))

# extract columns for final use
noncodingInput <- noncodingInput[, finalEpi]
apply(noncodingInput, 2, summary)

# standardize all epigenetic scores to have mean 0 and variance 1
standNoncoding <- noncodingInput
for (index in 1:ncol(noncodingInput)){
  x <- noncodingInput[, index]
	trainingIdx <- which(names(mean_scores_training_epigenetic) == colnames(noncodingInput)[index])
  standNoncoding[, index] <- (x - mean_scores_training_epigenetic[trainingIdx])/sd_scores_training_epigenetic[trainingIdx]
}

# add variant ID and evolutionary class
standNoncoding <- cbind(annotData[, c(1,2,3,5)], annotData[, evoCol], standNoncoding)
head(standNoncoding)
dim(standNoncoding)

# get rid of missing conservation scores SNPs
standNoncoding <- na.omit(standNoncoding)
dim(standNoncoding)
sum(is.na(standNoncoding))

# save
setwd(outputDir)
write.table(standNoncoding,file=paste0(outputFnameRoot, "_ready.txt"), quote=FALSE,row.names=FALSE,sep="\t")















