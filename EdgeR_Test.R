# ---- This exercise is from https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/rnaSeq_DE.pdf ---- #

# ---- Import data ---- #
library(edgeR)
load("~/Desktop/Counts.RData")
Counts <- tmp$counts
colnames(Counts) <- c("16N", "16T", "18N", "18T", "19N", "19T")
# ---- inspecting the data  ---- #

dim(Counts) # outputs dimensions

head(Counts) # outputs the first few lines of Counts

#---------------------------------------------------------------------------
# Identify differentially expressed genes using 'edgeR' specific commands
#---------------------------------------------------------------------------

# Creating a DGEList object
dgList <- DGEList(counts=Counts, genes=rownames(Counts)) #Creates dgList object
dgList # outputs some of dgList
dgList$samples # outputs the first 7 lines of the file???
head(dgList$counts) #see above
head(dgList$genes) #see above

## ----- Filtering the dataset for genes not expressed in both samples ---- #
countsPerMillion <- cpm(dgList)
summary(countsPerMillion) #'summary' is a useful function for exploring numeric data
countCheck <- countsPerMillion > 1 # creates an object 'countCheck' that has TRUE/FALSE data for what genes have more than 0 counts in each sample
head(countCheck)

keep <- which(rowSums(countCheck) >= 2) # I think this gets rid of those items in countCheck that have less than 2 samples in which they are expressed
dgList <- dgList[keep,]
summary(cpm(dgList))

# ---- Normalizing the data ---- #
# R uses the trimmed mean of M-values (TMM) method
?calcNormFactors
dgList <- calcNormFactors(dgList, method="TMM") # overwrites the dgList object to have normalized values instead of count data

# ---- Data Exploration
# Multidimensional scaling plot of dgList
plotMDS(dgList) # This plot should more or less separate the samples based on treatment and sample number

# ---- Building our model for statistics ---- #
# This matrix tells edgeR how the experiment was setup and how that setup is coded in the sample names
sampleType <- rep("N", ncol(dgList)) #N=normal; T=tumor
sampleType[grep("T", colnames(dgList))] <- "T" #'grep' is a string-matching function

sampleReplicate <- paste("S", rep(1:3, each=2), sep="")

designMat <- model.matrix(~sampleReplicate + sampleType)
designMat

# ---- Estimating Dispersions ---- #
# of the types of dispersion methods (common dispersion, trended dispersion, tagwise) we will use tagwise
# see .pdf for a longer definition of what we are doing here

dgList <- estimateGLMCommonDisp(dgList, design=designMat) #must run commondisp before tagwise
dgList <- estimateGLMTagwiseDisp(dgList, design=designMat)

plotBCV(dgList) # Biological coefficient of variation (BCV) is the sqare root of the dispersion parameter int he negative binomial model

# ---- Differential Expression calling ---- #
# after fitting the model, we can use functions (topTags) to explore the results and to define our thresholds for identifying DE genes

fit <- glmFit(dgList, designMat)
lrt <- glmLRT(fit, coef=4)

?topTags
edgeR_result <- topTags(lrt)

edgeRtable <- topTags(lrt,n=15000)$table
save(edgeRtable, file = 'edgeR_Results.RData')

# now we can plot the log-fold changes of all genes and highlight those that are differentially expressed
?decideTests
deGenes <- decideTestsDGE(lrt, p=0.001)
deGenes <- rownames(lrt)[as.logical(deGenes)]
plotSmear(lrt, de.tags = deGenes)
abline(h=c(-2,2), col=2)