## PM140316 - R3.2.2
##
## Analysis of predicted subtypes for rebuttal BMC Medical Genomics

library(purrr)
library(RColorBrewer)

load("Evaluation.RData")
# PredictorAnnotation: data frame that describes each of the 489 subtype predictors
# SA: character matrix describing the predicted subtypes for all samples and subtype predictors
# SampleAnnotation: data frame that describes each of the 5386 samples
all(PredictorAnnotation$Name == colnames(SA))
all(SampleAnnotation$SampleID == rownames(SA))

levels(as.factor(PredictorAnnotation$Train))
# Select the five datasets used as consensus sets
Training <- c("expO","Guedj","Sabatier","Bos","Li")

table(SampleAnnotation$Study)
# D1-D22: VDX = Wang + Yu
Datasets <- c("Bos","Dedeurwaerder","Desmedt","EXPO","Farmer","Guedj","Kao","Li","Lu","Miller","MSK","Pawitan","Richardson.1","Richardson.2","Sabatier","Schmidt","Shi","Symmans.1","Symmans.2","Symmans.3","UNT","Wang","Yu")
length(Datasets)

# We 'just' need to extract the right subsets from SA
SCM.cs <- grep("SCM\\.(Bos|expO|Guedj|Li|Sabatier)\\.(D|HK|W)\\.cs5$",colnames(SA),value=TRUE)
# What does 'mk' stand for?
SSP.cs <- grep("SSP\\.(Bos|expO|Guedj|Li|Sabatier)\\.(H|P|S)\\.mk\\.cs5$",colnames(SA),value=TRUE)
SSP.4s <- c("pam50.robust.4s","ssp2003.robust.4s","ssp2006.robust.4s")
SSP.5s <- c("pam50.robust","ssp2003.robust","ssp2006.robust")
SCM    <- c("scmgene.robust","scmod1.robust","scmod2.robust")

# Read in information about the rejected samples from QC.txt (Table S7, Additional File 1)
dat <- read.table("QC.txt",sep="&",header=TRUE,stringsAsFactors=FALSE,allowEscapes=TRUE)
indx <- match(gsub(" ","",as.character(dat$CEL)),rownames(SA))
indx 

# Leave out QC rejected samples
SA <- SA[-indx,]
SampleAnnotation <- SampleAnnotation[-indx,]

table(SampleAnnotation$Study)

calcFrequencies <- function(subtype.pred,subtypes=c("Basal","Her2","LumA","LumB")){
 # Use 'transpose' to create a nested list of subtype assignments per dataset, per predictor type
 ll <- transpose(apply(SA[,subtype.pred],2,split,as.factor(SampleAnnotation$Study)))
 # Average contingency tables per dataset
 mat <- sapply(ll,function(x){
   f <- factor((unlist(x)))
   levels(f) <- c(levels(f),setdiff(subtypes,levels(f)))
   table(factor(f,levels=subtypes))/length(subtype.pred)
 })
 mat <- mat[,Datasets]
 round(100*sweep(mat,2,colSums(mat),`/`),2)
}

subtypes.5s <- c("Basal","Her2","LumA","LumB","Normal")

freq.SCM.cs <- calcFrequencies(SCM.cs)
freq.SSP.cs <- calcFrequencies(SSP.cs)
freq.SCM    <- calcFrequencies(SCM)
freq.SSP.4s <- calcFrequencies(SSP.4s)
freq.SSP.5s <- calcFrequencies(SSP.5s,subtypes.5s)

createBarplot <- function(freq,subtype.pred,subtypes=c("Basal","Her2","LumA","LumB")){
  mycols <- brewer.pal(length(subtypes),"Blues")
  par(mar=c(5,8,4,2)+0.1)
  barplot(freq,horiz=TRUE,beside=TRUE,xlab="percentage",las=1,col=mycols,space=c(0,1.5),main=subtype.pred)
  legend("bottomright",legend=subtypes,fill=mycols)
}
pdf("subtypeFrequencies.pdf")
createBarplot(freq.SCM.cs,"SCM.cs")
createBarplot(freq.SSP.cs,"SSP.cs")
createBarplot(freq.SCM,"SCM")
createBarplot(freq.SSP.4s,"SSP.4s")
createBarplot(freq.SSP.5s,"SSP.5s",subtypes.5s)
dev.off()

# Calculate subtype percentages per predictor type in order to compare with
# Haibe-Kains, JNCI (2012), p.7
calcFrequenciesPredictor <- function(subtype.pred,subtypes=c("Basal","Her2","LumA","LumB")){
  SA.subset <- subset(SA,SampleAnnotation$Study %in% Datasets)
  SampleAnnotation.subset <- subset(SampleAnnotation,SampleAnnotation$Study %in% Datasets)
  # A nested list of subtype assignments per predictor type, per dataset
  ll <- apply(SA.subset[,subtype.pred],2,split,as.factor(SampleAnnotation.subset$Study))
  # Average contingency tables per predictor stype
  mat <- sapply(ll,function(x){
    f <- factor((unlist(x)))
    levels(f) <- c(levels(f),setdiff(subtypes,levels(f)))
    table(factor(f,levels=subtypes))/length(subtype.pred)
  })
  round(100*sweep(mat,2,colSums(mat),`/`),2)
}
apply(calcFrequenciesPredictor(SCM.cs),1,range)
apply(calcFrequenciesPredictor(SSP.cs),1,range)
apply(calcFrequenciesPredictor(SCM),1,range)
apply(calcFrequenciesPredictor(SSP.4s),1,range)
apply(calcFrequenciesPredictor(SSP.5s,subtypes.5s),1,range)
