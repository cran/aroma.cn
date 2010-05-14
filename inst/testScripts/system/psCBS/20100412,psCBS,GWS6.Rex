if (interactive()) savehistory();
library("aroma.cn");
library("psCBS");

log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

dataSet <- "psCBS,exampleData,paired";
chipType <- "HumanHap550";
chipType <- "GenomeWideSNP_6";

dsList <- list();
dsList$total <- AromaUnitTotalCnBinarySet$byName(dataSet, chipType=chipType);
dsList$fracB <- AromaUnitFracBCnBinarySet$byName(dataSet, chipType=chipType);

ugp <- getAromaUgpFile(dsList[["total"]]);
units <- getUnitsOnChromosome(ugp, chromosome=1);
pos <- getPositions(ugp, units=units);

theta <- extractMatrix(dsList$total, units=units);
beta <- extractMatrix(dsList$fracB, units=units);

CT <- 2*theta[,1]/theta[,2];
betaT <- beta[,1];
betaN <- beta[,2];

# SNPs only
isSNP <- is.finite(betaT);

subplots(2, ncol=1);
par(mar=c(3,5,1,1)+0.1);

cn <- RawCopyNumbers(cn=CT[isSNP], x=pos[isSNP], chromosome=1);
plot(cn, pch=".", ylim=c(0,5));
stext(side=3, pos=1, sprintf("Chr %d", cn$chromosome), cex=0.7);
fit <- segmentByCBS(cn);
cnr <- extractCopyNumberRegions(fit);
drawLevels(cnr, lwd=2, col="red");


cn <- RawCopyNumbers(cn=sqrt(CT[isSNP]), x=pos[isSNP], chromosome=1);
plot(cn, pch=".", ylim=c(0,5), ylab=substitute(sqrt(C)));
stext(side=3, pos=1, sprintf("Chr %d", cn$chromosome), cex=0.7);
fit <- segmentByCBS(cn);
cnr <- extractCopyNumberRegions(fit);
drawLevels(cnr, lwd=2, col="red");
