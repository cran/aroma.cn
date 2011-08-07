if (interactive()) savehistory();
library("aroma.cn");
library("psCBS");

log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

dataSet <- "psCBS,exampleData";
chipType <- "HumanHap550";
chipType <- "GenomeWideSNP_6";

dsList <- list();
dsList$total <- AromaUnitTotalCnBinarySet$byName(dataSet, chipType=chipType);
dsList$fracB <- AromaUnitFracBCnBinarySet$byName(dataSet, chipType=chipType);

filename <- ".average-signals-median-mad,16a9fce528d50d41b0144cc872558456.asb";
dfR <- AromaUnitTotalCnBinaryFile(filename, path=getPath(dsList$total));

ugp <- getAromaUgpFile(dsList[["total"]]);
units <- getUnitsOnChromosome(ugp, chromosome=1);
pos <- getPositions(ugp, units=units);

theta <- 2^extractMatrix(dsList$total, units=units);
beta <- extractMatrix(dsList$fracB, units=units);

CT <- 2*theta[,1]/theta[,2];
betaT <- beta[,1];
betaN <- beta[,2];

# SNPs only
isSNP <- is.finite(betaT);

cn <- RawCopyNumbers(cn=CT[isSNP], x=pos[isSNP], chromosome=1);
plot(cn, pch=".", ylim=c(0,5));

fit <- segmentByCBS(cn);
cnr <- extractCopyNumberRegions(fit);
drawLevels(cnr, lwd=2, col="red");
