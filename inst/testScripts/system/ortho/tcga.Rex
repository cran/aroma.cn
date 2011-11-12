library("aroma.cn");
library("PSCBS");

verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

rootPath <- "pscbsData";
rootPath <- Arguments$getWritablePath(rootPath);

imgPath <- file.path("figures,2010-10-27", chipType);
mkdirs(imgPath);

Clim <- c(0,6);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup paired PSCN data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
chipType <- "GenomeWideSNP_6";
dataSet <- "broad.mit.edu_GBM.Genome_Wide_SNP_6.1.5.0,ASCRMAv2";
#dataSet <- "broad.mit.edu_OV.Genome_Wide_SNP_6.12.6.0,ASCRMAv2";

#chipType <- "HumanHap550";
#dataSet <- "hudsonalpha.org_GBM.HumanHap550.2.3.0,XY";

#chipType <- "Human1M-Duo";
#dataSet <- "hudsonalpha.org_OV.Human1MDuo.1.1.0,XY";


dsTuple <- TGCAPairedTotalAndFracBDataTuple$byName(dataSet, chipType=chipType);
print(dsTuple);
sampleNames <- getNames(dsTuple);
sampleNames <- setdiff(sampleNames, c("TCGA-02-0002"));

#sampleNames <- c("TCGA-04-1336", "TCGA-04-1335", "TCGA-04-1332", "TCGA-04-1331");
#sampleNames <- c("TCGA-02-0028");
sampleNames <- c("TCGA-02-0011");


# Extract samples of interest
dsTuple <- extract(dsTuple, sampleNames);
print(dsTuple);
sampleNames <- getNames(dsTuple);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Save PairedPSCBS data to binary files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (FALSE) {
  rootPath <- "wholeGenomeData";
  path <- file.path(rootPath, dataSet, chipType);
  path <- Arguments$getWritablePath(path);
  for (ii in seq(along=sampleNames)) {
    sampleName <- sampleNames[ii]; 
    verbose && enter(verbose, sprintf("Sample #%d ('%s') of %d", ii, sampleName, length(sampleNames)));
  
    dfList <- extractSample(dsTuple, sampleName);
    fullname <- attr(dfList, "fullname");
    filename <- sprintf("%s,paired,data.xdr", fullname);
    pathname <- file.path(path, filename);
    if (!isFile(pathname)) {
      verbose && enter(verbose, "Extracting and saving data");
      data <- extractPairedPSCBSData(dsTuple, sampleName=sampleName);
      verbose && str(verbose, data);
      saveObject(data, file=pathname);
      verbose && exit(verbose);
    } else {
      verbose && cat(verbose, "Already exported. Skipping.");
    }
  
    verbose && exit(verbose);
  } # for (ii ...)
} # if (...)


for (ii in seq(along=sampleNames)) {
  sampleName <- sampleNames[ii]; 

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Paired PSCBS
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit <- NULL;
  tryCatch({
    fit <- doPairedPSCBS(dsTuple, sampleName=sampleName, flavor="tcn,dh,tcn", verbose=verbose);
  }, error = function(ex) {});
  if (is.null(fit)) {
    next;
  }

  print(head(as.data.frame(fit)));
  print(tail(as.data.frame(fit)));
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Plot genomewide (C1,C2)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  imgName <- sprintf("%s,C1C2,raw", sampleName);
  filename <- sprintf("%s.png", imgName);
  pathname <- file.path(imgPath, filename);
  devNew("png", pathname, width=640, height=640);
  plotC1C2Grid(fit, Clim=Clim);
  linesC1C2(fit);
  stext(side=3, pos=0, sampleName);
  title(main="(C1,C2) - raw");
  devDone();
} # for (ii ...)

stop()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calibrate
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
verbose && enter(verbose, "Calibrate in (C1,C2) space");
fitO <- calibrateC1C2(fit, verbose=verbose, debug=TRUE);
verbose && exit(verbose);

# Plot genomewide (C1,C2)
Clim <- c(-1,10);
devSet("C1C2,calibrated");
plotC1C2Grid(fitO, Clim=Clim);
linesC1C2(fitO);
stext(side=3, pos=0, sampleName);
title(main="(C1,C2) - calibrated");



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot tracks chromosome by chromosome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
devSet("C1C2,tracks,byChromosome");
par(ask=TRUE);
for (chr in 1:22) {
  chrTag <- sprintf("chr%02d", chr);

  fitListT <- list(
    extractByChromosome(fit, chr),
    extractByChromosome(fitO, chr)
  );

  subplots(2, ncol=1, byrow=FALSE);
  par(mar=c(1,3.5,1,0.5)+1);
  for (fit in fitListT) {
    plotTracks(fit, tracks="tcn,c1,c2", Clim=Clim, add=TRUE);
    stext(side=3, pos=0, sampleName);
    stext(side=1, pos=1, line=-1, chrTag);
  }
} # for (chr ...)
