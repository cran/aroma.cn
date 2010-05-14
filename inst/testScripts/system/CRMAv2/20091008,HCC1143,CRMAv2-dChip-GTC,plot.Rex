if (interactive()) savehistory();
library("aroma.affymetrix");
library("aroma.cn");
library("R.menu");

doSmooth <- TRUE;
rawCol <- "#000000";
if (doSmooth) {
  rawCol <- "#666666";
  rawCol <- "#999999";
  smoothTag <- "smooth";
} else {
  smoothTag <- "raw";
}

verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
fig <- 1;

methods <- c("CRMAv2", "GTC", "dChip");
trackCols <- c("blue", "red", "orange");
names(trackCols) <- methods;

dataSetFmt <- "GSE13372,HCC1143,onePair,%s,pairs";
dataSets <- sprintf(dataSetFmt, methods);
names(dataSets) <- methods;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load raw CN data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
dsList <- lapply(dataSets, FUN=function(dataSet) {
  AromaUnitTotalCnBinarySet$byName(dataSet, chipType="*");
});
# Keep only common samples (just in case)
names <- Reduce(intersect, lapply(dsList, FUN=getNames));
dsList <- lapply(dsList, FUN=extract, names);
dsList <- lapply(dsList, FUN=function(ds) {
  setFullNamesTranslator(ds, function(names, ...) {
    sprintf("%s,%s", getTags(ds)[1], names);
  });
  ds
});
print(dsList);


# Infer data set annotations
tags <- Reduce(intersect, lapply(dsList, FUN=getTags));
methods <- sapply(dsList, FUN=function(ds) setdiff(getTags(ds), tags));
platforms <- sapply(dsList, FUN=getPlatform);
chipTypes <- sapply(dsList, FUN=getChipType, fullname=FALSE);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Work with one array
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
array <- 1;
dfList <- lapply(dsList, getFile, array);
print(dfList);

ugpList <- lapply(dfList, getAromaUgpFile);
cdfList <- lapply(ugpList, getCdf);

verbose && enter(verbose, "Identifying all SNP units");
snpList <- lapply(cdfList, FUN=function(cdf) {
  which(getUnitTypes(cdf) == 2);
});
verbose && exit(verbose);

df <- dfList[[1]];
tags <- getTags(df);
label <- sprintf("%s (%s vs %s)", 
                getName(df), tags[1], gsub("ref=", "", tags[2]));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot all chromosomes
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
chromosomes <- getChromosomes(ugpList[[1]]);
if (interactive()) {
  chromosomes <- selectMenu(chromosomes);
}
nbrOfChromosomes <- length(chromosomes);


Clim <- c(0,6);
xScale <- 1e-6;

width <- 1280;
height <- 400;

# The name of the tumor/normal pair

figName <- sprintf("%s,tracks", label);
#devSet(figName);
#layout(seq(along=dfList));
#par(ask=(nbrOfChromosomes > 1));

for (cc in seq(length=nbrOfChromosomes)) {
  chr <- chromosomes[cc];
  chrTag <- sprintf("Chr%02d", chr);
  verbose && enter(verbose, sprintf("Chromosome #%d ('%s') of %d", 
                                     cc, chrTag, nbrOfChromosomes));

  verbose && enter(verbose, "Extracting signals");
  cnList <- list();
  for (kk in seq(along=dfList)) {
    df <- dfList[[kk]];
    snps <- snpList[[kk]];
    cn <- extractRawCopyNumbers(df, chromosome=chr, units=snps);
    cn$y <- 2*2^cn$y;
    cn$name <- gsub("log2", "", cn$name);
    cnList[[kk]] <- cn;
  } # for (kk ...)
  verbose && print(verbose, cnList);
  verbose && exit(verbose);

  xRange <- range(sapply(cnList, FUN=xRange));
  if (doSmooth) {
    verbose && enter(verbose, "Smoothing signals");
    cnSList <- lapply(cnList, FUN=function(cn) {
      binnedSmoothing(cn, from=xRange[1], to=xRange[2], by=200e3);
    });
    verbose && print(verbose, cnSList);
    verbose && exit(verbose);
  }

  xlim <- xRange;

  for (kk in seq(along=cnList)) {
    cn <- cnList[[kk]];
    if (doSmooth) {
      cnS <- cnSList[[kk]];
    }
    col <- trackCols[kk];

    fullname <- cn$name;
    method <- methods[kk];
    figPath <- filePath("reports", dataSets[kk], chipTypes[kk], smoothTag);
    figName <- sprintf("%s,%s,%s,tracks", fullname, chrTag, method);
    filename <- sprintf("%s.png", figName);
    pathname <- Arguments$getWritablePathname(filename, path=figPath);
    if (!isFile(pathname)) {
      devNew("png", pathname, width=width, height=height, label=figName);
      par(mar=c(3.2,4.3,2.5,2.1), mgp=c(2.8,2.0,0), cex.axis=3.5, cex.lab=4.0);
  
      plot(cn, col=rawCol, xlab="", ylab="", xlim=xlim*xScale, ylim=Clim, xScale=xScale);
      if (doSmooth) {
        points(cnS, cex=4, col="white", xScale=xScale);
        points(cnS, cex=1.4, col=col, xScale=xScale);
      }
      if (kk >= 1) {
        stext(side=3, pos=0, cex=2.5, label);
        stext(side=3, pos=1, cex=2.5, sprintf("%s", chrTag));
      }
      stext(side=4, line=-0.1, pos=0, cex=2.5, col="gray", method);
      stext(side=4, line=-0.1, pos=1, cex=2.5, col="gray", 
                                   sprintf("n=%d", nbrOfLoci(cn)));

      devDone();
    }
  } # for (kk ...)

  verbose && exit(verbose);
} # for (cc ...)
stop()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Extract subset of raw CNs across platforms
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Sample of interest
array <- 1;

# Chromosomal region of interest
chromosome <- 1;
region <- c(0,60)*1e6;

cnList <- lapply(dsList, FUN=function(ds) {
  df <- getFile(ds, array);
  cn <- extractRawCopyNumbers(df, chromosome=chromosome, region=region);
  cn$y <- 2^cn$y;
  cn;
});
#names(cnList) <- sprintf("%s\n%s\n%s", methods, platforms, chipTypes);
print(cnList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Smooth CNs using consecutive bins of width 100kb
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
xRange <- range(sapply(cnList, FUN=xRange));
cnSList <- lapply(cnList, FUN=function(cn) {
  binnedSmoothing(cn, from=xRange[1], to=xRange[2], by=100e3);
})
print(cnSList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot raw and binned CNs along chromosome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Clim <- c(0,6);
xlim <- xRange;

# The name of the tumor/normal pair
name <- getNames(dsList[[1]])[array];

devSet("panels");
layout(seq(along=cnList));
par(mar=c(4.2,4.2,1.3,2.1));
for (kk in seq(along=cnList)) {
  cn <- cnList[[kk]];
  cnS <- cnSList[[kk]];
  plot(cn, xlim=xlim, ylim=Clim);
  points(cnS, cex=1, col="white");
  points(cnS, cex=0.5, col="red");
  stext(side=3, pos=0, cex=0.8, name);
  stext(side=3, pos=1, cex=0.8, sprintf("Chr%02d (n=%d; s=%d)", 
                               chromosome, nbrOfLoci(cn), nbrOfLoci(cnS)));
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot pairwise binned CNs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
devSet("pairs");
Y <- sapply(cnSList, FUN=function(cs) getCNs(cs));
panel <- function(...) { abline(a=0, b=1, lty=2); points(...); }
pairs(Y, pch=20, col="red", lower.panel=NULL, upper.panel=panel, 
                                                     xlim=Mlim, ylim=Mlim);
