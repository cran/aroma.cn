setConstructorS3("TGCAPairedTotalAndFracBDataTuple", function(dsTotalT=NULL, dsTotalN=NULL, dsBetaT=NULL, dsBetaN=NULL, ...) {
  if (is.list(dsTotalT)) {
    dsList <- dsTotalT;
    dsTotalT <- dsList$total$T;
    dsTotalN <- dsList$total$N;
    dsBetaT <- dsList$fracB$T;
    dsBetaN <- dsList$fracB$N;
  }

  extend(Object(), "TGCAPairedTotalAndFracBDataTuple",
    dsTotalT = dsTotalT,
    dsTotalN = dsTotalN,
    dsBetaT = dsBetaT,
    dsBetaN = dsBetaN
  );
})

setMethodS3("getNames", "TGCAPairedTotalAndFracBDataTuple", function(this, ...) {
  dsList <- as.list(this);
  names <- NULL;
  for (kk in seq(along=dsList)) {
    ds <- dsList[[kk]];
    namesKK <- getNames(ds, ...);
    if (kk == 1) {
      names <- namesKK;
    } else {
      names <- intersect(names, namesKK);
    }
  } # for (kk ...)
  names;
})

setMethodS3("as.character", "TGCAPairedTotalAndFracBDataTuple", function(x, ...) {
  # To please R CMD check
  this <- x;

   s <- sprintf("%s:", class(this)[1])
   s <- c(s, sprintf("Names: %s", paste(getNames(this), collapse=", ")));
   class(s) <- "GenericSummary"
   s;
})


setMethodS3("getFullNamesTranslator", "TGCAPairedTotalAndFracBDataTuple", function(static, ...) {
  # Fullnames translator
  fnt <- function(names, ...) {
    pattern <- "^(TCGA-[0-9]{2}-[0-9]{4})-([0-9]{2}[A-Z])[-]*(.*)";
    gsub(pattern, "\\1,\\2,\\3", names);
  } # fnt()

  fnt;
}, static=TRUE, protected=TRUE);


setMethodS3("byName", "TGCAPairedTotalAndFracBDataTuple", function(static, dataSet, tags=NULL, chipType, ..., paths="totalAndFracBData") {
  dataSet <- Arguments$getCharacter(dataSet);
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
  }
  chipType <- Arguments$getCharacter(chipType);

  paths <- sapply(paths, FUN=function(path) {
    path <- Arguments$getReadablePath(path, mustExist=TRUE);
  });

  rootPath <- paths[1];

  fullname <- paste(c(dataSet, tags), collapse=",");
  path <- file.path(rootPath, fullname, chipType);
  path <- Arguments$getReadablePath(path, mustExist=TRUE);

  # Full names translator
  fnt <- TGCAPairedTotalAndFracBDataTuple$getFullNamesTranslator();

  # Load the total copy numbers
  ds <- AromaUnitTotalCnBinarySet$byName(dataSet, tags=tags, chipType=chipType);
  setFullNamesTranslator(ds, fnt);

  types <- sapply(ds, FUN=function(df) getTags(df)[1]);
  types <- gsub("([0-9]+)[A-Z]", "\\1", types);
  isT <- is.element(types, c("01", "02"));
  isN <- is.element(types, c("10", "11"));
  dsTotalT <- extract(ds, isT);
  dsTotalN <- extract(ds, isN);

  # Load the allele B fractions
  ds <- AromaUnitFracBCnBinarySet$byName(dataSet, tags=tags, chipType=chipType);
  setFullNamesTranslator(ds, fnt);

  types <- sapply(ds, FUN=function(df) getTags(df)[1]);
  types <- gsub("([0-9]+)[A-Z]", "\\1", types);
  isT <- is.element(types, c("01", "02"));
  isN <- is.element(types, c("10", "11"));
  dsBetaT <- extract(ds, isT);
  dsBetaN <- extract(ds, isN);

  res <- list(
    total = list(T=dsTotalT, N=dsTotalN),
    fracB = list(T=dsBetaT, N=dsBetaN)
  );

  res <- TGCAPairedTotalAndFracBDataTuple(res);

  res;
}, static=TRUE)


setMethodS3("as.list", "TGCAPairedTotalAndFracBDataTuple", function(this, ...) {
  res <- list(
    totalT = this$dsTotalT, 
    totalN = this$dsTotalN,
    fracBT = this$dsBetaT, 
    fracBN = this$dsBetaN
  );

  res;
})


setMethodS3("extract", "TGCAPairedTotalAndFracBDataTuple", function(this, ...) {
  res <- clone(this);
  
  fields <- c("dsTotalT", "dsTotalN", "dsBetaT", "dsBetaN");
  for (field in fields) {
    ds <- this[[field]];
    ds <- extract(ds, ...);
    res[[field]] <- ds;
  }

  res;
})

setMethodS3("extractSample", "TGCAPairedTotalAndFracBDataTuple", function(this, sampleName, ...) {
  sampleName <- Arguments$getCharacter(sampleName);

  dsList <- as.list(this);
  str(dsList);
  dsList <- lapply(dsList, FUN=function(ds) {
    cat("Data set: ", getFullName(ds), "\n", sep="");
    cat("Sample names:\n");
    print(getNames(ds));
    cat("Sample name: ", sampleName, "\n", sep="");
    idxs <- indexOf(ds, sampleName);
    getFile(ds, idxs);
  });

  name <- getName(dsList[[1]]);
  tagsT <- getTags(dsList[[1]]);
  tagsT <- setdiff(tagsT, "total");
  tagsN <- getTags(dsList[[2]]);
  tagsN <- setdiff(tagsN, "total");
  tags <- paste(tagsT, tagsN, sep="vs");
  fullname <- paste(c(name, tags), collapse=",");

  attr(dsList, "name") <- name;
  attr(dsList, "tags") <- tags;
  attr(dsList, "fullname") <- fullname;

  dsList;
}, protected=TRUE)


setMethodS3("extractPairedPSCBSData", "TGCAPairedTotalAndFracBDataTuple", function(this, ..., na.rm=TRUE) {
  dsList <- extractSample(this, ...);
  extractPairedPSCBSData(dsList, na.rm=na.rm);
})



setMethodS3("extractPairedPSCBSData", "list", function(dfList, na.rm=TRUE, ...) {
  dfTT <- dfList[[1]];
  dfTN <- dfList[[2]];
  dfBT <- dfList[[3]];
  dfBN <- dfList[[4]];

  # Sanity checks
  dfTT <- Arguments$getInstanceOf(dfTT, "AromaUnitTotalCnBinaryFile");
  dfTN <- Arguments$getInstanceOf(dfTN, "AromaUnitTotalCnBinaryFile");
  dfBT <- Arguments$getInstanceOf(dfBT, "AromaUnitFracBCnBinaryFile");
  dfBN <- Arguments$getInstanceOf(dfBN, "AromaUnitFracBCnBinaryFile");

  # Get the UGP annotation data
  ugp <- getAromaUgpFile(dfTT);
  gp <- readDataFrame(ugp);
  
  chromosome <- gp$chromosome;
  position <- gp$position;
  thetaT <- dfTT[,1,drop=TRUE];
  thetaN <- dfTN[,1,drop=TRUE];
  betaT <- dfBT[,1,drop=TRUE];
  betaN <- dfBN[,1,drop=TRUE];

  # Not needed anymore
  rm(ugp, gp);

  data <- data.frame(
    chromosome = chromosome,
    position = position,
    thetaT = thetaT,
    thetaN = thetaN,
    betaT = betaT,
    betaN = betaN
  );

  # Drop missing values?
  if (na.rm) {
    keep <- which(is.finite(chromosome) & is.finite(position));
    data <- data[keep,,drop=FALSE];

    keep <- which(is.finite(thetaT) & is.finite(thetaN));
    data <- data[keep,,drop=FALSE];

    # Not needed anymore
    rm(keep);
  }

  data$CT <- 2 * data$thetaT/data$thetaN;

  res <- list(
    name = attr(dfList, "name"),
    tags = attr(dfList, "tags"),
    fullname = attr(dfList, "fullname"),
    data = data
  );

  res;
})


setMethodS3("doPairedPSCBS", "TGCAPairedTotalAndFracBDataTuple", function(this, sampleName, ..., flavor=c("tcn,dh", "tcn,dh,tcn"), rootPath="pairedPSCBSData/", force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'sampleName':
  sampleName <- Arguments$getCharacter(sampleName);

  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'rootPath':
  rootPath <- Arguments$getWritablePath(rootPath);

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 
  

  verbose && enter(verbose, "PairedPSCBS");
  verbose && cat(verbose, "Sample name: ", sampleName);
  verbose && cat(verbose, "Flavor: ", flavor);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ds <- this$dsTotalT;
  dataSet <- getFullName(ds);
  chipTypeS <- getChipType(ds, fullname=FALSE);
  path <- file.path(rootPath, dataSet, chipTypeS);
  path <- Arguments$getWritablePath(path);
  verbose && cat(verbose, "Output path: ", path);


  tags <- c("pairedPSCBS");
  if (flavor == "tcn,dh,tcn") {
    tags <- c(tags, "psTCN");
  }
  fullname <- paste(c(sampleName, tags), collapse=",");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already done?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  filename <- sprintf("%s.xdr", fullname);
  pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=FALSE);
  if (!force && isFile(pathname)) {
    verbose && cat(verbose, "Loading already processed results.");
    verbose && cat(verbose, "Pathname: ", pathname);
    res <- loadObject(pathname);
    verbose && exit(verbose);
    return(res);
  }


  # Load "tcn,dh" data?
  res <- NULL;
  if (flavor == "tcn,dh,tcn") {
    fullname <- paste(c(sampleName, tags[-length(tags)]), collapse=",");
    filename <- sprintf("%s.xdr", fullname);
    pathnameT <- Arguments$getReadablePathname(filename, path=path, mustExist=FALSE);
    if (!force && isFile(pathnameT)) {
      verbose && cat(verbose, "Loading some of the results.");
      verbose && cat(verbose, "Pathname: ", pathnameT);
      res <- loadObject(pathnameT);
      verbose && exit(verbose);
    }
  }

  if (is.null(res)) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Load data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Loading data");
    verbose && cat(verbose, "Sample name: ", sampleName);
    data <- extractPairedPSCBSData(dsTuple, sampleName, na.rm=TRUE);
    verbose && str(verbose, data);
    verbose && exit(verbose);

    attachLocally(data);

    # AD HOC: Robustification
    CT[CT < 0] <- 0;
    CT[CT > 30] <- 30;

    # Chromosome to segment
    chromosomes <- sort(unique(chromosome));
    verbose && cat(verbose, "Chromosomes: ", seqToHumanReadable(chromosomes));
    chromosomes <- setdiff(chromosomes, 23:25);
    verbose && cat(verbose, "Filtered chromosomes: ", seqToHumanReadable(chromosomes));

    verbose && enter(verbose, "Segmenting data");
    res <- segmentByPairedPSCBS(CT, betaT=betaT, betaN=betaN, chromosome=chromosome, x=position, verbose=verbose);

    verbose && print(verbose, head(as.data.frame(res)));
    verbose && print(verbose, tail(as.data.frame(res)));

    verbose && exit(verbose);
  } # if (is.null(res))


  # Postsegment TCN?
  if (flavor == "tcn,dh,tcn") {
    verbose && enter(verbose, "Post-segmenting TCNs");
    res <- postsegmentTCN(res, verbose=verbose);
    verbose && str(verbose, res);
    verbose && print(verbose, head(as.data.frame(res)));
    verbose && exit(verbose);
  }

  verbose && enter(verbose, "Saving PairedPSCBS results");
  verbose && cat(verbose, "Pathname: ", pathname);
  saveObject(res, file=pathname);
  verbose && exit(verbose);

  res;
}) # doPairedPSCBS()



############################################################################
# HISTORY:
# 2010-10-18
# o Added argument 'na.rm' to extractPairedPSCBSData().
# o Updated doPairedPSCBS() to utilize the new multi-chromosome support
#   of segmentByPairedPSCBS().
# 2010-10-14
# o Added doPairedPSCBS().
# 2010-10-13
# o Created.
############################################################################  
