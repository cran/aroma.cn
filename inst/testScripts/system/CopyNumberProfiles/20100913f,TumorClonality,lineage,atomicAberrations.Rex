library("aroma.cn"); 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
xAxisBar <- function(x0, x1=NULL, y=par("usr")[3], lwd=4, col=par("col"), length=0.04, xpd=FALSE, ...) {
  if (is.null(x1)) {
    x1 <- x0[2];
    x0 <- x0[1];
  }
  arrows(x0=x0, x1=x1, y0=y, code=3, angle=90, length=length, lwd=lwd, col=col, xpd=xpd, lend=2, ...);
} # xAxisBar()

# Get apply
gapply <- function(X, ..., FUN="[[", useNames=FALSE) {
  res <- sapply(X=X, FUN=FUN, ...);
  if (!useNames) {
    names(res) <- NULL;
  }
  res;
} # gapply()


verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup normal and tumor profiles
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The normal/germline CN profile (TCN=2 everywhere)
pN <- cnr(1,1000, 2);
print(pN);

flavor <- "C**";
if (flavor == "A") {
  # Tumor clone #1
  pT1 <- pN +
         cnr(400,500) + 
         cnr(600,800);
  print(pT1);
  
  # Tumor clone #2
  pT2 <- pT1 +
         cnr(600,700) + 
         cnr(100,200) - cnr(850,900);
  print(pT2);
} else if (flavor == "A*") {
  # Tumor clone #1
  pT1 <- pN +
         cnr(400,500);
  print(pT1);
  
  # Tumor clone #2
  pT2 <- pT1 +
         cnr(100,200) - cnr(850,900);
  print(pT2);
} else if (flavor == "A**") {
  # Tumor clone #1
  pT1 <- pN +
         cnr(400,500) + 
         cnr(600,800);
  print(pT1);
  
  # Tumor clone #2
  pT2 <- pT1 +
         cnr(100,200) - cnr(850,900);
  print(pT2);
} else if (flavor == "B") {
  # Tumor clone #1
  pT1 <- pN +
         cnr(600,700) + 
         cnr(100,200) - cnr(850,900);
  print(pT1);
  
  # Tumor clone #2
  pT2 <- pT1 +
         cnr(400,500) + 
         cnr(600,800);
  print(pT2);
} else if (flavor == "C*") {
  # Tumor clone #1
  pT1 <- pN +
         cnr(400,500) + 
         cnr(300,600);
  print(pT1);
  
  # Tumor clone #2
  pT2 <- pT1 +
         cnr(100,200) - cnr(850,900);
  print(pT2);
} else if (flavor == "C**") {
  # Tumor clone #1
  pT1 <- pN +
         cnr(300,600);
  print(pT1);
  
  # Tumor clone #2
  pT2 <- pT1;
  print(pT2);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulate copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cnN <- simulateRawCopyNumbers(pN, n=2000, sd=1/3);
cnT1 <- simulateRawCopyNumbers(pT1, n=2000, sd=1/3);
cnT2 <- simulateRawCopyNumbers(pT2, n=2000, sd=1/3);


profiles <- list(
  "Tumor at time t=1"=list(cnr=pT1, cn=cnT1, col="blue"),
  "Tumor at time t=2"=list(cnr=pT2, cn=cnT2, col="magenta")
);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Different tumor clonalities
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
kappas <- c(0.5, 0.7);
ltys <- c(1,1,1);
for (kk in seq(along=kappas)) {
  kappa <- kappas[kk];
  key <- sprintf("%d%% Tumor t=1 + %d%% Tumor t=2", 
                   round(100*kappa), round(100*(1-kappa)));
  cnr <- kappa *  pT1 + (1-kappa)* pT2;
  cn  <- kappa * cnT1 + (1-kappa)*cnT2;
  profiles[[key]] <- list(cnr=cnr, cn=cn, col="red");
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
H <- 1;
profiles <- profiles[length(profiles)];
subplots(length(profiles), ncol=1);
par(mar=c(0.1,4,2,2)+0.1);
for (kk in seq(along=profiles)) {
  p <- profiles[[kk]];
  xlim <- c(-50, xMax(p$cn)/1e6);
  plot(p$cn, col="gray", xlim=xlim, ylim=c(0,4.5), axes=FALSE);
  axis(side=2);
  lines(p$cnr, col="white", lty=1, lwd=8);
  lines(p$cnr, col="black", lty=1, lwd=5);
  drawDensity(p$cnr, col="black", height=50, offset=5);
  stext(side=3, pos=0.5, names(profiles)[kk], cex=1);

  fit <- findAtomicAberrations(cnr=p$cnr, data=p$cn, H=H, verbose=verbose);
  df <- fit$res;

  for (kk in seq(length=nrow(df))) {
    dfKK <- df[kk,];
    segments <- as.integer(dfKK[,c("firstRegion", "lastRegion")]);
    segments <- segments[1]:segments[2];
    xRange <- as.double(dfKK[,c("start", "stop")]);
    cnrKK <- subset(p$cnr, subset=segments);
    drawLevels(cnrKK, col="white", lwd=6);
    drawLevels(cnrKK, col="red", lwd=5);
    x <- xRange/1e6;
    y <- y0 + 0.2*dy;
#    xAxisBar(x0=x[1], x1=x[2], y=y, col=col);
  } 
} # for (kk ...)
