library("aroma.cn"); 

verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup normal and tumor profiles
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The normal/germline CN profile (TCN=2 everywhere)
pN <- cnr(1,1000, 2);
print(pN);

# Build up tumor CN profile generation by generation
pT1 <- pN + cnr(250,450);
print(pT1);

pT2 <- pN - cnr(550,750);
print(pT2);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulate copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
set.seed(0xBEEF);
cnT1 <- simulateRawCopyNumbers(pT1, n=2000, sd=1/3);
cnT2 <- simulateRawCopyNumbers(pT2, n=2000, sd=1/3);
profiles <- list(
  "gain"=list(cnr=pT1, cn=cnT1, col="blue", lty=3),
  "loss"=list(cnr=pT2, cn=cnT2, col="red", lty=3)
);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot profiles with data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subplots(length(profiles), ncol=1);
par(mar=c(0.1,4,2,2)+0.1);
for (kk in seq(along=profiles)) {
  p <- profiles[[kk]];
  xlim <- c(-50, xMax(p$cn)/1e6);
  plot(p$cn, col="gray", xlim=xlim, ylim=c(0,5), axes=FALSE);
  axis(side=2);
  lines(p$cnr, col="white", lty=1, lwd=8);
  lines(p$cnr, col=p$col, lty=1, lwd=5);
  drawDensity(p$cnr, col=p$col, height=50, offset=5);
  title(main=names(profiles)[kk]);
} # for (kk ...)
