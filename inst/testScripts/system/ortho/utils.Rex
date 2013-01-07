library("R.methodsS3");




setMethodS3("callByC1C2Grid", "PairedPSCBS", function(fit, grid, ...) {
  # Argument 'grid':
  if (!is.list(grid)) {
    throw("Argument 'grid' is not a list: ", class(grid)[1]);
  }
  if (length(grid) != 2) {
    throw("Argument 'grid' is not a list of length 2: ", length(grid));
  }
  grid <- lapply(grid, FUN=function(x) { Arguments$getDoubles(x); });

  data <- extractC1C2(fit);
  c1c2 <- as.matrix(data[,1:2,drop=FALSE]);
  rownames(c1c2) <- NULL;
  n <- nrow(c1c2);

  # Call C1 and C2 to the estimated modes
  names <- colnames(c1c2);
  callList <- list();
  for (kk in 1:2) {
    c <- c1c2[,kk,drop=TRUE];
    mus <- grid[[kk]];
    calls <- rep(as.integer(NA), times=n);
    minDists <- rep(Inf, times=n);
    for (cc in seq(along=mus)) {
      mu <- mus[cc];
      dist <- abs(c - mu);
      idxs <- which(dist < minDists);
      if (length(idxs) > 0) {
        minDists[idxs] <- dist[idxs];
        calls[idxs] <- cc;
      }
    } # for (cc ...)
    callList[[kk]] <- calls;
  } # for (kk ...)
  calls <- as.data.frame(callList);
  colnames(calls) <- sprintf("%s.call", names(pList));

  C1C2 <- cbind(muC1=grid[[1]][calls[,1]], muC2=grid[[2]][calls[,2]]);
  dC1C2 <- c1c2 - C1C2;
  colnames(dC1C2) <- sprintf("d%s", colnames(dC1C2));
  dist <- sqrt(dC1C2[,1]^2+dC1C2[,2]^2);
  calls <- cbind(calls, c1c2, C1C2, dC1C2, dist=dist);

  # Add weights (for now)
  n <- data[,4, drop=TRUE];
  n <- sqrt(n);
  w <- n/sum(n, na.rm=TRUE);

  calls <- cbind(calls, weights=w);

  calls;
}) # callByC1C2Grid()

setMethodS3("collapseCalls", "data.frame", function(calls, maxDist=Inf, ...) {
  callData <- subset(calls, dist <= maxDist);
  callData <- subset(callData, is.finite(C1.call));

  # Joint weights
  dist <- callData$dist;
  dw <- sqrt(1/(dist+0.001));
  dw <- 1;
  w <- dw * callData$weights;
  w <- w / sum(w, na.rm=TRUE);
  callData$w <- w;

  # Identify all calls
  callTag <- callData[,1:2];
  callTag <- sprintf("%02d,%02d", callTag[,1], callTag[,2]);
  callData <- cbind(callData, tag=callTag);

  tags <- sort(unique(callTag));

  uCallData <- NULL;
  keep <- c("C1.call", "C2.call", "muC1", "muC2");
  for (kk in seq(along=tags)) {
    callDataKK <- subset(callData, tag == tags[kk]);
    wKK <- sum(callDataKK$w, na.rm=TRUE);
    callDataKK <- cbind(callDataKK[1,keep], weights=wKK);
    uCallData <- rbind(uCallData, callDataKK);
  }
  rownames(uCallData) <- NULL;

  uCallData;  
}) 

setMethodS3("calcC2Shift", "data.frame", function(calls, maxDist=0.04, minWeight=0.05, ...) {
  callsT <- collapseCalls(calls, maxDist=maxDist);
  callsT <- subset(callsT, weights >= minWeight);
  d <- callsT[,"muC2"] - callsT[,"muC1"];
  min(d, na.rm=TRUE);
})

