### Code for computing Density Recovery Profile (DRP).
### Wed 13 Feb 2002.


autodrp <- function(xs, ys, nbins, r, a=NULL) {
  ## Compute the autoDRP.  This is a simple wrapper around the cross DRP.
  crossdrp(xs, ys, xs, ys, nbins, r, a)
}

crossdrp <- function(xs1, ys1, xs2, ys2, nbins, r, a=NULL) {
  ## Compute the crossDRP.
  if (is.null(a)) {
    ## If a was not specified, we calculate bounds of dataset from the
    ## array.
    lims <- range(c(xs1, xs2));
    left <- lims[1]; right <-lims[2];

    lims <- range(c(ys1, ys2));
    bottom <- lims[1]; top <-lims[2];
    ##cat(paste(left, right, bottom, top, "\n"))
  }
  else {
    if (is.numeric(a)) {
      ## From the a vector, extract bounds of the area to study.
      left <- a[1]; right <- a[2]; bottom <- a[3]; top <- a[4];
    }
    else {
      stop("a is invalid - should be a vector of 4 numbers\n")
    }
  }
  
  l <- top-bottom; w <- right-left;

  ## now filter out points that are outside the bounding area A.
  subset1.x <- ((xs1 >= left) & (xs1 <= right))
  subset1.y <- ((ys1 >= bottom) & (ys1 <= top))
  subset1 <- which( subset1.x & subset1.y)
  xs1 <- xs1[subset1]; ys1 <- ys1[subset1]
  
  subset2.x <- ((xs2 > left) & (xs2 < right))
  subset2.y <- ((ys2 > bottom) & (ys2 < top))
  subset2 <- which( subset2.x & subset2.y)
  xs2 <- xs2[subset2]; ys2 <- ys2[subset2]
  
  ns <- binit2(xs1, ys1, xs2, ys2, nbins, r)
  rs <- r * (0:(nbins-1))                 #starting radius of each annulus.
  area <- l * w
  npts1 <- length(xs1); npts2 <- length(xs2);
  npts <- (npts1 + npts2)*0.5;
  ## since this could be either xs1 or xs2.

  fs <- drpcorrections(rs, l, w)

  density <- npts / area

  areas   <- pi * r^2 * ( (2*(1:nbins))-1);# Areas of each annulus (eq 1).
  areas   <- areas * fs                     #correction factor.
  lambdas <- npts * density * areas

  ds <- (ns / areas) / npts

  res <- drpeffrad(lambdas, ns,  npts, density)
  effrad <- res$r
  if (effrad < 0) {
    ## have a back-up in case we can't compute effective radius.
    effrad <- 1; p <- 0; maxr <- 0; k <- 0; Dc <- 0;
  } else {
    res <- drppackingfactor(effrad, density);
    p <- res$p; maxr <- res$maxr
    res <- drpreliability(density, area, r);
    k <- res$k; Dc <- res$Dc;
  }

  res <- list(effrad = effrad, p = p, maxr = maxr, k = k, Dc = Dc,
              ds =ds, density=density, n1=npts1, n2=npts2,
              nbins=nbins,r=r)
  class(res) <- "sjedrp"

  res
}

plot.sjedrp <- function (x) {
  ## Plot the results of the density recovery profile.
  hts <- x$ds
  last.bin <- x$nbins * x$r
  plot.label <- paste("eff rad", signif(x$effrad,3),
                      "packing", signif(x$p,3),
                      "reliability", signif(x$k,3))
  ##names(hts) <- x$rs
  barplot(hts, col="gray",space=0, width=x$r, xlim=c(0,last.bin),
          main=plot.label)
  lines( c(0,last.bin), c(x$density, x$density))
  axis(1, at=c(0, last.bin))
  lines( c(x$effrad, x$effrad), c(0, x$density))
  points( c(x$maxr), c(0),pch='|')
}

drppackingfactor <- function (effrad, d) {
  ## Return the packing factor P and the maximum radius MAX_R.
  maxr <-  sqrt(sqrt(4/3)/d)          # (eq 19)
  p <- (effrad/maxr) ^2               # (eq 21)
  if ( (p<0.0) ||(p>1.0))
    warning(paste("packing factor", p,"not in range [0,1]"))
  
  list (p = p, maxr = maxr)
}


drpeffrad <- function (lambdas, ns, n, d) {
  ## Return the effective radius.
  ## If the effective radius cannot be computed, it returns a negative value.
  diffs <- lambdas - ns;
  fracs <-  ns / lambdas;
  
  ##negs <- which(diffs < 0);
  
  ## Rather than finding first value when lambda - n < 0, we instead see
  ## when the value n/lambda > 1.  This should be equivalent, but
  ## sometims rounding errors mean that it never finds such a point.  In
  ## that case, we take the largest value of n/lambda and assume that
  ## bin is the one where the histogram first crosses the density line.
  
  negs <-  which(fracs > 1.0);
  ##browser()
  if (length(negs) == 0) {
    v <- max(fracs); firstnegative <- which(fracs == v)
    warning(paste("no negative elements in drpeffrad. Closest"
                  ,v, "at position", firstnegative))
    volume  <- sum( diffs[1:(firstnegative-1)])/n; # (eq 7)
    r <- sqrt(volume/(pi*d));	# (eq 9)
  }
  else {
    ## We have negative values, so we can just sum the difference
    ## values until just before the first negative entry.
    firstnegative <- negs[1];
    if (firstnegative == 1) {
      ## If the firstnegative element is in bin 1, this means that bin1
      ## is above average -- this can happen (and perhaps is indicative of
      ## clustering), in which case I think volume, and hence effective radius
      ## should be zero.
      ##warning("firstnegative should be greater than 1")
      volume <- 0
    } else {
      volume  <- sum( diffs[1:(firstnegative-1)])/n; # (eq 7)
    }
    r <- sqrt(volume/(pi*d));	# (eq 9)
  }
  
  list( r=r, first = firstnegative)
}


drpreliability <- function (density, area, r)
  {
    ## Return reliability factor K and critical density Dc.
    Dc <- 1.0 / (sqrt( area* pi) * r);    # (eq 11)
    k  <- density/ Dc;                    # (eq 13)
    list (k = k, Dc = Dc)
  }


binit2 <- function (xs1, ys1, xs2, ys2, nbins, r)
{
  ## Bin the cross-correlation distances.  This is a wrapper around
  ## a C routine that does all the hard work.
  npts1 <- length(xs1)
  npts2 <- length(xs2)
  z <- .C("drp_bin_it_r",
          as.double(xs1), as.double(ys1), as.integer(npts1),
          as.double(xs2), as.double(ys2), as.integer(npts2),
          as.integer(nbins), as.integer(r),
          ## create memory to store return values.
          ns = integer(nbins),
          PACKAGE = "sjedrp"
          )
  z$ns
}



drpcorrections <- function(rs, l, w)
{       
  ## Return the correction factors.
  ## TODO: check that region `a' is always non-negative, else these
  ## corrections shouldn't be used.
  ## (eq 32)
  fs <- 1 - (2*rs *(l+w)/( pi * l * w)) + ( (rs * rs)/ (pi * l * w))
  fs
}

drp.makestf <- function (x, y, file) {
  ## Make an STF file suitable for reading in with Bob Rodieck's program.
  npts <- length(x)
  stopifnot(npts == length(y))
  type <- numeric(length=npts)+1
  z <- numeric(length=npts)
  d <- data.frame(type=type, x=x, y=y, z=z, tag=z, br=z)
  write.table(d, file=file,
              eol="\r",                 #\r for the macintosh...
              sep="\t", quote=FALSE, row.names=FALSE)
}
