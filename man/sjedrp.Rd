\name{sjedrp}
\alias{autodrp}
\alias{crossdrp}
\alias{plot.sjedrp}
\alias{sjedrp}
\alias{drp.makestf}
\title{Density recovery profile}
\description{
  Compute Rodieck's density recovery profile.
}
\usage{
autodrp(xs, ys, nbins, r, a)
crossdrp(xs1, ys1, xs2, ys2, nbins, r, a,auto)
\method{plot}{sjedrp}(x, scale = 1, title = NULL, mirror = FALSE,
     show.title = TRUE, ylab = "density", xlab = "distance", ...)
drp.makestf(x, y, file)
}
\arguments{
  \item{xs,xs1,xs2,x}{vector of x coordinates}
  \item{ys,ys1,ys2,y}{vector of y coordinates}
  \item{nbins}{Number of bins}
  \item{r}{Width of each bin}
  \item{a}{Bounds of region.  Either NULL if you want to make the
    boundary exactly fit all the data points; or a vector of length 4
    with elements (lowx, highx, lowy, highy).  Cells outside the bounds
    are rejected.  If a is omitted, packing factors may be slightly
    higher than when giving bounds.}
  \item{auto}{A boolean, FALSE by default.  This is set to TRUE if
    autodrp has called crossdrp to do the correlations, so the user
    won't need to worry about setting this.}
  \item{file}{File output (with Macintosh end of lines) for Bob
    Rodieck's program.  This has extra columns, as needed by pre-2001
    versions of his program.}
  \item{scale}{Scaling factor for multiplying values on y-axis.  For
    example, if scale=1e6, this converts the density from $um^{-2}$ to
    $mm^{-2}$, which is more readable.}
  \item{mirror}{A boolean, FALSE by default.  If true, the DRP is
    plotted symmetrically around the y-axis.}
  \item{title}{Title to place on top of the plot.}
  \item{show.title}{A boolean, TRUE by default.  If TRUE, show the
    title.}
  \item{xlab}{Label to show on the x axis.}
  \item{ylab}{Label to show on the y axis.}
  \item{...}{Extra arguments to pass to plotting functions.}
}
\value{
  \item{res}{A list with the following components:}
  \item{effrad}{the effective radius}
  \item{p}{the packing factor}
  \item{maxr}{the maximum radius}
  \item{k}{reliability factor}
  \item{Dc}{critical density}
}
\details{
  
  Follow the Rodieck prescription for finding the effective radius and
  the packing factor.  Note this method is very close to that of the
  Ripley method for computing Khat.

  \code{drp.makestf} takes a set of data points and will output a text
  file (with Macintosh end of lines) that can then be input into the
  MacDRP program from Bob Rodieck.  This serves as a useful comparison.
  His 2001 version of the program now accepts just two-column files, so
  this function is no longer required.
}
\references{
  Rodieck RW (1991). The density recovery profile: A method for the
  analysis of points in the plane applicable to retinal studies. 
  Visual Neuroscience 6: 95-111.
}
\seealso{
  \code{\link[splancs]{khat}} and other functions by Ripley and Diggle.
  See also \code{\link{rodieckach}} for the data and code to generate the DRP
  from Rodieck (1992).  Various data files have been tested, see
  \code{\link{sjedrp.d10}} for details.
}

\examples{
data(sjedrp.d20)
bb <- c(0, 700, 0, 700)
res <- autodrp(sjedrp.d20$x, sjedrp.d20$y, nbins=20, r=10, a=bb)
par(mfrow=c(2,1), las=1)
plot(res, scale=1e6)
plot(res, mirror=TRUE, scale=1e6)
title(xlab=expression(paste("distance (", mu, "m)")),
      ylab=expression(paste("density (cells/",mm^2, ")")))
plot(res, scale=1e6, ylim=c(0, 800))
}

\author{Stephen Eglen}
\keyword{spatial}
