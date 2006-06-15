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
res <- autodrp(xs, ys, nbins, r, a)
res <- crossdrp(xs1, ys1, xs2, ys2, nbins, r, a,auto)
plot.sjedrp(res)
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
  \code{\link[spatial]{khat}} and other functions by Ripley and Diggle.
  See also \code{\link{rodieckach}} for the data and code to generate the DRP
  from Rodieck (1992).  Various data files have been tested, see
  \code{\link{sjedrp.d10}} for details.
}
\author{Stephen Eglen}
\keyword{spatial}
