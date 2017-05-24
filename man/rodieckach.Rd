\name{rodieckach}
\alias{rodieckach}
\title{Rodieck ACh retinal cells}
\description{
  This is the file "Cells.stf" that Bob Rodieck distributes with MacDRP,
  for computing the density recovery profile.  The method should produce
  a packing factor of 0.26, effective radius of 32.4 um.  This DRP
  should be the same as that given in the manual and in figure 7b of the
  1991 paper.}
\usage{data(rodieckach)}
\format{Rodieck's old STF format for files; only columns used here are
  those titled x and y.}

\source{Rodieck(1991).}

\references{
  Rodieck RW (1991). The density recovery profile: A method for the
  analysis of points in the plane applicable to retinal studies. 
  Visual Neuroscience 6: 95-111.
}

\seealso{autodrp}
\examples{
## Try the Rodieck data set.
data(rodieckach)
xs.all <- rodieckach$x; ys.all <- rodieckach$y
subset <- which( (ys.all >=3000) & (ys.all <= 4000)); length(subset)
xs <- xs.all[subset]; ys <- ys.all[subset]
plot(xs, ys, asp=1, main="Rodieck example ACh data")
temp <- cbind(xs, ys); print(temp[26:30,]);
res <- autodrp(xs, ys, 20, 10, c(-373, 624, 3000, 4000)); plot(res)
stopifnot(abs(0.27 - res$p)< 0.01)
stopifnot(abs(32.7 - res$effrad)< 0.1)
}
\author{Stephen Eglen}
\keyword{datasets}
