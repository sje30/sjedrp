#include <R.h>
#include <S.h>			/* for seed_in, seed_out */

#define SMALL_VALUE 0.00001

/* If we are doing an AUTODRP, then data sets 1 and 2 are the same points.
 * In this case, when i==j the distance will be zero.  But we don't want
 * to include those points, so we reject very small distances.
 */

/* Octave code was in ~/mosaics/code/drp */

void drp_bin_it_r(Sfloat *x1s, Sfloat *y1s, int *pn1,
		  Sfloat *x2s, Sfloat *y2s, int *pn2,
		  int *pnbins, int *pr,
		  int   *ns);



void drp_bin_it_r(Sfloat *x1s, Sfloat *y1s, int *pn1,
		  Sfloat *x2s, Sfloat *y2s, int *pn2,
		  int *pnbins, int *pr,
		  int   *ns)
{
  /* First data set is (x1s, y1s) of length n1.  Second data set is
   * (x2s, y2s) of length n2.  For each point in data set 1, we find
   * the distance to each point in data set 2.  This distance is then
   * binned into one of NBINS values in NS, spaced R um apart.
   *
   * This is used by the DRP procedure and ensures that it is _quick_,
   * say for greater than 200 cells.  Another speed, maybe in R, is to
   * sort the points by their x-coordinate, so that we only have to
   * measure distance to a subset of cells.  This will take more
   * effort though.  This seems to work fine.
   */
     
  int i,j;
  Sfloat x1, y1, x2, y2, dx, dy, dist;
  int bin;

  int nbins, r, n1, n2;
  nbins = *pnbins; n1 = *pn1; n2 = *pn2; r = *pr;
  /*Rprintf("nbins is %d\n", nbins); */

  for(i=0; i<nbins; i++)
    /* clear the bins first. */
    ns[i]=0;
  

  for(i=0; i<n1; i++) {
    x1 = x1s[i]; y1 = y1s[i];

    for (j=0; j<n2; j++) {
      x2 = x2s[j]; y2 = y2s[j];
      dx = x1 - x2; dy = y1 - y2;
      dist = sqrt( (double) ((dx*dx) + (dy*dy)));
      bin = (int)(dist/r);
      if ((bin < nbins) && (dist > SMALL_VALUE)) {
	ns[bin]++;
	/*if (bin == 0) printf("bin 0: distance %f i %d j %d\n", dist, i, j);*/
      }
    }
  }
}

