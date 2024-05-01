#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>	/* for DLL types */


    



/* Octave code was in ~/mosaics/code/drp */

void cross_corr_r(double *x1s, double *y1s, int *pn1,
		  double *x2s, double *y2s, int *pn2,
		  int *pnbins, double *pr, int *pauto,
		  double *dxs,  double *dys, int *pk);

void drp_bin_it_r(double *x1s, double *y1s, int *pn1,
		  double *x2s, double *y2s, int *pn2,
		  int *pnbins, double *pr, int *pauto,
		  int   *ns);

R_CMethodDef cMethods[] = {
  {"cross_corr_r", (DL_FUNC) &cross_corr_r, 12,
   (R_NativePrimitiveArgType[12])
   {REALSXP, REALSXP, INTSXP,
    REALSXP, REALSXP, INTSXP,
    INTSXP,  REALSXP, INTSXP,
    REALSXP, REALSXP, INTSXP}},
  {"drp_bin_it_r", (DL_FUNC) &drp_bin_it_r, 10,
   (R_NativePrimitiveArgType[10])
   {REALSXP, REALSXP, INTSXP,
    REALSXP, REALSXP, INTSXP,
    INTSXP, REALSXP, INTSXP,
    INTSXP}},
  {NULL, NULL, 0}
};

void R_init_sjedrp(DllInfo *info) {
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

void drp_bin_it_r(double *x1s, double *y1s, int *pn1,
		  double *x2s, double *y2s, int *pn2,
		  int *pnbins, double *pr, int *pauto,
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
  double x1, y1, x2, y2, dx, dy, dist, r;
  int bin;

  int nbins, n1, n2, do_auto;
  nbins = *pnbins; n1 = *pn1; n2 = *pn2; r = *pr; do_auto=*pauto;
  
  /*Rprintf("nbins is %d\n", nbins); */

  for(i=0; i<nbins; i++)
    /* clear the bins first. */
    ns[i]=0;
  

  for(i=0; i<n1; i++) {
    x1 = x1s[i]; y1 = y1s[i];

    for (j=0; j<n2; j++) {
      if (do_auto && i==j) {
	;
      } else {
	x2 = x2s[j]; y2 = y2s[j];
	dx = x1 - x2; dy = y1 - y2;
	dist = sqrt( (double) ((dx*dx) + (dy*dy)));
	bin = (int)(dist/r);
	if (bin < nbins)
	  ns[bin]++;
      }
    }
  }
}

void cross_corr_r(double *x1s, double *y1s, int *pn1,
		  double *x2s, double *y2s, int *pn2,
		  int *pnbins, double *pr, int *pauto,
		  double *dxs,  double *dys, int *pk)

{
  /* First data set is (x1s, y1s) of length n1.  Second data set is
   * (x2s, y2s) of length n2.  For each point in data set 1, we find
   * the offset to each cell in data set 2.
   */
     
  int i,j;
  double x1, y1, x2, y2, dx, dy, dist2, r, maxdist;

  int nbins, n1, n2, do_auto, k;
  nbins = *pnbins; n1 = *pn1; n2 = *pn2; r = *pr; do_auto=*pauto;
  k = 0;
  maxdist = (r*nbins)*(r*nbins);
  
  for(i=0; i<n1; i++) {
    x1 = x1s[i]; y1 = y1s[i];

    for (j=0; j<n2; j++) {
      if (do_auto && (i==j)) {	/* TODO: check autocorrelation problem! */
	;
      } else {
	x2 = x2s[j]; y2 = y2s[j];
	dx = x1 - x2; dy = y1 - y2;
	dist2 =  (double) ((dx*dx) + (dy*dy));
	if (dist2 <= maxdist) {
	  dxs[k] = dx; dys[k] = dy;
	  k++;
	}
      }
    }
  }
  *pk=k;
}

