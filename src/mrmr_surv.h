#include <iostream>
#include <string>
#include <math.h> 
using namespace std;
#include <algorithm>
#include <ctime>
#include <cmath>
#include <vector>
#include <map>

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>


/* Entry points called from the R functions */
extern "C" 
{
	SEXP get_concordanceIndex_onevariable(SEXP msurv, SEXP ustrat, SEXP x2,SEXP cl2,SEXP st, SEXP se, SEXP weights, SEXP strat, SEXP N, SEXP outx, SEXP lenS, SEXP lenU);
	SEXP mrmr_cIndex( SEXP data, SEXP namat, SEXP cIndex, SEXP nvar, SEXP nsample, SEXP threshold);

}
double get_correlation(double data [],int ind_x, int ind_y, int size);
void build_mim_subset(double mim[],double data[],int nvar,int nsample, int subset [],int size_subset);
double returnConcordanceIndexC(int *msurv, int *ustrat, double *x2, int *cl2, double *st, int *se, double *weights, int *strat, int *N, int *outx, int *lenS, int *lenU);
