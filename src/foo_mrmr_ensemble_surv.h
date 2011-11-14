#include <iostream>
#include <string>
#include <math.h> 
using namespace std;
#include <algorithm>
#include <ctime>
#include <cmath>
#include <vector>
#include <map>
#include "tree.h"

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

/* Entry points called from the R functions */
extern "C" 
{
//	SEXP get_concordanceIndex_onevariable(SEXP msurv, SEXP ustrat, SEXP x2,SEXP cl2,SEXP st, SEXP se, SEXP weights, SEXP strat, SEXP N, SEXP outx, SEXP lenS, SEXP lenU);
//	SEXP mrmr_cIndex( SEXP data, SEXP namat, SEXP nvar, SEXP nsample, SEXP threshold,SEXP msurv, SEXP ustrat,SEXP cl2,SEXP st, SEXP se, SEXP weights, SEXP strat, SEXP N, SEXP outx, SEXP lenS, SEXP lenU);
	SEXP mrmr_cIndex_ensemble_remove(SEXP data, SEXP namat, SEXP maxparents, SEXP nvar, SEXP nsample, SEXP predn, SEXP npredn, SEXP rep_boot, SEXP maxnsol, SEXP threshold,SEXP msurv, SEXP ustrat,SEXP cl2,SEXP st, SEXP se, SEXP weights, SEXP strat, SEXP N, SEXP outx, SEXP lenS, SEXP lenU);

}

double get_correlation_ensemble(double data [],int ind_x, int ind_y, int size);
void build_mim_cIndex_subset(double mim[],double data[],int nvar,int nsample, int subset [],int size_subset,int *msurv, int *ustrat, int *cl2, double *st, int *se, double *weights, int *strat, int *N, int *outx, int *lenU);
double returnconcordanceIndexC(int *msurv, int *ustrat, double *x2, int *cl2, double *st, int *se, double *weights, int *strat, int *N, int *outx, int lenS, int *lenU);
void remove_childless_nodes( tree<int>& res, tree<double>&res_mean, int max_elements_tmp);
void bootstrap_tree(tree<int>& res,tree<double>& res_mrmr, double data[],int namat[], int nsamples,int n, int rep_boot,  int *msurv, int *ustrat, double *x2, int *cl2,
					double *st, int *se, double *weights, int *strat, int *N, int *outx, int lenS, int *lenU);
void bootstrap_mrmr(double &mean, double &sd, double data[],int namat[],int size, int rep_boot, int size_boot,int nsamples, int var_target, int var_interest, int nprev_sel,int* var_ind, int *msurv, int *ustrat, double *x2, int *cl2,
					double *st, int *se, double *weights, int *strat, int *N, int *outx, int lenS, int *lenU);
double mrnet_onegene(double mim [], int size, int nbvar,int *var_ind,int var_target, int var_interest);
int verify_equivalentset_nparents (tree<int>& tr, tree<int>::pre_order_iterator it, tree<int>::pre_order_iterator end,tree<double>& tr_mrmr, int maxnsol);
