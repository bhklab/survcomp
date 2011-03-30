#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "survcomp.h"

R_CMethodDef cMethods[] = {
    {"concordanceIndexC", (DL_FUNC)&concordanceIndexC, 16},
    {NULL, NULL, 0}
};

//{INTSXP,INTSXP,REALSXP,INTSXP,REALSXP,INTSXP,REALSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP}

void R_init_concordanceIndexC(DllInfo *info) {
    R_registerRoutines(info,cMethods,NULL,NULL,NULL);
}

void R_unload_concordanceIndexC(DllInfo *info) { }

void concordanceIndexC(int *msurv, int *ustrat, double *x2, int *cl2,
		double *st, int *se, double *weights, int *strat, int *N, int *outx,
		int *ch, int *dh, int *uh, int *rph, int *lenS, int *lenU) {
	
	int lenUstrat = *lenU;
	int lenStrat = *lenS;
	int Ns_old = 0;
	int Ns = 0;	
    for(int s=0; s < lenUstrat; s++) {
    	int ixs[lenStrat];
    	for(int i =0; i < lenStrat; i++){
    		ixs[i] = 0;
    		if(strat[i] == ustrat[s]){
    			ixs[i] = 1;
    		} else {
    			ixs[i] = 0;
    		}
    	}
    	Ns_old += Ns;
    	Ns = 0;
        for(int i=0; i < lenStrat; i++){
        	if(ixs[i] == 1){
        		Ns++;
        	}
        }
        double xs[Ns];
        int c = 0;
        for(int i=0; i < lenStrat; i++){
        	if(ixs[i] == 1){
        		xs[c] = x2[i];
        		c++;
        	}
        }
        int cls[Ns];
        c = 0;
        for(int i=0; i < lenStrat; i++){
        	if(ixs[i] == 1){
        		cls[c] = cl2[i];
        		c++;
        	}
        }
        double sts[Ns];
        c = 0;
        for(int i=0; i < lenStrat; i++){
        	if(ixs[i] == 1){
        		sts[c] = st[i];
        		c++;
        	}
        }        
        int ses[Ns];
        c = 0;
        for(int i=0; i < lenStrat; i++){
        	if(ixs[i] == 1){
        		ses[c] = se[i];
        		c++;
        	}
        }
        double weightss[Ns];
        c = 0;
        for(int i=0; i < lenStrat; i++){
        	if(ixs[i] == 1){
        		weightss[c] = weights[i];
        		c++;
        	}
        }
        double chs[Ns];
        double dhs[Ns];
        double uhs[Ns];
        double rphs[Ns];
        for (int h=0; h < Ns; h++) {       	
        	double chsj, dhsj, uhsj, rphsj = 0;
            for (int j=0; j < Ns; j++) {
                double whj = weightss[h] * weightss[j];
                if ((*msurv == 1 && (sts[h] < sts[j] && ses[h] == 1)) || (*msurv == 0 && cls[h] > cls[j])) {
                  rphsj = rphsj + whj;
                  if (xs[h] > xs[j]) {  
                    chsj = chsj + whj;
                  } else if (xs[h] < xs[j]) {
                    dhsj = dhsj + whj;
                  }
                  else {
                    if (*outx == 1) {
                      uhsj = uhsj + whj;
                    } else {
                      dhsj = dhsj + whj;
                    }
                  }
                }    
                if ((*msurv == 1 && (sts[h] > sts[j] && ses[j] == 1)) || (*msurv == 0 && cls[h] < cls[j])) { 	
                  rphsj = rphsj + whj;
                  if (xs[h] < xs[j]) {
                    chsj = chsj + whj;
                  }
                  else if (xs[h] > xs[j]) {
                    dhsj = dhsj + whj;
                  }
                  else {
                    if (*outx == 1) {
                      uhsj = uhsj + whj;
                    } else {
                      dhsj = dhsj + whj;
                    }
                  }
                }
            }
            chs[h] = chsj;
            dhs[h] = dhsj;
            uhs[h] = uhsj;
            rphs[h] = rphsj;
            chsj = 0;
            dhsj = 0;
            uhsj = 0;
            rphsj = 0;
        }
        for(int i = 0; i < Ns; i++){
        	int pos = i + Ns_old;
        	ch[pos] = chs[i];
        	dh[pos] = dhs[i];
        	uh[pos] = uhs[i];
        	rph[pos] = rphs[i];
        }
    }
}
