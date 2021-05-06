#include <R_ext/RS.h>

void 
F77_SUB(additive)(int *ns0, int *ng0, int *nst0, int *nz0, double *xmin10, double *xmin20, 
                        double *tt00, double *tt10, int *ic0, int *groupe0, int *nva0,
                        int *str0, double *vax0, int *interaction, int *ag0, int *noVar, 
                        int *maxiter0, int *irep10, int *correl0, int *np, double *b, 
                        double *coef, double *varcoef, double *varcoef2, double *rhoEnd, 
                        double *covEnd, double *varcovEnd, double *varSigma2, double *varTau2, 
                        int *ni, double *res, double *LCV, double *k0, double *x1Out, double *lamOut, 
                        double *xSu1, double *suOut, double *x2Out, double *lam2Out, double *xSu2, 
                        double *su2Out, int *typeof0, int *equidistant, int *nbintervR0, int *mt, 
                        int *ier, double *ddl, int *istop, double *shapeweib, double *scaleweib, 
                        int *mt1, int *trunc, double *ziOut, double *time, double *Resmartingale, 
                        double *frailtypred, double *frailtypred2, double *frailtyvar, 
                        double *frailtyvar2, double *frailtycov, double *linearpred, double *EPS);
                        
void 
F77_SUB(cvpl)(int *nobs, int *nsujet, int *groupe0, int *c0, int *cdc0, int *nva10, 
                        int *nva20, double *ve0, double *vedc0, int *typeof0, 
                        int *nz0, double *zi0, double *ttt0, double *tttdc0, 
                        int *nbintervR0, int *nbintervDC0, int *np, double *b, 
                        double *H_1, double *t00, double *t10, double *t0dc0, 
                        double *t1dc0, int *nt, double *valT, double *rl_cond, 
                        double *epoir, double *contribt, double *atrisk);
                        
void 
F77_SUB(cvpl_logn)(int *nobs, int *nsujet, int *groupe0, int *c0, int *cdc0, int *nva10, 
                        int *nva20, double *ve0, double *vedc0, int *typeof0, 
                        int *nz0, double *zi0, double *ttt0, double *tttdc0, 
                        int *nbintervR0, int *nbintervDC0, int *np, double *b, 
                        double *H_1, double *t00, double *t10, double *t0dc0, 
                        double *t1dc0, int *nt, double *valT, double *rl_cond, 
                        double *epoir, double *contribt, double *atrisk);
                        
void 
F77_SUB(cvpl_long)(int *ng0, int *nsujet0, int *nsujety0, int *groupe0, int *groupey0,
                        int *c0, int *cdc0, double *Y0, int *nva10, int *nva20, int *nva30, 
                        int *nb10, int *netar0, int *netadc0, int *link0, double *ve0, 
                        double *vedc0,double *velong0, double *matzy0, double *s_cag0,
                        int *s_cag_id0, int *typeof0, int *typeJoint0, int *nz0, 
                        double *zi0, int *np0, double *b0, double *H_1, double *t00, 
                        double *t10, double *t0dc0, double *t1dc0, int *nt, 
                        double *valT, double *rl_cond, double *epoir, 
                        double *contribt, double *atrisk);
                        
void
F77_SUB(cvplnl)(int *ng0, int *nsujet0, int *nsujety0, int *groupe0, int *groupey0,
                    int *c0, int *cdc0, double *Y0, int *nva10, int *nva20, int *nva30,
                    int *nva40, int *nb10, int *which_random0, double *box_cox0, int *netar0,
                    int *netadc0, int *link0, double *ve0,
                    double *vedc0, double *velong0, double *matzy0, double *s_cag0,
                    int *s_cag_id0, int *typeof0, int *nz0, double *zi0, int *np0,
                    double *b0, double *H_1, double *t00, double *t10, double *t0dc0,
                    double *t1dc0, int *nt, double *valT, double *rl_cond, double *epoir,
                    double *contribt, double *atrisk, int *GH, double *paGH0, double *weights0,
                    double *nodes0, int *nnodes_all0);
                        
void 
F77_SUB(frailpenal)(int *nsujetAux, int *ngAux, int *icenAux, int *nstAux, int *effetAux,
                        int *nzAux, double *axT, double *tt0Aux, double *tt1Aux,int *icAux, 
                        int *groupeAux, int *nvaAux, double *strAux, double *vaxAux, int *AGAux, 
                        int *noVar, int *maxitAux, int *irep1, int *np, double *b, 
                        double *H_hessOut, double *HIHOut, double *resOut, double *LCV,
                        double *xTOut, double *lamTOut, double *xSuT, double *suTOut, 
                        int *typeof0, int *equidistant, int *nbintervR0, int *mt,
                        int *ni, int *cpt, int *ier, double *k0, double *ddl, int *istop, 
                        double *shapeweib, double *scaleweib, int *mt1, double *ziOut, 
                        double *Resmartingale, double *martingaleCox, double *frailtypred, 
                        double *frailtyvar, double *frailtysd, double *linearpred, 
                        double *time, int *intcensAux, double *ttUAux,int *logNormal0,
                        int *timedep0, int *nbinnerknots0, int *qorder0,int *filtretps0, 
                        double *BetaTpsMat, double *EPS, int *nbgh);
                        
void 
F77_SUB(frailpred_sha_nor_mc)(int *np0, double *frailtypred, double *sig20, double *res10, 
                        int *nig0);
                        
void 
F77_SUB(joint)(int *nsujet0, int *ngrp, int *strAux, int *lignedc0, int *nz0, 
                        double *axT, double *tt00, double *tt10, int *ic0, int *groupe0, 
                        int *groupe00, int *fam0, double *tt0dc0, double *tt1dc0, int *icdc0, 
                        double *tempdc, int *icdc00, int *nva10, double *vax0, int *nva20, double *vaxdc0, 
                        double *vaxdc00, int *noVar, double *wtsvec0, int *maxit0, int *np, double *b, 
                        double *H_hessOut, double *HIHOut, double *resOut, double *LCV, double *x1Out, 
                        double *lamOut, double *xSu1, double *suOut, double *x2Out, 
                        double *lam2Out, double *xSu2, double *su2Out, int *typeofequidist, 
                        int *nbinterv0, int *mtaille, int *counts, 
                        int *IerIstop, double *paraweib, double *MartinGales, 
                        double *linearpred, double *linearpreddc, double *ziOut, 
                        double *time, double *timedc, double *linearpredG, int *typeJoint0, 
                        int *intcens0, int *indices0, double *ttU0, int *ordretmp, int *initialize,
                        int *logNormal0, int *paratps, int *filtretps0, double *BetaTpsMat, 
                        double *BetaTpsMatDc, double *EPS, int *nbgauss);    
                
void 
F77_SUB(joint_longi)(const int *vectnsujet, const int *ngnzag, double *k0, 
                     double *tt00, double *tt10, const int *ic0, const int *groupe0, 
                     double *tt0dc0, double *tt1dc0, const int *icdc0, const int *link0, 
                     double *yy0, double *bb0, const int *groupey0, const int *groupeb0, 
                     const int *vectnb0, double *fixed_binary0, double *matzy0, double *matzb0, const double *cag0, 
                     const int *vectnvar, const double *vax0, const double *vaxdc0, 
                     const double *vaxy0, const double *vaxb0, const int *novar, int *maxit0, 
                     int *np, const int *neta0, double *b, double *h_hessout, double *hihout, 
                     double *resout, double *lcv, double *x1out, double *lamout, double *xsu1, 
                     double *suout, double *x2out, double *lam2out, double *xsu2, double *su2out, 
                     int *typeof0, const int *equidistant, const int *mtaille, int *counts, 
                     int *ier_istop, double *paraweib, double *martingales, double *reslongi, 
                     double *pred_y0, const int *glmlog0, const int *positionvartime, 
                     const int *numinterac, double *linearpred, double *linearpreddc, 
                     double *ziout, const int *paratps, const int *filtretps0, double *betatpsmat, 
                     double *betatpsmatdc, double *betatpsmaty, double *eps, const int *gh, 
                     const double *pagh);

void
F77_SUB(jointlonginl)(int *nsujet0, int *nsujety0, int *ng0, int *nz0, double *k0, double *tt00,
                        double *tt10, int *ic0, int *groupe0, double *tt0dc0, double *tt1dc0,
                        int *icdc0, int *link0, double *yy0, int *groupey0, int *nb0, int *which_random0,
                        double *box_cox0, double *matzy0, double *cag0, int *nva10, double *vax0,
                        int *nva20, double *vaxdc0, int *nva30, int *nva40, double *vaxy0, int *noVar,
                        int *ag0, int *maxit0, int *np, int *neta0, double *b, double *H_hessOut,
                        double *HIHOut, double *resOut, double *LCV, double *x1Out, double *lamOut,
                        double *xSu1, double *suOut, double *x2Out, double *lam2Out, double *xSu2,
                        double *su2Out, int *typeof0, int *equidistant, int *mtaille, int *counts,
                        int *ier_istop, double *paraweib, double *ziOut, double *EPS,
                        int *GH, double *paGH0, double *b_pred0, int *effet0, int *indic_alpha0,
                        double *weights0, double *nodes0, int *nnodes_all0, int *RE_which0);
                        
void 
F77_SUB(joint_multiv)(int *nobsEvent, int *nz0, double *k0, double *tt00, double *tt10, 
                        double *tt0meta0,double *tt1meta0, int *ic0, int *icmeta0,
                        int *groupe0, int *groupe0meta, int *groupe0dc, double *tt0dc0, 
                        double *tt1dc0, int *icdc0, int *nbvar, double *vax0, double *vaxmeta0, 
                        double *vaxdc0, int *noVarEvent, int *maxIteration, int *initialize, 
                        int *np, double *b, double *H_hessOut, double *HIHOut, double *resOut, 
                        double *LCV, int *critCV, double *x1Out, double *lamOut, double *xSu1, 
                        double *suOut, double *x2Out,    double *lam2Out, double *xSu2, 
                        double *su2Out, double *x3Out, double *lam3Out, double *xSu3, 
                        double *su3Out, int *typeof0, int *equidistant0, int *nbIntervEvent, 
                        int *mtEvent, int *ni, int *cptEvent, double *shape_weib, double *scale_weib, 
                        int *mt1Event, int *irep, int *ag0, double *ResMartingaleEvent, 
                        double *frailtyEstimates, double *linearpred, double *linearpreddc, 
                        double *linearpredM, double *ziOut1, double *ziOutdc, double *ziOutmeta, 
                        double *time, double *timedc, double *timeM);
                        
void
F77_SUB(longiuninl)(int *nsujety0, int *ng0, double *yy0, int *groupey0, int *nb0,
                        int *which_random0, double *box_cox0, double *matzy0, double *cag0,
                        int *nva30, int *nva40, double *vaxy0, int *noVar, int *maxit0,
                        int *np, double *b, double *H_hessOut, double *HIHOut, double *resOut,
                        double *LCV, int *counts, int *ier_istop, double *EPS,
                        int *GH, double *paGH, double *b_pred, double *weights0,
                        double *nodes0, int *nnodes_all0, int *initialGH, double *axT);
                        
void 
F77_SUB(nested)(int *ns0, int *ng0, int *nssgbyg0, int *nst0, int *nz0, double *axT, 
                        double *tt00, double *tt10, int *ic0, int *groupe0, int *ssgroupe0, 
                        int *nva0, double *str0, double *vax0, int *AG0, int *noVar, 
                        int *maxiter0, int *irep1, int *np, int *maxngg, double *b, 
                        double *H_hessOut, double *HIHOut, double *resOut, double *LCV, 
                        double *x1Out, double *lamOut, double *xSu1, double *suOut, 
                        double *x2Out, double *lam2Out, double *xSu2, double *su2Out, 
                        int *typeof0, int *equidistant, int *nbintervR0, int *mt, int *ni,
                        int *cpt, int *ier, double *k0, double *ddl, int *istop, 
                        double *shapeweib, double *scaleweib, int *mt1, double *ziOut, 
                        double *time, double *Resmartingale, double *frailtypred, 
                        double *frailtypredg, double *frailtyvar, double *frailtyvarg, 
                        double *frailtysd, double *frailtysdg, double *linearpred, double *EPS, int *nbgl);
                        
void 
F77_SUB(predict)(int *np, double *b, int *nz, int *nbintervR, int *nbintervDC, int *nva1, 
                        int *nva2, int *nst, int *typeof0, int *typevent, double *zi, 
                        double *HIHOut, double *time, double *timedc,int *ntimeAll, 
                        int *npred0, double *predTime, double *window, double *predtimerec, 
                        int *nrec0, double *vaxpred0, double *vaxdcpred0, double *predAll1, 
                        double *predAll2, double *predAll3, double *predAll1R, 
                        double *predAlllow1, double *predAllhigh1, double *predAlllow2, 
                        double *predAllhigh2, double *predAlllow3, double *predAllhigh3, 
                        double *predAlllow1R, double *predAllhigh1R, int *icproba, 
                        int *nsample, int *intcens, double *trunctime, double *lowertime, 
                        double *uppertime, int *movingwindow, double *timeAll, int *modeltype, int *indic_alpha);
                        
void 
F77_SUB(predict_biv)(int *np, double *b, int *nz, int *nva20, int *nva30, int *nb_re0, int *nzyd, 
                        int *link0, int *nst, int *typeof0, double *zi0, double *HIHOut,    
                        int *ntimeAll, int *npred0, double *predTime, double *window, 
                        int *nrec0, double *yy0, double *vaxdcpred0, double *vaxypred0, 
                        int *groupey, int *uniGroupey, int *nsujety, double *predAll1, 
                        double *predAlllow1, double *predAllhigh1, int *icproba, 
                        int *nsample, int *movingwindow, double *timeAll, int *s_cag_id0, double *s_cag0);
                        
void 
F77_SUB(predictfam)(int *np, double *b, int *nz, int *nva1, int *nva2, int *nst, int *typeof0, 
                        double *zi, double *HIHOut,    int *indID, double *tt1T, double *tt1dcT, 
                        int *icdctime, int *ntimeAll, int *nsujet, int *npred0, double *window, 
                        int *nrec0, int *nrec, int *nrecT,double *vaxpred0, double *vaxdcpred0, 
                        int *icproba, int *nsample,    double *predAll, double *predIClow, 
                        double *predIChigh, double *frailfam0, double *frailind0, double *pred, int *indic);
                        
void 
F77_SUB(predict_logn_sha)(int *npred0, double *surv_s, double *surv_t, double *betapred, 
                        double *sigma2, double *predAll, int *icproba, int *ntimeAll, 
                        int *nsample, double *sig2alea, double *surv_smc, double *surv_tmc,
                        double *betapredmc, double *predAlllow, double *predAllhigh);
                        
void 
F77_SUB(predict_recurr_sha)(int *LogN, int *npred0, double *surv_s, double *surv_t, 
                        double *surv_r, double *betapred, double *var, double *predAll, 
                        int *nreci, int *ntimeAll, int *icproba, int *nsample, 
                        double *varalea, double *surv_smc, double *surv_tmc, 
                        double *surv_rmc, double *betapredmc, double *predAlllow, 
                        double *predAllhigh);
                        
void 
F77_SUB(predict_tri)(int *np, double *b, int *nz, int *nva10, int *nva20, int *nva30, 
                        int *nb_re0, int *nzyr, int *nzyd, int *link0, int *nst, 
                        int *typeof0, double *zi0, double *HIHOut, int *ntimeAll, 
                        int *npred0, double *predTime, double *window, 
                        double *predtimerec, int *nrec0, int *nrecy0, double *yy0, 
                        double *vaxpred0, double *vaxdcpred0, double *vaxypred0, 
                        int *groupey, int *uniGroupey, int *nsujety, int *nsujet, 
                        double *predAll1, double *predAlllow1, double *predAllhigh1,
                        int *icproba, int *nsample, int *movingwindow, double *timeAll, 
                        int *s_cag_id0, double *s_cag0);
                        
void
F77_SUB(predicttrinl)(int *np, double *b, int *nz, int *nva10, int *nva20, int *nva30, 
                        int *nva40, int *nb_re0, int *random_which0, double *box_cox0,
                        int *nzyr, int *nzyd, int *link0, int *nst, 
                        int *typeof0, double *zi0, double *HIHOut, int *ntimeAll, 
                        int *npred0, double *predTime, double *window, 
                        double *predtimerec, int *nrec0, int *nrecy0, double *yy0, 
                        double *matzy0, double *vaxpred0, double *vaxdcpred0, double *vaxypred0, 
                        int *groupey, int *uniGroupey, int *nsujety, int *nsujet, 
                        double *predAll1, double *predAlllow1, double *predAllhigh1,
                        int *icproba, int *nsample, int *movingwindow, double *timeAll, 
                        int *s_cag_id0, double *s_cag0);
                        
void 
F77_SUB(risque2)(double *t, double *the0, int *nz, double *zi0, 
                        double *lam, int *nst);
void 
F77_SUB(survival_cpm)(double *t, double *b, int *nst, int *nbintervR, double *time,double *surv);

void 
F77_SUB(survival_cpm2)(double *t, double *b, int *nst, int *nbintervR, 
                        double *time, double *surv);
                        
void 
F77_SUB(survival_frailty)(double *t, double *the_s, double *the1_s, int *nz, double *zi_s, 
                        double *su, double *lam, int *nst);
                        
void 
F77_SUB(survival2)(double *t, double *the0, int *nz, double *zi0, double *su, int *nst);

void 
F77_SUB(survivalj_cpm2)(double *t, double *b, int *nst, int *nbintervR, int *nbintervDC, 
                        double *time, double *timedc, double *surv);
void F77_SUB(somme)(double *ab, double *s, int *nboot, int *nbr, double *vOut);
             
void
F77_SUB(jointsurrogate)(int *nsujet1, int *ng, int *ntrials1, int *maxiter, int *nst,
             int *nparamfrail, int *indice_a_estime, int *param_risque_base, int *nbrevar,
             int *filtre0, double *donnees, double *death, double *p, double *prop_i, int *n_sim1,
             double *EPS2, double *kappa0, double *vect_kappa, int *logNormal, int *nsim_node, int *Param_kendall_boot,
             int *vrai_val_init, double *param_init, double *revision_echelle,
             int *random_generator0, int *sujet_equi, double *prop_trait, double *paramSimul,
             int *autreParamSim, double *fichier_kendall, double *fichier_R2, double *param_estimes, int *sizeVect,
             double *b, double *H_hessOut, double *HIHOut, double *resOut, double *LCV,
             double *x1Out, double *lamOut, double *xSu1, double *suOut, double *x2Out,
             double *lam2Out, double *xSu2, double *su2Out, int *ni, int *ier, int *istop,
             double *ziOut, int *affiche_itter, double *Varcov, double *dataHessian, double *dataHessianIH, double *datab,
             double *vbetast, double *vbetastinit);

void
F77_SUB(surrosim) (double *don_simul, double *don_simulS1, int *n_obs, 
		  int *n_col, int *lognormal, int *affiche_stat,
		  double *vrai_theta, int *ng, int *ver, double *truealpha, 
		  double *propC, double *cens_A, double *gamma1, 
		  double *gamma2, double *theta2, double *lambda_S, 
		  double *nu_S, double *lambda_T, double *nu_T,
          double *betas, double *betat, int *n_essai, 
		  double *rsqrt, double *sigma_s, double *sigma_t, 
		  double *p, double *prop_i, double *gamma, double *alpha,
          int *frailt_base, int *random_generator0, int *aleatoire, 
		  int *nbre_sim, int *graine, int *nbre_don_non_cons, 
		  int *param_weibull0, double *thetacopule, int *filtre, 
          int *filtre2, int *type_joint_simul, int *pfs);

void
F77_SUB(jointsurrokendall)(
          double *theta, double *gamma, double *alpha, double *eta,
          int *adaptative, int *npg, int *ndim,
          double *ui_chap_Essai, double *invBi_chol_Essai_k,double *xx1, double *ww1,
          double *sigma_v, int *methodeInt, int *N_MC_kendall, int *method_int_kendal, 
          int *random_gen, int *aleatoire, int *nbre_sim, int *seed, double *tau_kendal_00, 
          double *tau_kendal_01, double *tau_kendal_10, double *tau_kendal_11, double *ss
  );
  

void F77_SUB(test)(double *donnee, int *nrow, double *som);


void
  F77_SUB(frailpenalgen)(int *nsujetAux, int *ngAux, int *icenAux, int *nstAux, int *effetAux,
          int *nzAux, double *axT, double *tt0Aux, double *tt1Aux,int *icAux,
          int *groupeAux, int *nvaAux, double *strAux, double *vaxAux, int *AGAux,
          int *noVar, int *maxitAux, int *irep1, int *np, double *b,
          double *H_hessOut, double *HIHOut, double *resOut, double *LCV,
          double *xTOut, double *lamTOut, double *xSuT, double *suTOut,
          int *typeof0, int *equidistant, int *nbintervR0, int *mt,
          int *ni, int *cpt, int *ier, double *k0, double *ddl, int *istop,
          double *shapeweib, double *scaleweib, int *mt1, double *ziOut,
          double *Resmartingale, double *martingaleCox, double *frailtypred,
          double *frailtyvar, double *frailtysd, double *linearpred,
          double *time, int *intcensAux, double *ttUAux,int *logNormal0,
          int *timedep0, int *nbinnerknots0, int *qorder0,int *filtretps0,
          double *BetaTpsMat, double *EPS, int *nbgh, int *familyrisk);
		  
void
F77_SUB(jointgen)(int *nsujet0, int *ngrp, int *strAux, int *lignedc0, int *nz0,
						double *axT, double *tt00, double *tt10, int *ic0, int *groupe0,
						int *groupe00, int *fam0, double *tt0dc0, double *tt1dc0, int *icdc0,
						double *tempdc, int *icdc00, int *nva10, double *vax0, int *nva20, double *vaxdc0,
						double *vaxdc00, int *noVar, double *wtsvec0, int *maxit0, int *np, double *b,
						double *H_hessOut, double *HIHOut, double *resOut, double *LCV, double *x1Out,
						double *lamOut, double *xSu1, double *suOut, double *x2Out,
						double *lam2Out, double *xSu2, double *su2Out, int *typeofequidist,
						int *nbinterv0, int *mtaille, int *counts,
						int *IerIstop, double *paraweib, double *MartinGales,
						double *linearpred, double *linearpreddc, double *ziOut,
						double *time, double *timedc, double *linearpredG, int *typeJoint0,
						int *intcens0, int *indices0, double *ttU0, int *ordretmp, 
						int *init_and_logN,
						int *paratps, int *filtretps0, double *BetaTpsMat,
						double *BetaTpsMatDc, double *EPS, int *nbgauss, 
						int*familyrisk);
						
