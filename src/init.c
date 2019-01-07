#include <stdlib.h> 
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#include "frailtypack.h"

static const R_FortranMethodDef FortEntries[] = {
    {"additive",             (DL_FUNC) &F77_SUB(additive),             62},
    {"cvpl",                 (DL_FUNC) &F77_SUB(cvpl),                 29},
    {"cvpl_logn",            (DL_FUNC) &F77_SUB(cvpl_logn),            29},
    {"cvpl_long",            (DL_FUNC) &F77_SUB(cvpl_long),            38},
	{"cvplnl",            	(DL_FUNC) &F77_SUB(cvplnl),                45},
    {"frailpenal",           (DL_FUNC) &F77_SUB(frailpenal),           58},
    {"frailpred_sha_nor_mc", (DL_FUNC) &F77_SUB(frailpred_sha_nor_mc),  5},
    {"joint",                (DL_FUNC) &F77_SUB(joint),                64},
    {"joint_longi",          (DL_FUNC) &F77_SUB(joint_longi),          62},
	{"jointlonginl",          (DL_FUNC) &F77_SUB(jointlonginl),	       62},
    {"joint_multiv",         (DL_FUNC) &F77_SUB(joint_multiv),         63},
	{"longiuninl",         (DL_FUNC) &F77_SUB(longiuninl),             31},
    {"nested",               (DL_FUNC) &F77_SUB(nested),               57},
    {"predict",              (DL_FUNC) &F77_SUB(predict),              44},
    {"predict_biv",          (DL_FUNC) &F77_SUB(predict_biv),          32},
	{"predictfam",          (DL_FUNC) &F77_SUB(predictfam),            31},
    {"predict_logn_sha",     (DL_FUNC) &F77_SUB(predict_logn_sha),     15},
    {"predict_recurr_sha",   (DL_FUNC) &F77_SUB(predict_recurr_sha),   19},
    {"predict_tri",          (DL_FUNC) &F77_SUB(predict_tri),          38},
	{"predicttrinl",	      (DL_FUNC) &F77_SUB(predicttrinl),  	   42},
    {"risque2",              (DL_FUNC) &F77_SUB(risque2),               6},
    {"survival_cpm",         (DL_FUNC) &F77_SUB(survival_cpm),          6},
    {"survival_cpm2",        (DL_FUNC) &F77_SUB(survival_cpm2),         6},
    {"survival_frailty",     (DL_FUNC) &F77_SUB(survival_frailty),      8},
    {"survival2",            (DL_FUNC) &F77_SUB(survival2),             6},
    {"survivalj_cpm2",       (DL_FUNC) &F77_SUB(survivalj_cpm2),        8},
	{"somme",                (DL_FUNC) &F77_SUB(somme),                 5},
    {"jointsurrogate",       (DL_FUNC) &F77_SUB(jointsurrogate),       55},
    {"test",                 (DL_FUNC) &F77_SUB(test),       3},
    {"surrosim",             (DL_FUNC) &F77_SUB(surrosim), 36},
    {"jointsurrokendall",    (DL_FUNC) &F77_SUB(jointsurrokendall), 24},
    {NULL, NULL, 0}
};

void 
attribute_visible 
R_init_frailtypack(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
	  R_forceSymbols(dll, TRUE);
}
