/* 
   Copyright (c) 2008 - Chris Buckley. 

   Permission is granted for use and modification of this file for
   research, non-commercial purposes. 
*/
#include "common.h"
#include "sysfunc.h"
#include "trec_eval.h"
#include "functions.h"
#include "trec_format.h"
double log2(double x);

static int 
te_calc_dcg_cut (const EPI *epi, const REL_INFO *rel_info,
		  const RESULTS *results, const TREC_MEAS *tm, TREC_EVAL *eval);
static long long_cutoff_array[] = {1, 5, 10, 15, 20, 30, 100, 200, 500, 1000};
static PARAMS default_dcg_cutoffs = {
    NULL, sizeof (long_cutoff_array) / sizeof (long_cutoff_array[0]),
    &long_cutoff_array[0]};

/* See trec_eval.h for definition of TREC_MEAS */
TREC_MEAS te_meas_dcg_cut =
    {"dcg_cut",
     "    Discounted Cumulative Gain at cutoffs.\n\
    Compute a traditional DCG measure according to Jarvelin and\n\
    Kekalainen (ACM ToIS v. 20, pp. 422-446, 2002) at cutoffs.\n\
    See comments for dcg.\n\
    Gain values are the relevance values in the qrels file.  For now, if you\n\
    want different gains, change the qrels file appropriately.\n\
    Cutoffs must be positive without duplicates\n\
    Default params: -m dcg_cut.5,10,15,20,30,100,200,500,1000\n\
    Based on an implementation by Yifan\n",
     te_init_meas_a_float_cut_long,
     te_calc_dcg_cut,
     te_acc_meas_a_cut,
     te_calc_avg_meas_a_cut,
     te_print_single_meas_a_cut,
     te_print_final_meas_a_cut,
     (void *) &default_dcg_cutoffs, -1};

static int 
te_calc_dcg_cut (const EPI *epi, const REL_INFO *rel_info,
		  const RESULTS *results, const TREC_MEAS *tm, TREC_EVAL *eval)
{
    long  *cutoffs = (long *) tm->meas_params->param_values;
    long cutoff_index = 0;
    RES_RELS res_rels;
    double gain, sum;
    long cur_lvl, lvl_count;
    long i;
   
    if (UNDEF == te_form_res_rels (epi, rel_info, results, &res_rels))
	return (UNDEF);

    sum = 0.0;
    for (i = 0; i < res_rels.num_ret; i++) {
        if (i == cutoffs[cutoff_index]) {
            /* Calculate previous cutoff threshold.
             Note i guaranteed to be positive by init_meas */
            eval->values[tm->eval_index + cutoff_index].value = sum / i;
            if (++cutoff_index == tm->meas_params->num_params)
                break;
	    if (epi->debug_level > 0) 
		printf("dcg_cut: cutoff %ld dcg %6.4f\n", i, sum);
        }
	gain = res_rels.results_rel_list[i];
	if (gain > 0) {
	    /* Note: i+2 since doc i has rank i+1 */
        if (i==0) sum = 1;
        else sum += 1 / log2((double) (i+1));
	    if (epi->debug_level > 1) 
		printf("dcg_cut:%ld %3.1f %6.4f\n", i, gain, sum);
	}
    }
    /* calculate values for those cutoffs not achieved */
    while (cutoff_index < tm->meas_params->num_params) {
	eval->values[tm->eval_index + cutoff_index].value = sum / cutoffs[cutoff_index];
	if (epi->debug_level > 0) 
	    printf("dcg_cut: cutoff %ld dcg %6.4f\n",
		   cutoffs[cutoff_index], sum);
        cutoff_index++;
    }
    

    return (1);
}
