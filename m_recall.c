#ifdef RCSID
static char rcsid[] = "$Header: /home/smart/release/src/libevaluate/trvec_trec_eval.c,v 11.0 1992/07/21 18:20:35 chrisb Exp chrisb $";
#endif

/* Copyright (c) 2008
*/

#include "common.h"
#include "sysfunc.h"
#include "trec_eval.h"
#include "functions.h"
#include "trec_format.h"

/*      "    Recall at cutoffs\n\
    Recall (relevant retrieved / relevant) measured at various doc level\n\
    cutoffs in the ranking. If the cutoff is larger than the number of docs\n\
    retrieved, then it is assumed nonrelevant docs fill in the rest.\n\
    Recall is a fine single topic measure, but does not average well.\n\
    Cutoffs must be positive without duplicates\n\
    Default param: trec_eval -m recall.5,10,15,20,30,100,200,500,1000\n",
*/

int 
te_calc_recall (const EPI *epi, const REL_INFO *rel_info,
		const RESULTS *results, const TREC_MEAS *tm, TREC_EVAL *eval)
{
    LONG_PARAMS *cutoffs = (LONG_PARAMS *) tm->meas_params;
    long cutoff_index = 0;
    long i;
    RANK_REL rank_rel;
    long rel_so_far = 0;

    if (UNDEF == form_ordered_rel (epi, rel_info, results, &rank_rel))
	return (UNDEF);

    for (i = 0; i < rank_rel.num_ret; i++) {
	if (i == cutoffs->param_values[cutoff_index]) {
	    /* Calculate previous cutoff threshold.
	     Note cutoffs guaranteed to be positive by init_meas */
	    eval->values[tm->eval_index + cutoff_index].value =
		(double) rel_so_far / (double) rank_rel.num_rel;
	    if (++cutoff_index == cutoffs->num_params)
		break;
	}
	if (rank_rel.results_rel_list[i] >= epi->relevance_level)
	    rel_so_far++;
    }
    /* calculate values for those cutoffs not achieved */
    while (cutoff_index < cutoffs->num_params) {
	eval->values[tm->eval_index + cutoff_index].value =
	    (double) rel_so_far / (double) rank_rel.num_rel;
	cutoff_index++;
    }
    return (1);
}