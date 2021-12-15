*! Resets existing results to null or initializes to correct null specs if no previous results.

program _ddml_reset_model_results
	version 13

	syntax [if] [in] , mname(name)

	// reset previous results
	// assumes kfolds has been set
	mata: reset_model_results(`mname')

end

mata:

struct eqnStruct init_eqnStruct()
{
	struct eqnStruct scalar		e
	return(e)
}

void reset_model_results(		struct ddmlStruct m
								)

{
	pointer(struct eqnStruct) scalar p
	
	m.crossfitted	= 0
	
	kfolds			= m.kfolds

	for (i=1; i<=cols(m.eqnlist); i++) {
		p					= m.eqnlist[i]
		(*p).MSE			= J(0,1,0)
		(*p).MSE_h			= J(0,1,0)
		(*p).MSE0			= J(0,1,0)
		(*p).MSE1			= J(0,1,0)
		(*p).N				= J(0,1,0)
		(*p).N_h			= J(0,1,0)
		(*p).N0				= J(0,1,0)
		(*p).N0				= J(0,1,0)

		(*p).MSE_folds		= J(0,kfolds,0)
		(*p).MSE_h_folds	= J(0,kfolds,0)
		(*p).MSE0_folds		= J(0,kfolds,0)
		(*p).MSE1_folds		= J(0,kfolds,0)
		(*p).N_folds		= J(0,kfolds,0)
		(*p).N_h_folds		= J(0,kfolds,0)
		(*p).N0_folds		= J(0,kfolds,0)
		(*p).N1_folds		= J(0,kfolds,0)
		
		(*p).crossfitted	= 0
	}

}

end

/*
	real colvector		MSE
	real matrix			MSE_folds		// MSE by fold; col=fold, row=resample
	real colvector      MSE_h 			// (intended for LIE)
	real matrix			MSE_h_folds		// (intended for LIE)
	real colvector 		MSE0
	real colvector 		MSE1
	real matrix			MSE0_folds		// MSE by fold; col=fold, row=resample
	real matrix			MSE1_folds		// MSE by fold; col=fold, row=resample
	real colvector		N
	real matrix			N_folds			// sample size by fold; col=fold, row=resample
	real colvector		N_h				// (intended for LIE)
	real matrix         N_h_folds		// (intended for LIE)
	real colvector		N0
	real colvector		N1
	real matrix			N0_folds		// sample size by fold; col=fold, row=resample
	real matrix			N1_folds		// sample size by fold; col=fold, row=resample
*/
