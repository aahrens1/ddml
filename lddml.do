* Locals used in whichddml; set when compiled
local stata_version `c(stata_version)'
local born_date `c(born_date)'
local current_date `c(current_date)'

version 13
mata:
mata clear

void whichddml()
{
""
"Mata library for ddml and related programs,"
"compiled `current_date' under Stata " + "`stata_version'" + " born " + "`born_date'."
"authors AA/MS"
st_sclear()
st_global("s(stata_born_date)","`born_date'")
st_global("s(stata_version)","`stata_version'")
st_global("s(compiled_date)","`current_date")
}

// uniquely identified by vname = dep var in the equation
// equation structure: one for E[y|x], E[d|x] etc
struct eStruct {
	string scalar					vname		// name of variable to be orthogonalized
	real matrix						vtlist		// list of orthogonalized (learner) variables
	string scalar					shortstack	// name of shortstack variable
	real scalar						nlearners	// number of learners
	real scalar						lieflag		// =1 if LIE spec with two estimation strings etc.
	class AssociativeArray scalar	lrnAA		// AssociativeArray with all learners //
												// (keys=vtilde,object)
	class AssociativeArray scalar	resAA		// AssociativeArray with all learner results //
												// (keys=vtilde,object,rep)

}

// should perhaps make vname a required argument for this function
struct eStruct init_eStruct()
{
	struct eStruct scalar			d
//	class AssociativeArray scalar	A2, A3

	d.vname			= ""
	d.vtlist		= J(1,0,"")
	d.shortstack	= ""
	d.nlearners		= 0
	d.lieflag		= 0
	
	(d.lrnAA).reinit("string",2)
	(d.resAA).reinit("string",3)
	(d.lrnAA).notfound(NULL)
	(d.resAA).notfound(NULL)
	
	return(d)
}

// clear results from eStruct
void clear_equation_results(struct eStruct e)
{
	(e.resAA).reinit("string",3)
	e.shortstack	= ""
}

// ddml model structure
struct mStruct {
	string scalar					model			// model; partial, iv, late, etc
	real colvector					id				// id variable (name in Stata will be modelname_id)
	real scalar						nreps			// number of resamplings
	real scalar						ncombos			// number of possible specifications (=0 if not yet estimated)
	real scalar						kfolds			// number of crossfitting folds
	class AssociativeArray scalar	eqnAA			// AssociativeArray with all equations
	class AssociativeArray scalar	estAA			// AA wth all estimation results
	string scalar					nameY			// dependent variable 
	string colvector				nameD			// treatment variable(s)
	string colvector				nameZ			// instrument(s)
	real scalar						ssflag			// flag for shortstacking
	string scalar					strDatavars		// string with expanded names of Stata variables
	real matrix						matDatavars		// matrix with values of Stata variables
	real scalar						crossfitted   	// =1 if crossvalidation has been done; 0 if not
}

struct mStruct init_mStruct()
{

	struct mStruct scalar			d
	
	d.model			= ""
	d.id			= J(0,1,.)
	d.nreps			= 0
	d.ncombos		= 0
	d.kfolds		= 0
	d.nameY			= ""
	d.nameD			= J(1,0,"")
	d.nameZ			= J(1,0,"")
	d.ssflag		= 0
	d.strDatavars	= ""
	d.matDatavars	= J(0,0,.)
	d.crossfitted	= 0
	
	(d.eqnAA).reinit("string",1)
	(d.eqnAA).notfound(NULL)
	
	(d.estAA).reinit("string",2)
	(d.estAA).notfound(NULL)
	
	// initialize counters used for default learner names
	(d.estAA).put(("ycounter","all"),1)
	(d.estAA).put(("dcounter","all"),1)
	(d.estAA).put(("zcounter","all"),1)
	
	return(d)
}

// clear results from all eStructs in mStruct
void clear_model_results(struct mStruct d)
{
	// if not crossfitted, nothing to clear
	if (d.crossfitted==1) {
		struct eStruct scalar			e
		class AssociativeArray scalar	A3
		
		eqnlist = (d.nameY, d.nameD, d.nameZ)
		for (i=1; i<=cols(eqnlist); i++) {
			clear_equation_results((d.eqnAA).get(eqnlist[i]))
		}
		
		// also clear:
		d.ncombos		= 0
		d.crossfitted	= 0
		d.ssflag		= 0
		d.strDatavars	= ""
		d.matDatavars	= J(0,0,.)
		(d.estAA).reinit("string",2)
		(d.estAA).notfound(NULL)
	}
}

// add item about learner to eStruct
void add_learner_item(				struct eStruct e,
									string scalar key1,
									string scalar key2,
									transmorphic scalar s
									)
{
	(e.lrnAA).put((key1,key2),s)
}

// retrieve item about learner from eStruct
transmorphic return_learner_item(	struct eStruct e,
									string scalar key1,
									string scalar key2
									)
{
	return((e.lrnAA).get((key1,key2)))
}

// add result from learner/resample to eStruct
void add_result_item(				struct eStruct e,
									string scalar key1,
									string scalar key2,
									string scalar rep,
									s						//  s can be anything
									)
{
	(e.resAA).put((key1,key2,rep),s)
}

// retrieve result from learner/resample from eStruct
transmorphic return_result_item(	struct eStruct e,
									string scalar key1,
									string scalar key2,
									string scalar rep
									)
{
	return((e.resAA).get((key1,key2,rep)))
}


mata mlib create lddml, dir(PERSONAL) replace
mata mlib add lddml *()
mata mlib index
mata describe using lddml

end

/*
// OLD VERSION //
// some of the string matrices are actually vectors
struct ddmlStruct {
	string scalar		model			// model; partial, iv, late, etc
	real colvector		id				// id variable (name in Stata will be modelname_id)
	real matrix			idFold			// col 1 = id, col 2 = fold identifier (can probably drop this)
	real matrix			idSample		// col 1 = id, col 2 = sample indicator (can probably drop this)
	real scalar			nreps			// number of resamplings
	real scalar			kfolds			// number of crossfitting folds
	string scalar		strDatavars		// string with expanded names of Stata variables
	string scalar		nameY			// dependent variable 
	string colvector	nameYtilde		// names of orthogonalized variables
	string colvector	nameY0tilde		// names of orthogonalized variables
	string colvector	nameY1tilde		// names of orthogonalized variables
	string colvector	nameYopt 		// name of optimal orthog. Y variable
	string colvector	nameY0opt		// name of optimal orthog. Y variable E[Y|D=0]
	string colvector	nameY1opt		// name of optimal orthog. Y variable E[Y|D=1]
	string colvector	nameD			// name of treatment variable(s)
	string matrix		nameDtilde		// names of orthogonalized treatment variables OR name of optimal instrument
	string matrix		nameD0tilde		// names of orthogonalized treatment variables OR name of optimal instrument
	string matrix		nameD1tilde		// names of orthogonalized treatment variables OR name of optimal instrument
	string matrix		nameDopt		// name of optimal orthog. D variable(s) (partial linear model)
	string matrix		nameD0opt		// name of optimal orthog. D variable(s) E[D|Z=0]
	string matrix		nameD1opt		// name of optimal orthog. D variable(s) E[D|Z=1]
	string scalar		nameDH			// treatment variable 
	string colvector	nameDHtilde		// E[D|X,Z]
	string scalar		nameDHopt 		// E[D|X,Z] optimal 
	string colvector	nameZ			// name of instrument(s)
	string matrix		nameZtilde		// names of orthogonalized instruments
	string matrix		nameZopt		// names of optimal orthog. instruments
	pointer rowvector	eqnlist			// list of all orthog equations
	string rowvector	eqnlistNames	// names of corresponding Vtilde variable
	real scalar 		crossfitted   	// =1 if crossvalidation has been done; 0 if not
}
*/

/*
// to add: boolean to indicate min MSE / optimal orthogonalized var
struct eqnStruct {
	string scalar		eqntype			// yeq, deq, zeq or dheq
	real scalar			interactive 	// 
	string scalar		Vname			//  
	string scalar		Vtilde			// also acts as equation name  
	string scalar		Vtilde0			//  
	string scalar		Vtilde1			//  
	string scalar 		Vtilde_h 		// intended for LIE
	real scalar 		resid 			// 
	real matrix			idVtilde		// col 1 = id, col 2 = orthogonalized
	string scalar		eststring
	string scalar 		eststring_h  	// secondary estimation string indended for LIE
	string scalar		command
	string scalar 		command_h 		// (intended for LIE)
	string scalar 		vtype 			// type of variable that is generated by -predict-
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
	// possibly drop, or make a colvector corresponding to the resampling number
	real scalar 		crossfitted   	// =1 if crossvalidation has been done; 0 if not
	real matrix 		stack_weights   // weights from use of pystacked
	real matrix 		stack_weights_h   // weights from use of pystacked
	real matrix 		stack_weights0   // weights from use of pystacked
	real matrix 		stack_weights1   // weights from use of pystacked
}
*/

/*
struct ddmlStruct init_ddmlStruct()
{
	struct ddmlStruct scalar	d

	d.eqnlist		= J(1,0,NULL)
	d.nameY			= ""
	d.nameYtilde	= J(1,0,"")
	d.nameY0tilde	= J(1,0,"")
	d.nameY1tilde	= J(1,0,"")
	d.nameYopt		= ""
	d.nameY0opt		= ""
	d.nameY1opt		= ""
	d.nameD			= J(1,0,"")
	d.nameDtilde	= J(1,0,"")
	d.nameD0tilde	= J(1,0,"")
	d.nameD1tilde	= J(1,0,"")
	d.nameDopt		= J(1,0,"")
	d.nameD0opt		= J(1,0,"")
	d.nameD1opt		= J(1,0,"")
	d.nameDH		= J(1,0,"")
	d.nameDHtilde	= J(1,0,"")
	d.nameDHopt		= J(1,0,"")
	d.nameZ			= J(1,0,"")
	d.nameZtilde	= J(1,0,"")
	d.nameZopt		= J(1,0,"")
	return(d)
}
*/
