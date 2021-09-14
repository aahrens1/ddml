* last edited: 18 jun 2021

* notes:
* why is crossfitted field set in additive code but not here?
* check it's correct that interactive-type estimation always goes with reporting mse0 and mse1
* are multiple Ds allowed?

*** ddml cross-fitting for the interactive model & LATE
program _ddml_crossfit_interactive, eclass sortpreserve

	syntax [anything] [if] [in] , /// 
							[ kfolds(integer 5)		///
							NOIsily					///
							debug					/// 
							Robust					///
							TABFold					///
							foldlist(numlist)		///
							mname(name)				///
							reps(integer 0)			/// if reps specified then overwrite any existing fold vars
							]

	// no checks included yet
	// no marksample yet

	local debugflag		= "`debug'"~=""
	if ("`noisily'"=="") local qui qui
		
	*** extract details of estimation
	
	// model
	mata: st_local("model",`mname'.model)
	di "Model: `model'"
	mata: st_local("numeqns",strofreal(cols(`mname'.eqnlistNames)))
	mata: st_local("numeqnsY",strofreal(cols(`mname'.nameYtilde)))
	mata: st_local("numeqnsD",strofreal(cols(`mname'.nameDtilde)))
	mata: st_local("numeqnsZ",strofreal(cols(`mname'.nameZtilde)))
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("listYtilde",invtokens(`mname'.nameYtilde))
	di "Number of Y estimating equations: `numeqnsY'"
	if `numeqnsD' {
		mata: st_local("listD",invtokens(`mname'.nameD))
		mata: st_local("listDtilde",invtokens(`mname'.nameDtilde))
		di "Number of D estimating equations: `numeqnsD'"
	}
	if `numeqnsZ' {
		mata: st_local("listZ",invtokens(`mname'.nameZ))
		mata: st_local("listZtilde",invtokens(`mname'.nameZtilde))
		di "Number of Z estimating equations: `numeqnsZ'"
	}		 
	// initialize opt lists with void matrices of correct dimension
	mata: st_local("numvarsD",strofreal(cols(`mname'.nameD)))
	mata: st_local("numvarsZ",strofreal(cols(`mname'.nameZ)))
	mata: `mname'.nameY0opt = J(0,1,"")
	mata: `mname'.nameY1opt = J(0,1,"")
	mata: `mname'.nameDopt = J(0,`numvarsD',"")
	mata: `mname'.nameD0opt = J(0,`numvarsD',"")
	mata: `mname'.nameD1opt = J(0,`numvarsD',"")
	mata: `mname'.nameZopt = J(0,`numvarsZ',"")
	
	// folds and resamplings
	mata: st_local("reps",strofreal(`mname'.nreps))
	mata: st_local("kfolds",strofreal(`mname'.kfolds))

	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eqnStruct()

	forvalues m=1/`reps' {
	
		di
		di as text "Starting cross-fitting (sample = `m')"

		*** initialize tilde variables
		forvalues i=1/`numeqns' {
			mata: `eqn'=*(`mname'.eqnlist[1,`i'])
			mata: st_local("vtilde",`eqn'.Vtilde)
			/// macro m is resampling counter
			cap drop `vtilde'_`m'
			qui gen double `vtilde'_`m'=.
		}
			
		*** do cross-fitting
		
		forvalues i=1/`numeqns' {
		
			// initialize prior to calling crossfit
			mata: `eqn'=*(`mname'.eqnlist[1,`i'])
			mata: st_local("vtilde",`eqn'.Vtilde)
			mata: st_local("vname",`eqn'.Vname)
			mata: st_local("eststring",`eqn'.eststring)
			mata: st_local("eqntype",`eqn'.eqntype)
			// seems to be unused
			// mata: st_local("vtype",`eqn'.vtype)
			local touse `mname'_sample
			// always request fitted values
			local resid
			mata: `eqn'.resid = 0

			di as text "Cross-fitting equation `i' (`vname', `vtilde')" _c
	
			/* why not used here but used in interactive code?
			// has the equation already been crossfitted?
			mata: st_numscalar("cvdone",`eqn'.crossfitted)
			if ("`cvdone'"=="1") continue
			*/
			
			if ("`eqntype'"=="yeq") & ("`model'"=="interactive") {
				local treatvar	`listD'
			}
			else if ("`eqntype'"=="yeq") & ("`model'"=="late") {
				local treatvar	`listZ'
			}
			else if ("`eqntype'"=="deq") & ("`model'"=="interactive") {
				local treatvar
			}
			else if ("`eqntype'"=="deq") & ("`model'"=="late") {
				local treatvar	`listZ'
			}
			else if ("`eqntype'"=="zeq") {
				local treatvar
			}
			else {
				di as err "Unknown equation type `eqntype'"
				exit 198
			}
	
			crossfit if `touse',					///
				eststring(`eststring')				///
				kfolds(`kfolds')					///
				foldvar(`mname'_fid_`m')			///
				vtilde(`vtilde'_`m')				///
				vname(`vname')						///
				treatvar(`treatvar')				///
				`resid'
			
			// store MSE and sample size
			// assumes needed results from crossfit are in r(.) macros
			if ("`eqntype'"=="yeq") {
				mata: add_to_eqn01(`mname',`i', "0")
				mata: add_to_eqn01(`mname',`i', "1")
			}
			else if ("`eqntype'"=="deq") & ("`model'"=="interactive") {
				mata: add_to_eqn(`mname',`i')
			}
			else if ("`eqntype'"=="deq") {
				mata: add_to_eqn01(`mname',`i', "0")
				mata: add_to_eqn01(`mname',`i', "1")
			}
			else if ("`eqntype'"=="zeq") {
				mata: add_to_eqn(`mname',`i')
			}
	
		}

		// for each equation: save names of tilde vars with smallest MSE
		// in each resample m and for each variable in nameY/listD/etc.
		// subroutine will handle case of empty lists
		if "`model'"=="interactive" {
			// dependent variable
			_ddml_crossfit_update_optlist `mname', etype(yeq) vlist(`nameY') m(`m') zett(0)
			_ddml_crossfit_update_optlist `mname', etype(yeq) vlist(`nameY') m(`m') zett(1)
			// D variable
			_ddml_crossfit_update_optlist `mname', etype(deq) vlist(`listD') m(`m')
		}
		// late model
		else if "`model'"=="late" {
			// dependent variable
			_ddml_crossfit_update_optlist `mname', etype(yeq) vlist(`nameY') m(`m') zett(0)
			_ddml_crossfit_update_optlist `mname', etype(yeq) vlist(`nameY') m(`m') zett(1)
			// D variable
			_ddml_crossfit_update_optlist `mname', etype(deq) vlist(`listD') m(`m') zett(0)
			_ddml_crossfit_update_optlist `mname', etype(deq) vlist(`listD') m(`m') zett(1)
			// Z variable
			_ddml_crossfit_update_optlist `mname', etype(zeq) vlist(`listZ') m(`m')
		}
		else {
			di as err "internal error - unknown model `model'"
			exit 198
		}
	}

	// set crossfitted field to 1	
	mata: `mname'.crossfitted = 1

	// report results by equation type with resamplings grouped together
	di
	di as res "Reporting crossfitting results:"
	_ddml_crossfit_report `mname'

end


mata:

struct eqnStruct init_eqnStruct()
{
	struct eqnStruct scalar		e
	return(e)
}

void add_to_eqn(					struct ddmlStruct m,
									real scalar eqnumber)

{
	pointer(struct eqnStruct) scalar p

	cmd 			= st_global("r(cmd)")
	mse				= st_numscalar("r(mse)")
	mse_folds		= st_matrix("r(mse_folds)")
	n				= st_numscalar("r(N)")
	n_folds			= st_matrix("r(N_folds)")
	p				= m.eqnlist[1,eqnumber]
	(*p).MSE		= ((*p).MSE \ mse)
	(*p).N			= ((*p).N \ n)
	(*p).command	= cmd

	if (cmd == "pystacked") {
		(*p).stack_weights = st_matrix("r(pysw)")		 
	}

	// MSE by fold list should be initialized to void 0-by-k matrix
	// (otherwise concat fails because of conformability)
	(*p).MSE_folds	= ((*p).MSE_folds \ mse_folds)
	(*p).N_folds	= ((*p).N_folds \ n_folds)
	
	// set crossfitted flag = 1
	(*p).crossfitted	= 1

}

void add_to_eqn01(					struct ddmlStruct m,
									real scalar eqnumber,
									string scalar Z)
{
	pointer(struct eqnStruct) scalar p

	cmd 			= st_global("r(cmd)")
	mse				= st_numscalar("r(mse" + Z + ")")
	mse_folds		= st_matrix("r(mse" + Z + "_folds)")
	n				= st_numscalar("r(N" + Z + ")")
	n_folds			= st_matrix("r(N" + Z + "_folds)")

	p				= m.eqnlist[1,eqnumber]
	if (Z=="0") {
		(*p).MSE0		= ((*p).MSE0 \ mse)
		(*p).N0			= ((*p).N0 \ n)
		// MSE by fold list should be initialized to void 0-by-k matrix
		// (otherwise concat fails because of conformability)
		(*p).MSE0_folds	= ((*p).MSE0_folds \ mse_folds)
		(*p).N0_folds	= ((*p).N0_folds \ n_folds)
		if (cmd == "pystacked") {
			(*p).stack_weights0 = st_matrix("r(pysw0)")		 
		}
	}
	else {
		(*p).MSE1		= ((*p).MSE1 \ mse)
		(*p).N1			= ((*p).N1 \ n)
		// MSE by fold list should be initialized to void 0-by-k matrix
		// (otherwise concat fails because of conformability)
		(*p).MSE1_folds	= ((*p).MSE1_folds \ mse_folds)
		(*p).N1_folds	= ((*p).N1_folds \ n_folds)
		if (cmd == "pystacked") {
			(*p).stack_weights1 = st_matrix("r(pysw1)")		 
		}
	}
	
	// set crossfitted flag = 1
	(*p).crossfitted	= 1

}

end
