*** ddml cross-fitting

* notes:
* why is vtype field needed?
* crossfitting code in separate program
* reporting code in separate subroutine below
* debug reporting in separate subroutine below
* number of resamples set here in reps(.) option
* in eqn struct, (*p).idVtilde is a problem ... but is it needed? currently commented out.
* add noisily option

program _ddml_crossfit_additive, eclass sortpreserve

	syntax [anything] ,								/// 
							[						///
							NOIsily					///
							debug					/// 
							Robust					///
							mname(name)				///
							]

	// no checks included yet

	local debugflag		= "`debug'"~=""
			
	*** extract details of estimation
	
	// model
	mata: st_local("model",`mname'.model)
	mata: st_local("numeqns",strofreal(cols(`mname'.eqnlistNames)))
	mata: st_local("numeqnsY",strofreal(cols(`mname'.nameYtilde)))
	mata: st_local("numeqnsD",strofreal(cols(`mname'.nameDtilde)))
	mata: st_local("numeqnsDH",strofreal(cols(`mname'.nameDHtilde)))
	mata: st_local("numeqnsZ",strofreal(cols(`mname'.nameZtilde)))
	di "Model: `model'"
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("listYtilde",invtokens(`mname'.nameYtilde))
	di "Number of Y estimating equations: `numeqnsY'"
	if `numeqnsD' {
		mata: st_local("listD",invtokens(`mname'.nameD))
		mata: st_local("listDtilde",invtokens(`mname'.nameDtilde))
		di "Number of D estimating equations: `numeqnsD'"
	}
	if `numeqnsDH' {
		mata: st_local("listDH",invtokens(`mname'.nameDH))
		mata: st_local("listDHtilde",invtokens(`mname'.nameDHtilde))
		di "Number of DH estimating equations: `numeqnsDH'"
	}
	if `numeqnsZ' {
		mata: st_local("listZ",invtokens(`mname'.nameZ))
		mata: st_local("listZtilde",invtokens(`mname'.nameZtilde))
		di "Number of Z estimating equations: `numeqnsZ'"
	}
	// initialize opt lists with void matrices of correct dimension
	mata: st_local("numvarsD",strofreal(cols(`mname'.nameD)))
	mata: st_local("numvarsZ",strofreal(cols(`mname'.nameZ)))
	mata: `mname'.nameYopt = J(0,1,"")
	mata: `mname'.nameDopt = J(0,`numvarsD',"")
	mata: `mname'.nameDHopt = J(0,`numvarsD',"")
	mata: `mname'.nameZopt = J(0,`numvarsZ',"")
	// folds and resamplings
	mata: st_local("reps",strofreal(`mname'.nreps))
	mata: st_local("kfolds",strofreal(`mname'.kfolds))

	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eqnStruct()


	if `debugflag' {
		*** report estimates for full sample for debugging purposes
		report_debugging `mname'
	}
	else {
		forvalues m=1/`reps' {
		
			di
			di as text "Starting cross-fitting (sample = `m')"

			*** initialize tilde variables
			forvalues i=1/`numeqns' {
				mata: `eqn'=*(`mname'.eqnlist[1,`i'])
				mata: st_local("vtilde",`eqn'.Vtilde)
				mata: st_local("cvdone",strofreal(`eqn'.crossfitted))
				if (`cvdone'!=1) {
					/// macro m is resampling counter
					cap drop `vtilde'_`m'
					qui gen double `vtilde'_`m'=.
				}
			}

			*** do cross-fitting
			
			forvalues i=1/`numeqns' {
	
				// initialize prior to calling crossfit
				mata: `eqn'=*(`mname'.eqnlist[1,`i'])
				mata: st_local("vtilde",`eqn'.Vtilde)
				mata: st_local("vname",`eqn'.Vname)
				mata: st_local("eststring",`eqn'.eststring)
				mata: st_local("eqntype",`eqn'.eqntype)

				if ("`model'"=="optimaliv") {
					mata: st_local("vtilde_h",`eqn'.Vtilde_h)
					mata: st_local("eststring_h",`eqn'.eststring_h)					
				}

				// seems to be unused
				// mata: st_local("vtype",`eqn'.vtype)
				local touse `mname'_sample

				// request residuals unless optimal IV model & deq or dheq
				if ~("`model'"=="optimaliv"&("`eqntype'"=="deq"|"`eqntype'"=="dheq")) {
					local resid 
					mata: `eqn'.resid = 0
				}
				else {
					local resid resid
					mata: `eqn'.resid = 1
				}
				
				di as text "Cross-fitting equation `i' (`vname', `vtilde')" _c
				
				crossfit if `touse',					///
					eststring(`eststring')				///
					kfolds(`kfolds')					///
					foldvar(`mname'_fid_`m')			/// macro m is resampling counter
					vtilde(`vtilde'_`m')				///
					vtildeh(`vtilde_h'_`m')				/// LIE only
					eststringh(`eststring_h')			/// LIE only
					vname(`vname')						///
					`resid'

				// store MSE and sample size; also set eqn crossfitted flag = 1
				// assumes needed results from crossfit are in r(.) macros
				mata: add_to_eqn(`mname',`i')
				if ("`eqntype'"=="deq"&"`model'"=="optimaliv") {
					mata: add_to_eqn_h(`mname',`i')	
				}	
			}
	
			// for each equation: save names of tilde vars with smallest MSE
			// in each resample m and for each variable in nameY/listD/etc.
			// subroutine will handle case of empty lists
			_ddml_crossfit_update_optlist `mname', etype(yeq) vlist(`nameY') m(`m') model(`model')
			_ddml_crossfit_update_optlist `mname', etype(deq) vlist(`listD') m(`m') model(`model')
			_ddml_crossfit_update_optlist `mname', etype(dheq) vlist(`listDH') m(`m') model(`model')
			_ddml_crossfit_update_optlist `mname', etype(zeq) vlist(`listZ') m(`m') model(`model')
		
		}	// end crossfitting block
	}	// end resampling block
	
	// report results by equation type with resamplings grouped together
	di
	di as res "Reporting crossfitting results:"
	_ddml_crossfit_report `mname'

end

program report_debugging

	syntax name(name=mname)

	// blank eqn - declare this way so that it's a struct and not transmorphic
	tempname eqn
	mata: `eqn' = init_eqnStruct()

	mata: st_local("numeqns",strofreal(cols(`mname'.eqnlistNames)))
	mata: st_local("eqnlist",invtokens(`mname'.eqnlistNames))
	
	*** report estimates for full sample for debugging purposes
	forvalues i=1/`numeqns' {
		mata: `eqn'=*(`mname'.eqnlist[1,`i'])
		mata: st_local("vtilde",`eqn'.Vtilde)
		local do_eqn	: list posof "`vtilde'" in eqnlist

		if `do_eqn' {
			mata: st_local("vname",`eqn'.Vname)
			mata: st_local("eststring",`eqn'.eststring)
			local 0 "`eststring'"
			syntax [anything] , [*]
			local est_main `anything'
			local est_options `options'
			if "`foldlist'"=="" {
				di
				di as res "Estimating equation `i', `vname'/`vtilde':"
				di as res "(full sample, for debugging; no crossfit)"
				// estimate
				`est_main' if `mname'_sample, `est_options'
			}
			else {
				di
				di as res "Estimating equation `i', `vname'/`vtilde':"
				foreach k of numlist `foldlist' {
					di
					di as res "Fold=`k' (for debugging; no crossfit)"
					// estimate
					`est_main' if `mname'_fid!=`k' & `mname'_sample, `est_options'
				}
			}
		}
	}	

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

void add_to_eqn_h(					struct ddmlStruct m,
									real scalar eqnumber)
{
	pointer(struct eqnStruct) scalar p

	cmd 			= st_global("r(cmd_h)")
	mse_h			= st_numscalar("r(mse_h)")
	mse_h_folds		= st_matrix("r(mse_h_folds)")
	n_h				= st_numscalar("r(N_h)")
	n_h_folds		= st_matrix("r(N_h_folds)")
	p				= m.eqnlist[1,eqnumber]
	(*p).MSE_h		= ((*p).MSE_h \ mse_h)
	(*p).N_h		= ((*p).N_h \ n_h)
	(*p).command_h	= cmd

	if (cmd == "pystacked") {
		(*p).stack_weights_h = st_matrix("r(pysw_h)")		 
	}

	// MSE by fold list should be initialized to void 0-by-k matrix
	// (otherwise concat fails because of conformability)
	(*p).MSE_h_folds= ((*p).MSE_h_folds \ mse_h_folds)
	(*p).N_h_folds	= ((*p).N_h_folds \ n_h_folds)

}

// function to set crossfit dummy indicating whether crossfit has been done already
void set_crossfit(					struct ddmlStruct m,
									real scalar eqnumber,
									real scalar cf)
{
	pointer(struct eqnStruct) scalar p
	p				= m.eqnlist[1,eqnumber]
	(*p).crossfitted = cf
}

end
