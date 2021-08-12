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
							[ kfolds(integer 5)		///
							NOIsily					///
							debug					/// 
							Robust					///
							TABFold					///
							foldlist(numlist)		///
							mname(name)				///
							/* eqnlist(namelist) */	/// not in use
							reps(integer 0)			///
							]

	// no checks included yet

	local debugflag		= "`debug'"~=""

	/* not in use
	// if eqnlist is empty, populate with full list of eqn names
	if "`eqnlist'"=="" {
		mata: st_local("eqnlist",invtokens(`mname'.eqnlistNames))
	}
	*/
			
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
	
	// folds and fold IDs
	mata: st_local("hasfoldvars",strofreal(cols(`mname'.idFold)))
	// if reps is specified then overwrite any existing fold vars.
	if `reps'>0 {
		local hasfoldvars = 0
	}
	else {
		// default reps = 1
		local reps = 1
	}
	// store number of resamplings
	mata: `mname'.nreps = `reps'

	// if empty:
	// add fold IDs to model struct (col 1 = id, col 2 = fold id 1, col 3 = fold id 2 etc.)
	// first initialize with id
	
	if `hasfoldvars'==0 {
		mata: `mname'.idFold = st_data(., ("`mname'_id"))
	}

	forvalues m=1/`reps' {
	
		if `hasfoldvars'==0 {
			*** gen folds
			cap drop `mname'_fid_`m'
			tempvar uni cuni
			qui gen double `uni' = runiform() if `mname'_sample
			qui cumul `uni' if `mname'_sample, gen(`cuni')
			qui gen int `mname'_fid_`m' = ceil(`kfolds'*`cuni') if `mname'_sample
			// add fold id to model struct (col 1 = id, col 2 = fold id)
			mata: `mname'.idFold = (`mname'.idFold , st_data(., ("`mname'_fid_`m'")))
		}
		if ("`tabfold'"!="") {
			di
			di "Overview of frequencies by fold (sample `m'):"
			tab `mname'_fid_`m' if `mname'_sample
			di
		}
	
	}

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
	
				// has the equation already been crossfitted?
				mata: st_numscalar("cvdone",`eqn'.crossfitted)
				if ("`cvdone'"=="1") continue
	
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
				
				if ~("`model'"=="optimaliv"&("`eqntype'"=="deq"|"`eqntype'"=="dheq")) {
					// request residuals unless optimal IV model & deq or dheq
					local resid resid
				}
				else {
					// default is predicted values
					local resid
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
				
				// store MSE and sample size
				mata: add_to_eqn(`mname',`i',"`mname'_id `vtilde'", `r(mse)',`r(N)',"`r(cmd)'")
				if ("`eqntype'"=="deq"&"`model'"=="optimaliv") {
					mata: add_to_eqn_h(`mname',`i',"`mname'_id `vtilde'", `r(mse_h)',`r(N_h)',"`r(cmd_h)'")	
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

	// set crossfitted field to 1	
	mata: `mname'.crossfitted = 1
	
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
									real scalar eqnumber,
									string scalar vnames,
									real scalar mse,
									real scalar n,
									string scalar cmd)
{
	pointer(struct eqnStruct) scalar p
	p				= m.eqnlist[1,eqnumber]
	//(*p).idVtilde	= st_data(., tokens(vnames))
	(*p).MSE		= ((*p).MSE \ mse)
	(*p).N			= n
	(*p).command	= cmd
}

void add_to_eqn_h(					struct ddmlStruct m,
									real scalar eqnumber,
									string scalar vnames,
									real scalar mse,
									real scalar n,
									string scalar cmd_h)
{
	pointer(struct eqnStruct) scalar p
	p				= m.eqnlist[1,eqnumber]
	//(*p).idVtilde	= st_data(., tokens(vnames))
	(*p).MSE_h		= ((*p).MSE_h \ mse)
	(*p).N_h		= n
	(*p).command_h	= cmd_h
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
