*** ddml cross-fitting

* notes:
* why is vtype field needed?
* crossfitting code in separate program
* reporting code in separate subroutine below
* debug reporting in separate subroutine below
* number of resamples set here in reps(.) option
* in eqn struct, (*p).idVtilde is a problem ... but is it needed? currently commented out.

program _ddml_crossfit_additive, eclass sortpreserve

	syntax [anything] ,								/// 
							[ kfolds(integer 2)		///
							NOIsily					///
							debug					/// 
							Robust					///
							TABFold					///
							mname(name)				///
							eqnlist(namelist)		///
							foldlist(numlist)		///
							reps(integer 1)			///
							]

	// no checks included yet

	local debugflag		= "`debug'"~=""

	// if eqnlist is empty, populate with full list of eqn names
	if "`eqnlist'"=="" {
		mata: st_local("eqnlist",invtokens(`mname'.eqnlistNames))
	}
			
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
	
	// folds and fold IDs
	mata: st_local("hasfoldvars",strofreal(cols(`mname'.idFold)))

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
	
				di as text "Cross-fitting equation `i'" _c
	
				// has the equation already been crossfitted?
				mata: st_numscalar("cvdone",`eqn'.crossfitted)
				if ("`cvdone'"=="1") continue
				
				// initialize prior to calling crossfit
				mata: `eqn'=*(`mname'.eqnlist[1,`i'])
				mata: st_local("vtilde",`eqn'.Vtilde)
				mata: st_local("vname",`eqn'.Vname)
				mata: st_local("eststring",`eqn'.eststring)
				mata: st_local("eqntype",`eqn'.eqntype)
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
				
				crossfit if `touse',					///
					eststring(`eststring')				///
					kfolds(`kfolds')					///
					foldvar(`mname'_fid_`m')			/// macro m is resampling counter
					vtilde(`vtilde'_`m')				///
					vname(`vname')						///
					`resid'
				
				// store MSE and sample size
				mata: add_to_eqn(`mname',`i',"`mname'_id `vtilde'", `r(mse)',`r(N)')
				
				/*
				*** old code replaced by crossfit subcommand
	
				forvalues k = 1(1)`kfolds' {
				
					// ML is applied to I^c sample (all data ex partition k)
					qui {
						
							mata: `eqn'=*(`mname'.eqnlist[1,`i'])
							mata: st_local("vtilde",`eqn'.Vtilde)
							mata: st_local("vname",`eqn'.Vname)
							mata: st_local("eststring",`eqn'.eststring)
							mata: st_local("eqntype",`eqn'.eqntype)
							mata: st_local("vtype",`eqn'.vtype)
							local 0 "`eststring'"
							syntax [anything] , [*]
							local est_main `anything'
							local est_options `options'
							di as res "Estimating equation `i':"
							di as res "  est_main: `est_main'"
							di as res "  est_options: `est_options'"
							
							tempvar vtilde_i
	
							// estimate excluding kth fold
							`est_main' if `mname'_fid!=`k' & `mname'_sample, `est_options'
							// get fitted values and residuals for kth fold	
							qui predict `vtype' `vtilde_i' if `mname'_fid==`k' & `mname'_sample
	
	
							if ("`model'"=="optimaliv"&("`eqntype'"=="deq"|"`eqntype'"=="dheq")) {
								// get predicted values if optimal IV model & deq or dheq
								qui replace `vtilde' = `vtilde_i' if `mname'_fid==`k' & `mname'_sample
							} 
							else {
								// get residuals
								qui replace `vtilde' = `vname' - `vtilde_i' if `mname'_fid==`k' & `mname'_sample
							}
						}
					}
					mata: set_crossfit(`mname',`i',1) // indicate that cross-validation has been done
				*/				
					
			}
	
		
			/*
			*** old code replaced by crossfit subcommand
			*** calculate MSE, store orthogonalized variables, etc.
			forvalues i=1/`numeqns' {
				mata: `eqn'=*(`mname'.eqnlist[1,`i'])
				mata: st_local("vtilde",`eqn'.Vtilde)
				mata: st_local("vname",`eqn'.Vname)
				tempvar vtilde_sq
				if ("`model'"=="optimaliv"&("`eqntype'"=="deq"|"`eqntype'"=="dheq")) {
					// need to use residuals here
					qui gen double `vtilde_sq' = (`vname'-`vtilde')^2 if `mname'_sample
					qui sum `vtilde_sq' if `mname'_sample, meanonly
					mata: add_to_eqn(`mname',`i',"`mname'_id `vtilde'", `r(mean)',`r(N)')
				} 
				else {
					qui gen double `vtilde_sq' = `vtilde'^2 if `mname'_sample
					qui sum `vtilde_sq' if `mname'_sample, meanonly
					mata: add_to_eqn(`mname',`i',"`mname'_id `vtilde'", `r(mean)',`r(N)')
				}
			}
			*/
				
			// loop through equations, display results, and save names of tilde vars with smallest MSE
			// subroutine will handle case of empty lists
			// note that vtype must match exactly
			// note that subroutine uses field names of struct
			di
			di as res "Reporting crossfitting results (sample=`m')
			report_crossfit_result `mname', vtype(y|X) vlist(`nameY') m(`m')
			report_crossfit_result `mname', vtype(D|X) vlist(`listD') m(`m')
			report_crossfit_result `mname', vtype(D|X,Z) vlist(`listDH') m(`m')
			report_crossfit_result `mname', vtype(Z|X) vlist(`listZ') m(`m')
			
			// dep var
			/*
			di
			di as res "Mean-squared error for y|X:"
			di _col(2) "Name" _c
			di _col(20) "Orthogonalized" _c
			di _col(40) "Command" _c
			di _col(54) "N" _c
			di _col(65) "MSPE"
			di "{hline 75}"
			_ddml_display_mspe `mname', vname(`nameY')
			mata: `mname'.nameYopt		= "`r(optname)'"
		
			// loop through D vars (if any)
			if `numeqnsD' {
				di
				di as res "Mean-squared error for D|X:"
				di _col(2) "Name" _c
				di _col(20) "Orthogonalized" _c
				di _col(40) "Command" _c
				di _col(54) "N" _c
				di _col(65) "MSPE"
				di "{hline 75}"
				// clear opt list
				mata: `mname'.nameDopt = J(1,0,"") 
				foreach var of varlist `listD' {
					_ddml_display_mspe `mname', vname(`var')
					mata: `mname'.nameDopt = (`mname'.nameDopt, "`r(optname)'")
				}
			}
			
			// loop through DH vars (if any)
			if `numeqnsDH' {
				di
				di as res "Mean-squared error for D|X,Z:"
				di _col(2) "Name" _c
				di _col(20) "Orthogonalized" _c
				di _col(40) "Command" _c
				di _col(54) "N" _c
				di _col(65) "MSPE"
				di "{hline 75}"
				// clear opt list
				mata: `mname'.nameDHopt = J(1,0,"") 
				foreach var of varlist `listDH' {
					_ddml_display_mspe `mname', vname(`var')
					mata: `mname'.nameDHopt = (`mname'.nameDHopt, "`r(optname)'")
				}
			}
	
			// loop through Z vars (if any)
			if `numeqnsZ' {
				di
				di as res "Mean-squared error for Z|X:"
				di _col(2) "Name" _c
				di _col(20) "Orthogonalized" _c
				di _col(40) "Command" _c
				di _col(54) "N" _c
				di _col(65) "MSPE"
				di "{hline 75}"
				// clear opt list
				mata: `mname'.nameZopt = J(1,0,"") 
				foreach var of varlist `listZ' {
					_ddml_display_mspe `mname', vname(`var')
					mata: `mname'.nameZopt = (`mname'.nameZopt, "`r(optname)'")
				}
			}
			*/
			
		}	// end crossfitting block
	}	// end resampling block

end

program report_crossfit_result
	syntax name(name=mname), vtype(string) [ vlist(string) m(integer 1) ]

	// set struct field name
	if "`vtype'"=="y|X" {
		local optname nameYopt
	}
	else if "`vtype'"=="D|X" {
		local optname nameDopt
	}
	else if "`vtype'"=="D|X,Z" {
		local optname nameDHopt
	}
	else if "`vtype'"=="Z|X" {
		local optname nameZopt
	}
	
	// may be called with empty list (e.g. if no endog regressors)
	local numeqns	: word count `vlist'
	if `numeqns' > 0 {
		
		di
		di as res "Mean-squared error for `vtype':"
		di _col(2) "Name" _c
		di _col(20) "Orthogonalized" _c
		di _col(40) "Command" _c
		di _col(54) "N" _c
		di _col(65) "MSPE"
		di "{hline 75}"
		// clear opt list
		// mata: `mname'.`optname' = J(1,0,"")
		foreach var of varlist `vlist' {
			// m is the rep number
			_ddml_display_mspe `mname', vname(`var') m(`m')
			local optlist `optlist' `r(optname)'
		}
		if `m'==1 {
			mata: `mname'.`optname' = tokens("`optlist'")'
		}
		else {
			mata: `mname'.`optname' = (`mname'.`optname' \ tokens("`optlist'")')
		}

	}
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
									real scalar n)
{
	pointer(struct eqnStruct) scalar p
	p				= m.eqnlist[1,eqnumber]
	//(*p).idVtilde	= st_data(., tokens(vnames))
	(*p).MSE		= ((*p).MSE \ mse)
	(*p).N			= n
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
