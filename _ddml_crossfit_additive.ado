*** ddml cross-fitting

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
	
	*** gen folds
	// create foldvar if not empty
	// Stata name will be mname_fid
	cap count if `mname'_fid < .
	if _rc > 0 {
		// fold var does not exist or is not a valid identifier
		cap drop `mname'_fid
		tempvar uni cuni
		qui gen double `uni' = runiform() if `mname'_sample
		qui cumul `uni' if `mname'_sample, gen(`cuni')
		qui gen int `mname'_fid = ceil(`kfolds'*`cuni') if `mname'_sample
		// add fold id to model struct (col 1 = id, col 2 = fold id)
		mata: `mname'.idFold = st_data(., ("`mname'_id", "`mname'_fid"))
	}
	if ("`tabfold'"!="") {
		di
		di "Overview of frequencies by fold:"
		tab `mname'_fid if `mname'_sample
		di
	}
	//

	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eqnStruct()


	if `debugflag' {
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
	}
	else {
		
		*** initialize tilde variables
		forvalues i=1/`numeqns' {
			mata: `eqn'=*(`mname'.eqnlist[1,`i'])
			mata: st_local("vtilde",`eqn'.Vtilde)
			mata: st_local("cvdone",strofreal(`eqn'.crossfitted))
			if (`cvdone'!=1) {
				cap drop `vtilde'
				qui gen double `vtilde'=.
			}
		}

		*** do cross-fitting
		di
		di as text "Cross-fitting equation " _c
		forvalues i=1/`numeqns' {

			if (`i'==`numeqns') {
				di as text "`i'"
			}
			else {
				di as text "`i' " _c
			}

			// has the equation already been crossfitted?
			mata: st_numscalar("cvdone",`eqn'.crossfitted)
			if ("`cvdone'"=="1") continue

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
		}
	
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
	
		// loop through equations, display results, and save names of tilde vars with smallest MSE
		// dep var
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
	}	// end crossfitting block

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
	(*p).idVtilde	= st_data(., tokens(vnames))
	(*p).MSE		= mse
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
