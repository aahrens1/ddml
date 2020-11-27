*** ddml cross-fitting
program _ddml_crossfit_partial, eclass sortpreserve

	syntax [anything] [if] [in] , /// 
							[ kfolds(integer 2) ///
							NOIsily ///
							debug /// 
							Robust ///
							TABFold ///
							yrclass ///
							drclass /// 
							mname(name)	///
							]

	// no checks included yet
	// no marksample yet

	local debugflag		= "`debug'"~=""
		
	*** extract details of estimation
	
	// model
	mata: st_local("model",`mname'.model)
	mata: st_local("numeqnsY",strofreal(cols(`mname'.eqnlistY)))
	mata: st_local("numeqnsD",strofreal(cols(`mname'.eqnlistD)))
	mata: st_local("numeqnsZ",strofreal(cols(`mname'.eqnlistZ)))
	di "Model: `model'"
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
	
	*** gen folds
	// create foldvar if not empty
	// Stata name will be mname_fid
	cap count if `mname'_fid < .
	if _rc > 0 {
		// fold var does not exist or is not a valid identifier
		cap drop `mname'_fid
		tempvar uni cuni
		qui gen double `uni' = runiform()
		qui cumul `uni', gen(`cuni')
		qui gen int `mname'_fid = ceil(`kfolds'*`cuni')
	}
	// add fold id to model struct (col 1 = id, col 2 = fold id)
	mata: `mname'.idFold = st_data(., ("`mname'_id", "`mname'_fid"))
	if ("`tabfold'"!="") {
		di
		di "Overview of frequencies by fold:"
		tab `mname'_fid
		di
	}
	//

	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eqnStruct()

	*** initialize tilde variables
	// Y equations
	forvalues i=1/`numeqnsY' {
		mata: `eqn'=*(`mname'.eqnlistY[1,`i'])
		mata: st_local("vtilde",`eqn'.Vtilde)
		cap drop `mname'_`vtilde'
		qui gen double `mname'_`vtilde'=.
	}
	// D equations
	forvalues i=1/`numeqnsD' {
		mata: `eqn'=*(`mname'.eqnlistD[1,`i'])
		mata: st_local("vtilde",`eqn'.Vtilde)
		cap drop `mname'_`vtilde'
		qui gen double `mname'_`vtilde'=.
	}
	// Z equations
	if ("`model'"=="iv") {
		forvalues i=1/`numeqnsZ' {
			mata: `eqn'=*(`mname'.eqnlistZ[1,`i'])
			mata: st_local("vtilde",`eqn'.Vtilde)
			cap drop `mname'_`vtilde'
			qui gen double `mname'_`vtilde'=.
		}
	}

	*** estimate equations that do not require crossfitting
	// also report estimates for full sample for debugging purposes
	// Y equations
	tempname crossfit
	forvalues i=1/`numeqnsY' {
		mata: `eqn'=*(`mname'.eqnlistY[1,`i'])
		mata: st_numscalar("`crossfit'",`eqn'.crossfit)
		if `crossfit'==0 | `debugflag' {
			mata: st_local("vtilde",`eqn'.Vtilde)
			mata: st_local("vname",`eqn'.Vname)
			mata: st_local("eststring",`eqn'.eststring)
			local 0 "`eststring'"
			syntax [anything] , [*]
			local est_main `anything'
			local est_options `options'
			if `debugflag' {
				if `crossfit'==0 {
					di
					di "Y estimating equation `i' (full sample, for debugging; no crossfit):"
				}
				else {
					di
					di "Y estimating equation `i' (no crossfit):"
				}
				di "  est_main: `est_main'"
				di "  est_options: `est_options'"
			}
			else {
				// set quietly flag
				local yquietly quietly
			}
			// estimate
			`yquietly' `est_main', `est_options'
			// get fitted values and residuals for no-crossfit case
			if `crossfit'==0 {
				tempvar vtilde_i
				qui predict double `vtilde_i'
				qui replace `mname'_`vtilde' = `vname' - `vtilde_i'
			}
		}
	}
	forvalues i=1/`numeqnsD' {
		mata: `eqn'=*(`mname'.eqnlistD[1,`i'])
		mata: st_numscalar("`crossfit'",`eqn'.crossfit)
		if `crossfit'==0 | `debugflag' {
			mata: st_local("vtilde",`eqn'.Vtilde)
			mata: st_local("vname",`eqn'.Vname)
			mata: st_local("eststring",`eqn'.eststring)
			local 0 "`eststring'"
			syntax [anything] , [*]
			local est_main `anything'
			local est_options `options'
			if `debugflag' {
				if `crossfit'==0 {
					di
					di "D estimating equation `i' (full sample, for debugging; no crossfit):"
				}
				else {
					di
					di "D estimating equation `i' (no crossfit):"
				}
				di "  est_main: `est_main'"
				di "  est_options: `est_options'"
			}
			else {
				// set quietly flag
				local dquietly quietly
			}
			// estimate
			`dquietly' `est_main', `est_options'
			// get fitted values and residuals for no-crossfit case
			if `crossfit'==0 {
				tempvar vtilde_i
				qui predict double `vtilde_i'
				qui replace `mname'_`vtilde' = `vname' - `vtilde_i'
			}
		}
	}
	forvalues i=1/`numeqnsZ' {
		mata: `eqn'=*(`mname'.eqnlistZ[1,`i'])
		mata: st_numscalar("r(crossfit)",`eqn'.crossfit)
		if r(crossfit)==0 | `debugflag' {
			mata: st_local("vtilde",`eqn'.Vtilde)
			mata: st_local("vname",`eqn'.Vname)
			mata: st_local("eststring",`eqn'.eststring)
			local 0 "`eststring'"
			syntax [anything] , [*]
			local est_main `anything'
			local est_options `options'
			if `debugflag' {
				if `crossfit'==0 {
					di
					di "Z estimating equation `i' (full sample, for debugging; no crossfit):"
				}
				else {
					di
					di "Z estimating equation `i' (no crossfit):"
				}
				di "  est_main: `est_main'"
				di "  est_options: `est_options'"
			}
			else {
				// set quietly flag
				local zquietly quietly
			}
			// estimate
			`zquietly' `est_main', `est_options'
			// get fitted values and residuals
			if `crossfit'==0 {
				tempvar vtilde_i
				qui predict double `vtilde_i'
				qui replace `mname'_`vtilde' = `vname' - `vtilde_i'
			}
		}
	}

	
	*** do cross-fitting
	di
	di as text "Cross-fitting fold " _c
	forvalues k = 1(1)`kfolds' {
	
		if (`k'==`kfolds') {
			di as text "`k'"
		}
		else {
			di as text "`k' " _c
		}
		// ML is applied to I^c sample (all data ex partition k)
		qui {
		
			// Y equations
			forvalues i=1/`numeqnsY' {
				mata: `eqn'=*(`mname'.eqnlistY[1,`i'])
				mata: st_numscalar("r(crossfit)",`eqn'.crossfit)
				if r(crossfit) {
					mata: st_local("vtilde",`eqn'.Vtilde)
					mata: st_local("vname",`eqn'.Vname)
					mata: st_local("eststring",`eqn'.eststring)
					local 0 "`eststring'"
					syntax [anything] , [*]
					local est_main `anything'
					local est_options `options'
					di "Y estimating equation `i':"
					di "  est_main: `est_main'"
					di "  est_options: `est_options'"
	
					// estimate excluding kth fold
					`est_main' if `mname'_fid!=`k', `est_options'
					// get fitted values and residuals for kth fold	
					tempvar vtilde_i
					qui predict double `vtilde_i' if `mname'_fid==`k' 
					qui replace `mname'_`vtilde' = `vname' - `vtilde_i' if `mname'_fid==`k'
				}
			}
			// D equations
			forvalues i=1/`numeqnsD' {
				mata: `eqn'=*(`mname'.eqnlistD[1,`i'])
				mata: st_numscalar("r(crossfit)",`eqn'.crossfit)
				if r(crossfit) {
					mata: st_local("vtilde",`eqn'.Vtilde)
					mata: st_local("vname",`eqn'.Vname)
					mata: st_local("eststring",`eqn'.eststring)
					local 0 "`eststring'"
					syntax [anything] , [*]
					local est_main `anything'
					local est_options `options'
					di "D estimating equation `i':"
					di "  est_main: `est_main'"
					di "  est_options: `est_options'"
	
					// estimate excluding kth fold
					`est_main' if `mname'_fid!=`k', `est_options'
					// get fitted values and residuals for kth fold	
					tempvar vtilde_i
					qui predict double `vtilde_i' if `mname'_fid==`k' 
					qui replace `mname'_`vtilde' = `vname' - `vtilde_i' if `mname'_fid==`k'
				}
			}
			// Z equations
			forvalues i=1/`numeqnsZ' {
				mata: `eqn'=*(`mname'.eqnlistZ[1,`i'])
				mata: st_numscalar("r(crossfit)",`eqn'.crossfit)
				if r(crossfit) {
					mata: st_local("vtilde",`eqn'.Vtilde)
					mata: st_local("vname",`eqn'.Vname)
					mata: st_local("eststring",`eqn'.eststring)
					local 0 "`eststring'"
					syntax [anything] , [*]
					local est_main `anything'
					local est_options `options'
					di "Z estimating equation `i':"
					di "  est_main: `est_main'"
					di "  est_options: `est_options'"
	
					// estimate excluding kth fold
					`est_main' if `mname'_fid!=`k', `est_options'
					// get fitted values and residuals for kth fold	
					tempvar vtilde_i
					qui predict double `vtilde_i' if `mname'_fid==`k' 
					qui replace `mname'_`vtilde' = `vname' - `vtilde_i' if `mname'_fid==`k'
				}
			}
		}
	}	

	*** calculate MSE, store orthogonalized variables, etc.
	forvalues i=1/`numeqnsY' {
		mata: `eqn'=*(`mname'.eqnlistY[1,`i'])
		mata: st_local("vtilde",`eqn'.Vtilde)
		tempvar vtilde_sq
		qui gen double `vtilde_sq' = `mname'_`vtilde'^2
		qui sum `vtilde_sq', meanonly
		mata: add_to_Yeqn(`mname',`i',"`mname'_id `mname'_`vtilde'", `r(mean)',`r(N)')
	}
	forvalues i=1/`numeqnsD' {
		mata: `eqn'=*(`mname'.eqnlistD[1,`i'])
		mata: st_local("vtilde",`eqn'.Vtilde)
		tempvar vtilde_sq
		qui gen double `vtilde_sq' = `mname'_`vtilde'^2
		qui sum `vtilde_sq', meanonly
		mata: add_to_Deqn(`mname',`i',"`mname'_id `mname'_`vtilde'", `r(mean)',`r(N)')
	}
	if ("`model'"=="iv") {
		forvalues i=1/`numeqnsZ' {
			mata: `eqn'=*(`mname'.eqnlistZ[1,`i'])
			mata: st_local("vtilde",`eqn'.Vtilde)
			tempvar vtilde_sq
			qui gen double `vtilde_sq' = `mname'_`vtilde'^2
			qui sum `vtilde_sq', meanonly
			mata: add_to_Zeqn(`mname',`i',"`mname'_id `mname'_`vtilde'", `r(mean)',`r(N)')
		}
	}
	
	// loop through equations and pick out those with smallest MSE

	// dep var
	di
	di as res "Mean-squared error for y|X:"
	di _col(2) "Name" _c
	di _col(20) "Orthogonalized" _c
	di _col(40) "Command" _c
	di _col(54) "N" _c
	di _col(65) "MSPE"
	di "{hline 75}"
	// initialize
	local yminmse = .
	forvalues i=1/`numeqnsY' {
		mata: `eqn'=*(`mname'.eqnlistY[1,`i'])
		mata: st_local("vname",`eqn'.Vname)
		mata: st_local("vtilde",`eqn'.Vtilde)
		mata: st_local("command",`eqn'.command)
		mata: st_numscalar("r(MSE)",`eqn'.MSE)
		mata: st_numscalar("r(N)",`eqn'.N)
		local MSE = `r(MSE)'
		local N = `r(N)'
		if `MSE' < `yminmse' {
			local yopt `vtilde'
			local yminmse `MSE'
		}
		di _col(2) "`vname'" _c
		di _col(20) "`vtilde'" _c
		di _col(40) "`command'" _c
		di _col(50) %6.0f `N' _c
		di _col(60) %10.6f `MSE'
	}
	mata: `mname'.nameYopt		= "`yopt'"

	// D vars
	di
	di as res "Mean-squared error for D|X:"
	di _col(2) "Name" _c
	di _col(20) "Orthogonalized" _c
	di _col(40) "Command" _c
	di _col(54) "N" _c
	di _col(65) "MSPE"
	di "{hline 75}"
	// initialize
	local first_dopt	=1
	foreach var of varlist `listD' {
		// initialize
		local dminmse = .
		forvalues i=1/`numeqnsD' {
			mata: `eqn'=*(`mname'.eqnlistD[1,`i'])
			mata: st_local("vname",`eqn'.Vname)
			mata: st_local("vtilde",`eqn'.Vtilde)
			mata: st_local("command",`eqn'.command)
			mata: st_numscalar("r(MSE)",`eqn'.MSE)
			mata: st_numscalar("r(N)",`eqn'.N)
			local MSE = `r(MSE)'
			local N = `r(N)'
			if "`var'"=="`vname'" {
				if `MSE' < `dminmse' {
					local dopt `vtilde'
					local dminmse `MSE'
				}
				di _col(2) "`vname'" _c
				di _col(20) "`vtilde'" _c
				di _col(40) "`command'" _c
				di _col(50) %6.0f `N' _c
				di _col(60) %10.6f `MSE'
			}
		}
		// if nameDopt already has vname in it, then add it to the list
		if `first_dopt' {
			mata: `mname'.nameDopt	= "`dopt'"
			local first_dopt	=0
		}
		else {
			mata: `mname'.nameDopt	= (`mname'.nameDopt, "`dopt'")
		}

	}

	// Z vars
	if ("`model'"=="iv") {
		di
		di as res "Mean-squared error for Z|X:"
		di _col(2) "Name" _c
		di _col(20) "Orthogonalized" _c
		di _col(40) "Command" _c
		di _col(54) "N" _c
		di _col(65) "MSPE"
		di "{hline 75}"
		// initialize
		local first_zopt	=1
		foreach var of varlist `listZ' {
			// initialize
			local zminmse = .
			forvalues i=1/`numeqnsZ' {
				mata: `eqn'=*(`mname'.eqnlistZ[1,`i'])
				mata: st_local("vname",`eqn'.Vname)
				mata: st_local("vtilde",`eqn'.Vtilde)
				mata: st_local("command",`eqn'.command)
				mata: st_numscalar("r(MSE)",`eqn'.MSE)
				local MSE = `r(MSE)'
				mata: st_numscalar("r(N)",`eqn'.N)
				local MSE = `r(MSE)'
				local N = `r(N)'
				if "`var'"=="`vname'" {
					if `MSE' < `zminmse' {
						local zopt `vtilde'
						local zminmse `MSE'
					}
					di _col(2) "`vname'" _c
					di _col(20) "`vtilde'" _c
					di _col(40) "`command'" _c
					di _col(50) %6.0f `N' _c
					di _col(60) %10.6f `MSE'
				}
			}
			// if nameZopt already has vname in it, then add it to the list
			if `first_zopt' {
				mata: `mname'.nameZopt	= "`zopt'"
				local first_zopt	=0
			}
			else {
				mata: `mname'.nameZopt	= (`mname'.nameZopt, "`zopt'")
			}

		}
	}
		

/*
	*** save all 
	ereturn clear
	ereturn scalar crossfit = 1
	ereturn scalar yest = `yestn'
	ereturn scalar dest = `destn'
	ereturn local cmd ddml_crossfit
	ereturn local depvar `yvar'
	ereturn local dvar `dvar'
	ereturn local model "`model'"
	ereturn scalar crossfit = 1

	* return variable and command names
	forvalues i = 1(1)`yestn' {
		ereturn local y`i' `yname`i''
		ereturn local ycmd`i' `ycmd`i''
	}
	forvalues i = 1(1)`destn' {
		ereturn local d`i' `dname`i''
		ereturn local dcmd`i' `dcmd`i''
	}
	if ("`model'"=="iv"|"`model'"=="late") {
		ereturn scalar zest = `zestn'
		ereturn local zvar `zvar'
		forvalues i = 1(1)`zestn' {
			ereturn local z`i' `zname`i''
			ereturn local zcmd`i' `zcmd`i''
		}
	}

	* return MSE and opt-ID for Y
	if ("`model'"=="partial"|"`model'"=="iv") {
		ereturn scalar yoptid = `yminmseid'
		ereturn matrix ymse = `MSEy'
	}
	if ("`model'"=="late"|"`model'"=="interactive") {
		ereturn scalar y0optid = `yminmseid0'
		ereturn scalar y1optid = `yminmseid1'
		ereturn matrix y0mse = `MSEy0'
		ereturn matrix y1mse = `MSEy1'
	}
	*** return MSE and opt-ID for D
	if ("`model'"=="partial"|"`model'"=="iv"|"`model'"=="interactive") {
		ereturn scalar doptid = `dminmseid'
		ereturn matrix dmse = `MSEd'
	} 
	if ("`model'"=="late") {
		ereturn scalar d1optid = `dminmseid1'
		ereturn scalar d0optid = `dminmseid0'
		ereturn matrix d1mse = `MSEd1'
		ereturn matrix d0mse = `MSEd0'
	}
	*** return MSE and opt-ID for Z
	if ("`model'"=="late"|"`model'"=="iv") {
		ereturn scalar zoptid = `zminmseid'
		ereturn matrix zmse = `MSEz'
	}
*/
end

mata:

struct eqnStruct init_eqnStruct()
{
	struct eqnStruct scalar		e
	return(e)
}

void add_to_Yeqn(					struct ddmlStruct m,
									real scalar eqnumber,
									string scalar vnames,
									real scalar mse,
									real scalar n)
{
	pointer(struct eqnStruct) scalar p
	p				= m.eqnlistY[1,eqnumber]
	(*p).idVtilde	= st_data(., tokens(vnames))
	(*p).MSE		= mse
	(*p).N			= n
}

void add_to_Deqn(					struct ddmlStruct m,
									real scalar eqnumber,
									string scalar vnames,
									real scalar mse,
									real scalar n)
{
	pointer(struct eqnStruct) scalar p

	p				= m.eqnlistD[1,eqnumber]
	(*p).idVtilde	= st_data(., vnames)
	(*p).MSE		= mse
	(*p).N			= n
}

void add_to_Zeqn(					struct ddmlStruct m,
									real scalar eqnumber,
									string scalar vnames,
									real scalar mse,
									real scalar n)
{
	pointer(struct eqnStruct) scalar p

	p				= m.eqnlistZ[1,eqnumber]
	(*p).idVtilde	= st_data(., vnames)
	(*p).MSE		= mse
	(*p).N			= n
}


end
