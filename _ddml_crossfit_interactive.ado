*** ddml cross-fitting for the interactive model & LATE
program _ddml_crossfit_interactive, eclass sortpreserve

	syntax [anything] [if] [in] , /// 
							[ kfolds(integer 2) ///
							NOIsily ///
							debug /// 
							Robust ///
							TABFold ///
							foldvar(name) ///
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

	*** initialize tilde variables
	forvalues i=1/`numeqns' {
		mata: `eqn'=*(`mname'.eqnlist[1,`i'])
		mata: st_local("vtilde",`eqn'.Vtilde)
		cap drop `mname'_`vtilde'
		qui gen double `mname'_`vtilde'=.
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
		
			// loop over equations
			forvalues i=1/`numeqns' {
				mata: `eqn'=*(`mname'.eqnlist[1,`i'])
				mata: st_local("eqntype",`eqn'.eqntype)
				if ("`eqntype'"=="yeq") {
					// Y equations
					mata: st_local("vtilde",`eqn'.Vtilde)
					mata: st_local("vname",`eqn'.Vname)
					mata: st_local("eststring",`eqn'.eststring)
					mata: st_local("vtype",`eqn'.vtype)
					local 0 "`eststring'"
					syntax [anything] , [*]
					local est_main `anything'
					local est_options `options'
					di "Y estimating equation `i':"
					di "  est_main: `est_main'"
					di "  est_options: `est_options'"
					
					// for D = 1
					// estimate excluding kth fold
					`est_main' if `mname'_fid!=`k' & `listD' == 1 & `mname'_sample, `est_options'
					// get fitted values and residuals for kth fold	
					tempvar vtilde_i
					qui predict `vtype'`vtilde_i' if `mname'_fid==`k' & `listD' == 1 & `mname'_sample
					qui replace `mname'_`vtilde' = `vname' - `vtilde_i' if `mname'_fid==`k' & `listD' == 1 & `mname'_sample

					// for D = 0
					// estimate excluding kth fold
					`est_main' if `mname'_fid!=`k' & `listD' == 0 & `mname'_sample, `est_options'
					// get fitted values and residuals for kth fold	
					tempvar vtilde_i
					qui predict `vtype' `vtilde_i' if `mname'_fid==`k' & `listD' == 0 & `mname'_sample
					qui replace `mname'_`vtilde' = `vname' - `vtilde_i' if `mname'_fid==`k' & `listD' == 0 & `mname'_sample
				
				}
				else if ("`eqntype'"=="deq") {
					// D equations
					mata: `eqn'=*(`mname'.eqnlist[1,`i'])
					mata: st_local("vtilde",`eqn'.Vtilde)
					mata: st_local("vname",`eqn'.Vname)
					mata: st_local("eststring",`eqn'.eststring)
					mata: st_local("vtype",`eqn'.vtype)
					local 0 "`eststring'"
					syntax [anything] , [*]
					local est_main `anything'
					local est_options `options'
					di "D estimating equation `i':"
					di "  est_main: `est_main'"
					di "  est_options: `est_options'"
						
					if ("`model'"=="interactive") {
						// estimate excluding kth fold
						`est_main' if `mname'_fid!=`k' & `mname'_sample, `est_options'
						// get fitted values and residuals for kth fold	
						tempvar vtilde_i
						qui predict `vtype' `vtilde_i' if `mname'_fid==`k' & `mname'_sample
						qui replace `mname'_`vtilde' = `vname' - `vtilde_i' if `mname'_fid==`k' & `mname'_sample
					}
					else {

						// for Z = 0
						// estimate excluding kth fold
						`est_main' if `mname'_fid!=`k' & `listZ' == 0 & `mname'_sample, `est_options'
						// get fitted values and residuals for kth fold	
						tempvar vtilde_i
						qui predict `vtype' `vtilde_i' if `mname'_fid==`k' & `listZ' == 0 & `mname'_sample
						qui replace `mname'_`vtilde' = `vname' - `vtilde_i' if `mname'_fid==`k'  & `listZ' == 0 & `mname'_sample

						// for Z = 1
						// estimate excluding kth fold
						`est_main' if `mname'_fid!=`k' & `listZ' == 1 & `mname'_sample, `est_options'
						// get fitted values and residuals for kth fold	
						tempvar vtilde_i
						qui predict `vtype' `vtilde_i' if `mname'_fid==`k' & `listZ' == 1 & `mname'_sample
						qui replace `mname'_`vtilde' = `vname' - `vtilde_i' if `mname'_fid==`k' & `listZ' == 1 & `mname'_sample
					}
				}
				else if ("`eqntype'"=="zeq") {

					// Z equations
					mata: `eqn'=*(`mname'.eqnlist[1,`i'])
					mata: st_local("vtilde",`eqn'.Vtilde)
					mata: st_local("vname",`eqn'.Vname)
					mata: st_local("eststring",`eqn'.eststring)
					mata: st_local("vtype",`eqn'.vtype)
					local 0 "`eststring'"
					syntax [anything] , [*]
					local est_main `anything'
					local est_options `options'
					di "Z estimating equation `i':"
					di "  est_main: `est_main'"
					di "  est_options: `est_options'"
		
					// estimate excluding kth fold
					`est_main' if `mname'_fid!=`k' & `mname'_sample, `est_options'
					// get fitted values and residuals for kth fold	
					tempvar vtilde_i
					qui predict `vtype' `vtilde_i' if `mname'_fid==`k' & `mname'_sample
					qui replace `mname'_`vtilde' = `vname' - `vtilde_i' if `mname'_fid==`k' & `mname'_sample
				}
			}
		}
 	}	

	*** calculate MSE (or any other postestimation calculations)
	if ("`model'"=="interactive") {
		*** calculate MSE, store orthogonalized variables, etc.
		forvalues i=1/`numeqns' {
			mata: `eqn'=*(`mname'.eqnlist[1,`i'])
			mata: st_local("eqntype",`eqn'.eqntype)
			mata: st_local("vtilde",`eqn'.Vtilde)
			tempvar vtilde_sq
			qui gen double `vtilde_sq' = `mname'_`vtilde'^2 if `mname'_sample
			if ("`eqntype'"=="yeq") {
				// for D=1
				sum `vtilde_sq' if `listD'==0 & `mname'_sample, meanonly
				mata: add_to_eqn01(`mname',`i',"`mname'_id `mname'_`vtilde'", `r(mean)',`r(N)',0)
				// for D=1
				qui sum `vtilde_sq' if `listD'==1 & `mname'_sample, meanonly
				mata: add_to_eqn01(`mname',`i',"`mname'_id `mname'_`vtilde'", `r(mean)',`r(N)',1)
			} 
			else if ("`eqntype'"=="deq") {
				qui sum `vtilde_sq' if `mname'_sample, meanonly
				mata: add_to_eqn(`mname',`i',"`mname'_id `mname'_`vtilde'", `r(mean)',`r(N)')
			}
		}
	}
	if ("`model'"=="late") {
		*** calculate MSE, store orthogonalized variables, etc.
		forvalues i=1/`numeqns' {
			mata: `eqn'=*(`mname'.eqnlist[1,`i'])
			mata: st_local("eqntype",`eqn'.eqntype)
			mata: st_local("vtilde",`eqn'.Vtilde)
			tempvar vtilde_sq
			qui gen double `vtilde_sq' = `mname'_`vtilde'^2 if `mname'_sample
			if ("`eqntype'"=="yeq") {
				// for Z=1
				qui sum `vtilde_sq' if `listZ'==0 & `mname'_sample, meanonly
				mata: add_to_eqn01(`mname',`i',"`mname'_id `mname'_`vtilde'", `r(mean)',`r(N)',0)
				// for Z=1
				qui sum `vtilde_sq' if `listZ'==1 & `mname'_sample, meanonly
				mata: add_to_eqn01(`mname',`i',"`mname'_id `mname'_`vtilde'", `r(mean)',`r(N)',1)
			} 
			else if ("`eqntype'"=="deq") {
				// for Z = 0
				sum `vtilde_sq' if `listZ'==0 & `mname'_sample, meanonly
				mata: add_to_eqn01(`mname',`i',"`mname'_id `mname'_`vtilde'", `r(mean)',`r(N)',0)
				// for Z = 1
				sum `vtilde_sq' if `listZ'==1 & `mname'_sample, meanonly
				mata: add_to_eqn01(`mname',`i',"`mname'_id `mname'_`vtilde'", `r(mean)',`r(N)',1)
			}
			else if ("`eqntype'"=="zeq") {
				qui sum `vtilde_sq' if `mname'_sample, meanonly
				mata: add_to_eqn(`mname',`i',"`mname'_id `mname'_`vtilde'", `r(mean)',`r(N)')
			}
		}
	}
 
	*** print results & find optimal model

	// interactive model
	if "`model'"=="interactive" {

		// dependent variable
		_ddml_display_header , str(y|X,D=0)
		_ddml_display_mspe `mname', vname(`nameY') zett(0)
		mata: `mname'.nameY1opt		= "`r(optname)'"
		_ddml_display_header , str(y|X,D=1)
		_ddml_display_mspe `mname', vname(`nameY') zett(1)
		mata: `mname'.nameY0opt		= "`r(optname)'"

		// D variable
		_ddml_display_header , str(D|X)
		foreach var of varlist `listD' {
			_ddml_display_mspe `mname', vname(`var') 
			mata: `mname'.nameDopt		= "`r(optname)'"
		}
	}

	// late model
	if "`model'"=="late" {

		// dependent variable
		_ddml_display_header , str(y|X,Z=0)
		_ddml_display_mspe `mname', vname(`nameY') zett(0)
		mata: `mname'.nameY1opt		= "`r(optname)'"
		_ddml_display_header , str(y|X,Z=1)
		_ddml_display_mspe `mname', vname(`nameY') zett(1)
		mata: `mname'.nameY0opt		= "`r(optname)'"

		// D variable
		_ddml_display_header , str(D|X,Z=0)
		foreach var of varlist `listD' {
			_ddml_display_mspe `mname', vname(`var') zett(0)
			mata: `mname'.nameD0opt		= "`r(optname)'"
		}
		_ddml_display_header , str(D|X,Z=1)
		foreach var of varlist `listD' {
			_ddml_display_mspe `mname', vname(`var') zett(1)
			mata: `mname'.nameD1opt		= "`r(optname)'"
		}

		// Z variable
		foreach var of varlist `listZ' {
			_ddml_display_mspe `mname', vname(`var')
			mata: `mname'.nameZopt = (`mname'.nameZopt, "`r(optname)'")
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

void add_to_eqn01(					struct ddmlStruct m,
									real scalar eqnumber,
									string scalar vnames,
									real scalar mse,
									real scalar n, 
									real scalar Z)
{
	pointer(struct eqnStruct) scalar p
	p				= m.eqnlist[1,eqnumber]
	(*p).idVtilde	= st_data(., tokens(vnames))
	if (Z==0) {
		(*p).MSE0		= mse
		(*p).N0			= n
	}
	else {
		(*p).MSE1		= mse
		(*p).N1			= n
	}

}

end
