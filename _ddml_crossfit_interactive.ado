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
	mata: st_numscalar("r(numeqns)",cols(`mname'.eqnlistY))
	local numeqnsY	= `r(numeqns)'
	mata: st_numscalar("r(numeqns)",cols(`mname'.eqnlistD))
	local numeqnsD	= `r(numeqns)'
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("listYtilde",invtokens(`mname'.nameYtilde))
	di "Number of Y estimating equations: `numeqnsY'"
	if `numeqnsD' {
		mata: st_local("listD",invtokens(`mname'.nameD))
		mata: st_local("listDtilde",invtokens(`mname'.nameDtilde))
		di "Number of D estimating equations: `numeqnsD'"
	}
	mata: st_numscalar("r(numeqns)",cols(`mname'.eqnlistZ))
	local numeqnsZ	= `r(numeqns)'
	if `numeqnsZ' {
		mata: st_local("listZ",invtokens(`mname'.nameZ))
		mata: st_local("listZtilde",invtokens(`mname'.nameZtilde))
		di "Number of Z estimating equations: `numeqnsZ'"
	}		 
	
	*** gen folds
	// use foldvar if not empty
	// populate provided foldvar if empty
	tempvar uni cuni
	qui gen double `uni' = runiform()
	qui cumul `uni', gen(`cuni')
	if "`foldvar'"~="" {
		// foldvar name provided so store on the model struct; can overwrite
		mata: `mname'.foldvar	= "`foldvar'"
	}
	else {
		// foldvar name not provided in foldvar(.) option so use default name
		// store on model struct and drop the variable if it already exists
		cap drop ddmlfold
		mata: `mname'.foldvar	= "ddmlfold"
	}
	// does variable exist and is it numeric?
	mata: st_local("kid",`mname'.foldvar)
	cap sum `kid'
	if _rc==0 {
		// variable exists and is numeric
		// does it have a valid number of groups?
		qui tab `kid'
		if r(r) < 2 {
			di as err "error - invalid fold variable"
			exit 1
		}
	}
	else {
		// create variable with provided name
		cap drop `kid'
		qui gen `kid' =ceil(`kfolds'*`cuni')
	}
	if ("`tabfold'"!="") {
		di ""
		di "Overview of frequencies by fold:"
		tab `kid' `listZ'
		tab `kid' `listD'
		di ""
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
		mata: st_local("vtilde",`eqn'.vtilde)
		cap drop `vtilde'
		qui gen double `vtilde'=.
	}
	// D equations
	forvalues i=1/`numeqnsD' {
		mata: `eqn'=*(`mname'.eqnlistD[1,`i'])
		mata: st_local("vtilde",`eqn'.vtilde)
		cap drop `vtilde'
		qui gen double `vtilde'=.
	}
	// Z equations
	if ("`model'"=="late") {
		forvalues i=1/`numeqnsZ' {
			mata: `eqn'=*(`mname'.eqnlistZ[1,`i'])
			mata: st_local("vtilde",`eqn'.vtilde)
			cap drop `vtilde'
			qui gen double `vtilde'=.
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
					mata: st_local("vtilde",`eqn'.vtilde)
					mata: st_local("vname",`eqn'.vname)
					mata: st_local("eststring",`eqn'.eststring)
					local 0 "`eststring'"
					syntax [anything] , [*]
					local est_main `anything'
					local est_options `options'
					di "Y estimating equation `i':"
					di "  est_main: `est_main'"
					di "  est_options: `est_options'"
					
					// for D = 1
					// estimate excluding kth fold
					`est_main' if `kid'!=`k', `est_options'
					// get fitted values and residuals for kth fold	
					tempvar vtilde_i
					qui predict double `vtilde_i' if `kid'==`k' & `listD' == 1
					qui replace `vtilde' = `vname' - `vtilde_i' if `kid'==`k' & `listD' == 1

					// for D = 0
					// estimate excluding kth fold
					`est_main' if `kid'!=`k', `est_options'
					// get fitted values and residuals for kth fold	
					tempvar vtilde_i
					qui predict double `vtilde_i' if `kid'==`k' & `listD' == 0
					qui replace `vtilde' = `vname' - `vtilde_i' if `kid'==`k' & `listD' == 0
				}
			}
			// D equations
			forvalues i=1/`numeqnsD' {
				mata: `eqn'=*(`mname'.eqnlistD[1,`i'])
				mata: st_numscalar("r(crossfit)",`eqn'.crossfit)
				if r(crossfit) {
					mata: st_local("vtilde",`eqn'.vtilde)
					mata: st_local("vname",`eqn'.vname)
					mata: st_local("eststring",`eqn'.eststring)
					local 0 "`eststring'"
					syntax [anything] , [*]
					local est_main `anything'
					local est_options `options'
					di "D estimating equation `i':"
					di "  est_main: `est_main'"
					di "  est_options: `est_options'"
					
					if ("`model'"=="interactive") {
						// estimate excluding kth fold
						`est_main' if `kid'!=`k', `est_options'
						// get fitted values and residuals for kth fold	
						tempvar vtilde_i
						qui predict double `vtilde_i' if `kid'==`k' 
						qui replace `vtilde' = `vname' - `vtilde_i' if `kid'==`k'
					}
					else {

						// for Z = 0
						// estimate excluding kth fold
						`est_main' if `kid'!=`k', `est_options'
						// get fitted values and residuals for kth fold	
						tempvar vtilde_i
						qui predict double `vtilde_i' if `kid'==`k' & `listZ' == 0
						qui replace `vtilde' = `vname' - `vtilde_i' if `kid'==`k' & `listZ' == 0

						// for Z = 1
						// estimate excluding kth fold
						`est_main' if `kid'!=`k', `est_options'
						// get fitted values and residuals for kth fold	
						tempvar vtilde_i
						qui predict double `vtilde_i' if `kid'==`k' & `listZ' == 1
						qui replace `vtilde' = `vname' - `vtilde_i' if `kid'==`k' & `listZ' == 1
					}
				}
			}
			// Z equations
			forvalues i=1/`numeqnsZ' {
				mata: `eqn'=*(`mname'.eqnlistZ[1,`i'])
				mata: st_numscalar("r(crossfit)",`eqn'.crossfit)
				if r(crossfit) {
					mata: st_local("vtilde",`eqn'.vtilde)
					mata: st_local("vname",`eqn'.vname)
					mata: st_local("eststring",`eqn'.eststring)
					local 0 "`eststring'"
					syntax [anything] , [*]
					local est_main `anything'
					local est_options `options'
					di "Z estimating equation `i':"
					di "  est_main: `est_main'"
					di "  est_options: `est_options'"
	
					// estimate excluding kth fold
					`est_main' if `kid'!=`k', `est_options'
					// get fitted values and residuals for kth fold	
					tempvar vtilde_i
					qui predict double `vtilde_i' if `kid'==`k' 
					qui replace `vtilde' = `vname' - `vtilde_i' if `kid'==`k'
				}
			}
		}
	}	

	*** calculate MSE (or any other postestimation calculations)
	
	// insert stats for each estimation

	if ("`model'"=="interactive") {
		forvalues i=1/`numeqnsY' {
			mata: `eqn'=*(`mname'.eqnlistY[1,`i'])
			mata: st_local("vtilde",`eqn'.vtilde)
			tempvar vtilde_sq
			qui gen double `vtilde_sq' = `vtilde'^2
			// for D=0
			qui sum `vtilde_sq' if `listD'==0, meanonly
			mata: add_stats_Y01(`mname',`i',`r(mean)',`r(N)',0)
			// for D=1
			qui sum `vtilde_sq' if `listD'==1, meanonly
			mata: add_stats_Y01(`mname',`i',`r(mean)',`r(N)',1)
		}
		forvalues i=1/`numeqnsD' {
			mata: `eqn'=*(`mname'.eqnlistD[1,`i'])
			mata: st_local("vtilde",`eqn'.vtilde)
			tempvar vtilde_sq
			qui gen double `vtilde_sq' = `vtilde'^2
			qui sum `vtilde_sq', meanonly
			mata: add_stats_D(`mname',`i',`r(mean)',`r(N)')
		}
	}
	if ("`model'"=="late") {
		forvalues i=1/`numeqnsY' {
			mata: `eqn'=*(`mname'.eqnlistY[1,`i'])
			mata: st_local("vtilde",`eqn'.vtilde)
			tempvar vtilde_sq
			qui gen double `vtilde_sq' = `vtilde'^2
			// for D=0
			qui sum `vtilde_sq' if `listZ'==0, meanonly
			mata: add_stats_Y01(`mname',`i',`r(mean)',`r(N)',0)
			// for D=1
			qui sum `vtilde_sq' if `listZ'==1, meanonly
			mata: add_stats_Y01(`mname',`i',`r(mean)',`r(N)',1)
		}
		forvalues i=1/`numeqnsD' {
			mata: `eqn'=*(`mname'.eqnlistD[1,`i'])
			mata: st_local("vtilde",`eqn'.vtilde)
			tempvar vtilde_sq
			qui gen double `vtilde_sq' = `vtilde'^2
			// for Z = 0
			sum `vtilde_sq' if `listZ'==0, meanonly
			mata: add_stats_D01(`mname',`i',`r(mean)',`r(N)',0)
			// for Z = 1
			sum `vtilde_sq' if `listZ'==1, meanonly
			mata: add_stats_D01(`mname',`i',`r(mean)',`r(N)',1)
		}
		forvalues i=1/`numeqnsZ' {
			mata: `eqn'=*(`mname'.eqnlistZ[1,`i'])
			mata: st_local("vtilde",`eqn'.vtilde)
			tempvar vtilde_sq
			qui gen double `vtilde_sq' = `vtilde'^2
			qui sum `vtilde_sq', meanonly
			mata: add_stats_Z(`mname',`i',`r(mean)',`r(N)')
		}
	}
	
	*** print results & find optimal model

	// Y|X,D=j or Y|X,Z=j
	forvalues j= 0/1 {
		di
	    if "`model'"=="interactive" di as res "Mean-squared error for y|X,D=`j':"
	    else  di as res "Mean-squared error for y|X,Z=`j':"
		di _col(2) "Name" _c
		di _col(20) "Orthogonalized" _c
		di _col(40) "Command" _c
		di _col(54) "N" _c
		di _col(65) "MSPE"
		di "{hline 75}"
		// initialize
		local yminmse`j' = .
		local yopt`j'
		forvalues i=1/`numeqnsY' {
			mata: `eqn'=*(`mname'.eqnlistY[1,`i'])
			mata: st_local("vname",`eqn'.vname)
			mata: st_local("vtilde",`eqn'.vtilde)
			mata: st_local("command",`eqn'.command)
			mata: st_numscalar("r(MSE)",`eqn'.MSE`j')
			mata: st_numscalar("r(N)",`eqn'.N`j')
			local MSE = `r(MSE)'
			local N = `r(N)'
			if `MSE' < `yminmse`j'' {
				local yopt`j' `vtilde'
				local yminmse`j' `MSE'
			}
			di _col(2) "`vname'" _c
			di _col(20) "`vtilde'" _c
			di _col(40) "`command'" _c
			di _col(50) %6.0f `N' _c
			di _col(60) %10.6f `MSE'
		}
		mata: `mname'.nameY`j'opt	= "`yopt`j''"
	}

	// D vars
	if "`model'"=="interactive" {
		di
		di as res "Mean-squared error for D|X:"
		di _col(2) "Name" _c
		di _col(20) "Orthogonalized" _c
		di _col(40) "Command" _c
		di _col(54) "N" _c
		di _col(65) "MSPE"
		di "{hline 75}"
		// initialize
		foreach var of varlist `listD' {
			// initialize
			local dminmse = .
			local dopt
			forvalues i=1/`numeqnsD' {
				mata: `eqn'=*(`mname'.eqnlistD[1,`i'])
				mata: st_local("vname",`eqn'.vname)
				mata: st_local("vtilde",`eqn'.vtilde)
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
			mata: `mname'.nameDopt	= "`dopt'"
		}
	}
	if "`model'"=="late" {
		forvalues j= 0/1 {
			di
			di as res "Mean-squared error for D|X,Z=`j':"
			di _col(2) "Name" _c
			di _col(20) "Orthogonalized" _c
			di _col(40) "Command" _c
			di _col(54) "N" _c
			di _col(65) "MSPE"
			di "{hline 75}"
			// initialize
			local dminmse`j' = .
			local dopt`j'
			forvalues i=1/`numeqnsD' {
				mata: `eqn'=*(`mname'.eqnlistD[1,`i'])
				mata: st_local("vname",`eqn'.vname)
				mata: st_local("vtilde",`eqn'.vtilde)
				mata: st_local("command",`eqn'.command)
				mata: st_numscalar("r(MSE)",`eqn'.MSE`j')
				mata: st_numscalar("r(N)",`eqn'.N`j')
				local MSE = `r(MSE)'
				local N = `r(N)'
				if `MSE' < `dminmse`j'' {
					local dopt`j' `vtilde'
					local dminmse`j' `MSE'
				}
				di _col(2) "`vname'" _c
				di _col(20) "`vtilde'" _c
				di _col(40) "`command'" _c
				di _col(50) %6.0f `N' _c
				di _col(60) %10.6f `MSE'
			}
			mata: `mname'.nameD`j'opt	= "`dopt`j''"
		}
	}

	// Z vars
	if ("`model'"=="late") {
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
			local zopt
			forvalues i=1/`numeqnsZ' {
				mata: `eqn'=*(`mname'.eqnlistZ[1,`i'])
				mata: st_local("vname",`eqn'.vname)
				mata: st_local("vtilde",`eqn'.vtilde)
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
			mata: `mname'.nameZopt	= "`zopt'"
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


void add_stats_Y01(					struct ddmlStruct m,
									real scalar eqnumber,
									real scalar mse,
									real scalar n,
									real scalar D)
{
	pointer(struct eqnStruct) scalar p

	p = m.eqnlistY[1,eqnumber]
	if (D==0) {
		(*p).MSE0	= mse
		(*p).N0		= n
	}
	else {
		(*p).MSE1	= mse
		(*p).N1		= n
	}
}

void add_stats_D(					struct ddmlStruct m,
									real scalar eqnumber,
									real scalar mse,
									real scalar n)
{
	pointer(struct eqnStruct) scalar p

	p = m.eqnlistD[1,eqnumber]
	(*p).MSE	= mse
	(*p).N		= n
}

void add_stats_D01(					struct ddmlStruct m,
									real scalar eqnumber,
									real scalar mse,
									real scalar n,
									real scalar Z)
{
	pointer(struct eqnStruct) scalar p

	p = m.eqnlistD[1,eqnumber]
	if (Z==0) {
		(*p).MSE0	= mse
		(*p).N0		= n
	}
	else {
		(*p).MSE1	= mse
		(*p).N1		= n
	}
}

void add_stats_Z(					struct ddmlStruct m,
									real scalar eqnumber,
									real scalar mse,
									real scalar n)
{
	pointer(struct eqnStruct) scalar p

	p = m.eqnlistZ[1,eqnumber]
	(*p).MSE	= mse
	(*p).N		= n
}


end
