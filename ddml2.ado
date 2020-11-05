
program ddml2, eclass

	version 13

	syntax [anything] , [*]

	local subcmd: word 1 of `anything'

	*** get latest version
	if "`subcmd'"=="update" {
		net install ddml, from(https://raw.githubusercontent.com/aahrens1/ddml/master/)
	} 
	
	*** describe model
	if "`subcmd'"=="desc" {
		local 0 ", `options'"
		// mname is required; could make optional with a default name
		// remaining args are temporary and for debugging only
		syntax , mname(name)

		mata: st_global("r(mstring)",`mname'.model)
		di "Model: `r(mstring)'"
		mata: st_global("r(mstring)",`mname'.nameY)
		di "Dependent variable (Y): `r(mstring)'"
		mata: st_global("r(mstring)",invtokens(`mname'.nameYtilde))
		di "Dependent variable (orthogonalized): `r(mstring)'"
		mata: st_global("r(mstring)",invtokens(`mname'.nameYopt))
		di "Minimum MSE orthogonalized dep var: `r(mstring)'"
		mata: st_global("r(mstring)",invtokens(`mname'.nameD))
		di "Causal variable (D): `r(mstring)'"
		mata: st_global("r(mstring)",invtokens(`mname'.nameDtilde))
		di "Causal variable (orthogonalized): `r(mstring)'"
		mata: st_global("r(mstring)",invtokens(`mname'.nameDopt))
		di "Minimum MSE orthogonalized causal var: `r(mstring)'"
		// List equations
		mata: st_numscalar("r(numeqns)",cols(`mname'.eqnlist))
		local numeqns	= `r(numeqns)'
		di "Number of estimating equations: `numeqns'"
		
		// blank eqn - declare this way so that it's a struct and not transmorphic
		tempname eqn
		mata: `eqn' = init_eqnStruct()
		forvalues i=1/`numeqns' {
			mata: `eqn'=*(`mname'.eqnlist[1,`i'])
			di "Estimating equation `i': " _c
			mata: st_numscalar("r(MSE)",`eqn'.MSE)
			di "MSE = " %10.6f `r(MSE)'
			mata: st_global("r(estring)",`eqn'.YorD)
			di "  Type: `r(estring)'" _c
			mata: st_global("r(estring)",`eqn'.vname)
			di _col(15) "Variable: `r(estring)'" _c
			mata: st_global("r(estring)",`eqn'.vtilde)
			di _col(40) "Orthogonalized: `r(estring)'"
			mata: st_global("r(estring)",`eqn'.eststring)
			di "  Command: `r(estring)'"
		}
		// clear this global from Mata
		mata: mata drop `eqn'

	}
	
	*** initialize new estimation
	if "`subcmd'"=="init" {
		local model: word 2 of `anything'
		if ("`model'"!="partial"&"`model'"!="iv"&"`model'"!="interactive"&"`model'"!="late"&"`model'"!="optimaliv") {
			di as err "no or wrong model specified." 
			exit 1
		}
		local 0 ", `options'"
		// mname is required; could make optional with a default name
		// remaining args are temporary and for debugging only
		syntax , mname(name)
		mata: `mname'=init_ddmlStruct()
		// fill by hand
		mata: `mname'.model		= "`model'"

	} 

	*** add equation  
	if "`subcmd'"=="yeq"|"`subcmd'"=="deq"|"`subcmd'"=="zeq" {

		** check that equation is consistent with model
		if ("`subcmd'"=="yeq"&"`model'"=="optimaliv") {
			di as err "not allowed; yeq not allowed with `model'"
		}
		if ("`subcmd'"=="zeq"&("`model'"=="optimaliv"|"`model'"=="partial"|"`model'"=="interactive")) {
			di as err "not allowed; deq not allowed with `model'"
		}

		** check that ddml has been initialized
		// to add

		** parsing
		local 0 ", `options'"
		syntax ,	mname(name)		///
					eqn(string)		///
					vname(name)		///
					gen(name)

		mata: add_eqn(`mname', "`subcmd'", "`vname'", "`gen'", "`eqn'")
		if "`subcmd'"=="yeq" {
			// check if nameY is already there; if it is, must be identical to vname here
			mata: st_global("r(vname)",`mname'.nameY)
			if "`r(vname)'"=="" {
				mata: `mname'.nameY		= "`vname'"
			}
			else if "`r(vname)'"~="`vname'" {
				di as err "error - incompatible y variables"
				exit 198
			}

		}
		if "`subcmd'"=="deq" {
			// check if nameD already has vname; if not, add it to the list
			mata: st_global("r(vname)",invtokens(`mname'.nameD))
			if "`r(vname)'"=="" {
				mata: `mname'.nameD		= "`vname'"
			}
			else {
				local dlist `r(vname)' `vname'
				local dlist : list uniq dlist
				mata: `mname'.nameD		= tokens("`dlist'")
			}
		}

	}

	*** cross-fitting
	if "`subcmd'" =="crossfit" {
		_ddml_crossfit, `options'
	}

	*** estimate
	if "`subcmd'" =="estimate" {
	
		local mstruct: word 2 of `anything'
		
		// check that mstruct is the name of a Mata ddmlStruct
		mata: st_global("r(structname)",structname(`mstruct'))
		if ("`r(structname)'" ~= "ddmlStruct") {
			di as err "you need to provide the name of the ddmlStruct with the crossfit results"
			exit 198
		}

		mata: st_global("r(model)",`mstruct'.model)

		if ("`r(model)'"=="partial") {
			_ddml_estimate_partial `mstruct', `options'
		}
		if ("`r(model)'"=="iv") {
			_ddml_estimate_iv, `options'
		}
		if ("`r(model)'"=="interactive") {
			_ddml_estimate_interactive, `options'
		}
		if ("`r(model)'"=="late") {
			_ddml_estimate_late, `options'
		}
		if ("`r(model)'"=="optimaliv") {
			_ddml_estimate_optimaliv, `options'
		}
	}
end

*** ddml cross-fitting
program _ddml_crossfit, eclass sortpreserve

	syntax [anything] [if] [in] , /// 
							[ kfolds(integer 2)  ///
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
	
	*** extract details of estimation
	
	// model
	mata: st_local("model",`mname'.model)
	mata: st_numscalar("r(numeqns)",cols(`mname'.eqnlist))
	local numeqns	= `r(numeqns)'
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("listYtilde",invtokens(`mname'.nameYtilde))
	mata: st_local("listD",invtokens(`mname'.nameD))
	mata: st_local("listDtilde",invtokens(`mname'.nameDtilde))
	di "Model: `model'"
	di "Number of estimating equations: `numeqns'"
	
	*** gen folds
	tempvar kid uni cuni
	gen double `uni' = runiform()
	cumul `uni', gen(`cuni')
	gen `kid' =ceil(`kfolds'*`cuni')
	if ("`tabfold'"!="") {
		di ""
		di "Overview of frequencies by fold:"
		tab `kid'
		di ""
	}
	//

	*** do cross-fitting
	di "Cross-fitting fold " _c
	forvalues k = 1(1)`kfolds' {
	
		if (`k'==`kfolds') {
			di as text "`k'"
		}
		else {
			di as text "`k' " _c
		}
		// ML is applied to I^c sample (all data ex partition k)
		qui {
		
			// blank eqn - declare this way so that it's a struct and not transmorphic
			tempname eqn
			mata: `eqn' = init_eqnStruct()
			forvalues i=1/`numeqns' {
				di "Estimating equation `i':"
				mata: `eqn'=*(`mname'.eqnlist[1,`i'])
				mata: st_local("vtilde",`eqn'.vtilde)
				mata: st_local("vname",`eqn'.vname)
				mata: st_local("eststring",`eqn'.eststring)
				local 0 "`eststring'"
				syntax [anything] , [*]
				local est_main `anything'
				local est_options `options'
				di "  est_main: `est_main'"
				di "  est_options: `est_options'"

				// estimate excluding kth fold
				cap drop `vtilde'
				qui gen double `vtilde'=.
				`est_main' if `kid'!=`k', `est_options'
				// get fitted values and residuals for kth fold	
				tempvar vtilde_i
				qui predict double `vtilde_i' if `kid'==`k' 
				qui replace `vtilde' = `vname' - `vtilde_i' if `kid'==`k'

			}
			// clear this global from Mata
			mata: mata drop `eqn'

		}
	}	

	*** calculate MSE 
	if ("`model'"=="partial"|"`model'"=="iv") {

		// blank eqn - declare this way so that it's a struct and not transmorphic
		tempname eqn
		mata: `eqn' = init_eqnStruct()
		
		// insert MSE for each estimation
		forvalues i=1/`numeqns' {
			mata: `eqn'=*(`mname'.eqnlist[1,`i'])
			mata: st_local("vtilde",`eqn'.vtilde)
			tempvar vtilde_sq
			qui gen double `vtilde_sq' = `vtilde'^2
			qui sum `vtilde_sq', meanonly
			mata: add_mse(`mname',`i',`r(mean)')
		}
		
		// loop through equations and pick out those with smallest MSE

		// dep var
		di
		di as res "Mean-squared error for y|X:"
		di _col(2) _c
		di _col(10) "Name" _c
		di _col(30) "Command" _c
		di _col(50) "MSPE"
		di "{hline 65}"
		// initialize
		local yminmse = .
		forvalues i=1/`numeqns' {
			mata: `eqn'=*(`mname'.eqnlist[1,`i'])
			mata: st_local("vname",`eqn'.vname)
			mata: st_local("vtilde",`eqn'.vtilde)
			mata: st_local("YorD",`eqn'.YorD)
			mata: st_local("command",`eqn'.command)
			mata: st_numscalar("r(MSE)",`eqn'.MSE)
			local MSE = `r(MSE)'
			if "`YorD'"=="yeq" {
				if `MSE' < `yminmse' {
					local yopt `vtilde'
					local yminmse `MSE'
				}
				di _col(10) "`vname'" _c
				di _col(30) "`command'" _c
				di _col(45) %10.6f `MSE'
			}
		}
		mata: `mname'.nameYopt		= "`yopt'"

		// D vars
		di
		di as res "Mean-squared error for y|X:"
		di _col(2) _c
		di _col(10) "Name" _c
		di _col(30) "Command" _c
		di _col(50) "MSPE"
		di "{hline 65}"
		// initialize
		local first_dopt	=1
		foreach var of varlist `listD' {
			// initialize
			local dminmse = .
			forvalues i=1/`numeqns' {
				mata: `eqn'=*(`mname'.eqnlist[1,`i'])
				mata: st_local("vname",`eqn'.vname)
				mata: st_local("vtilde",`eqn'.vtilde)
				mata: st_local("command",`eqn'.command)
				mata: st_numscalar("r(MSE)",`eqn'.MSE)
				local MSE = `r(MSE)'
				if "`var'"=="`vname'" {
					if `MSE' < `dminmse' {
						local dopt `vtilde'
						local dminmse `MSE'
					}
					di _col(10) "`vname'" _c
					di _col(30) "`command'" _c
					di _col(45) %10.6f `MSE'
				}
			}

			// if nameDopt already has vname in it, then add it to the list
			if `first_dopt' {
				mata: `mname'.nameDopt	= "`dopt'"
				local first_dopt	=0
			}
			else {
				mata: `mname'.nameDopt	= tokens("`r(vname)' `dopt'")
			}
		}
		
		
	}

/*
		tempname MSEy
		local yminmse = .
		local yminmseid = 1
		mat `MSEy' = J(1,`yestn',.)
		forvalues i = 1(1)`yestn' {
			qui replace `yname`i'' = `ytilde`i''
			tempname ytilde_sq`i'
			qui gen double `ytilde_sq`i'' = (`ytilde`i'')^2
			qui sum `ytilde_sq`i'' , meanonly
			mat `MSEy'[1,`i'] = r(mean)
			local newyminmse = r(mean)
			if (`newyminmse'<`yminmse') {
				local yminmse = `newyminmse'
				local yminmseid = `i'
			}
		}
*/
/*
		di " "
		di as res "Mean-squared error for y|X:"
		di _col(2) _c
		di _col(10) "Name" _c
		di _col(30) "Command" _c
		di _col(50) "MSPE"
		di "{hline 65}"
		forvalues i = 1(1)`yestn' {
			di _col(2) "`i'" _c
			di _col(10) "`yname`i''" _c
			di _col(30) "`ycmd`i''" _c
			if (`yminmseid'==`i') {
				di _col(50) `MSEy'[1,`i'] _c
				di "*"
			} 
			else {
				di _col(50) `MSEy'[1,`i'] 
			}
		}
*/

//	di as text "* indicates model with minimum MSE."


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

*** ddml estimation: partial linear model
program _ddml_estimate_partial, eclass sortpreserve

	syntax namelist(name=mstruct) [if] [in] , /// 
								[  ///
								ROBust ///
								show(string) /// dertermines which to post
								clear /// deletes all tilde-variables (to be implemented)
								avplot ///
								* ]

/*

	*** set defaults
	if ("`show'"=="") {
		local show opt
	}
	
	*** save everything that is needed in locals
	local yestn = e(yest)
	local destn = e(dest)
	local zestn = e(zest)
	local yoptid = e(yoptid)
	local doptid = e(doptid)
	local zoptid = e(zoptid)
	local yvar = e(depvar)
	local dvar = e(dvar)
	local zvar = e(zvar)
	
	*** retrieve variable names
	forvalues i = 1(1)`yestn' {
		local ytilde`i' `e(y`i')'
		local yname`i' `e(y`i')'
		local ycmd`i' `e(ycmd`i')'
	}
	forvalues i = 1(1)`destn' {
		local dtilde`i' `e(d`i')'
		local dname`i' `e(d`i')'
		local dcmd`i' `e(dcmd`i')'
	}

	*** do estimation
	if ("`show'"=="all") {
		forvalues i = 1(1)`yestn' {
			forvalues j = 1(1)`destn' {
				if (`i'==`yoptid' & `j'==`doptid') {
					// do nothing: optimal model always comes last 
					// and is estimated below
					di "" _c
				}
				else {
					di as text "DML with `ycmd`i'' (`yname`i'') and `dcmd`j'' (`dname`j''):"
					qui reg `ytilde`i'' `dtilde`j'', nocons `robust' noheader
					// display
					tempname b
					tempname V 
					mat `b' = e(b)
					mat `V' = e(V)
					matrix colnames `b' = "`dvar'"
					matrix rownames `b' = "`yvar'"
					matrix colnames `V' = "`dvar'"
					matrix rownames `V' = "`dvar'"
					ereturn clear
					ereturn post `b' `V' 
					ereturn display
				}
			}
		}
	}
*/

   	mata: st_global("r(Yopt)",`mstruct'.nameYopt)
	local Yopt `r(Yopt)'
   	mata: st_global("r(Dopt)",`mstruct'.nameDopt)
	local Dopt `r(Dopt)'
   	mata: st_global("r(nameY)",`mstruct'.nameY)
	local nameY `r(nameY)'
   	mata: st_global("r(nameD)",`mstruct'.nameD)
	local nameD `r(nameD)'

	*** estimate best model
	di as res "Optimal model: DML with optimal Y~ `r(Yopt)' and optimal D~ `r(Dopt)':"
	qui reg `Yopt' `Dopt', nocons `robust' noheader

	// plot
	if ("`avplot'"!="") {
	   twoway (scatter `ytilde`yoptid'' `dtilde`doptid'') (lfit `ytilde`yoptid'' `dtilde`doptid'')
	}

	// display
	tempname b
	tempname V 
	mat `b' = e(b)
	mat `V' = e(V)
	matrix colnames `b' = `nameD'
	matrix rownames `b' = `nameY'
 	matrix colnames `V' = `nameD'
	matrix rownames `V' = `nameD'
	ereturn clear
	ereturn post `b' `V'
	if "`robust'"~="" {
		ereturn local vcetype	robust
	}
	ereturn display

end



********************************************************************************
*** Mata section															 ***
********************************************************************************

mata:

struct ddmlStruct init_ddmlStruct()
{
	struct ddmlStruct scalar	d
	d.eqnlist		= J(0,0,NULL)
	d.nameY			= ""
	d.nameYtilde	= ""
	d.nameYopt		= ""
	d.nameD			= ""
	d.nameDtilde	= ""
	d.nameDopt		= ""
	return(d)
}

struct eqnStruct init_eqnStruct()
{
	struct eqnStruct scalar		e
	return(e)
}

void add_mse(						struct ddmlStruct m,
									real scalar eqnumber,
									real scalar mse)
{
	pointer(struct eqnStruct) scalar p

	p = m.eqnlist[1,eqnumber]
	(*p).MSE	= mse

}

void add_eqn(						struct ddmlStruct m,
									string scalar YorD,
									string scalar vname,
									string scalar vtilde,
									string scalar estcmd)
{
	struct eqnStruct scalar		e
	e.YorD		= YorD
	e.vname		= vname
	e.vtilde	= vtilde
	e.eststring	= estcmd
	e.command	= tokens(estcmd)[1,1]

	if (cols(m.eqnlist)==0) {
		m.eqnlist	= &e
	}
	else {
		m.eqnlist	= (m.eqnlist, &e)
	}

	if (YorD=="yeq") {
		if (m.nameYtilde=="") {
			m.nameYtilde	= vtilde
		}
		else {
			m.nameYtilde	= (m.nameYtilde, vtilde)
		}
	}
	else if (YorD=="deq") {
		if (m.nameDtilde=="") {
			m.nameDtilde	= vtilde
		}
		else {
			m.nameDtilde	= (m.nameDtilde, vtilde)
		}
	}
}

void ATE(   string scalar yvar,		// Y
			string scalar dvar,	 // D
			string scalar y0tilde,  // E[Y|X,D=0]
			string scalar y1tilde,  // E[Y|X,D=1]
			string scalar dtilde,   // E[D|X]
			string scalar sample,   // sample
			string scalar outate,   // output: name of matrix to store b
			string scalar outatese  // output: name of matrix to store V
			)
{
	st_view(my_d0x,.,y0tilde,sample)
	st_view(my_d1x,.,y1tilde,sample)
	st_view(md_x,.,dtilde,sample)
	st_view(d,.,dvar,sample)
	st_view(y,.,yvar,sample)

	n = rows(y)

	my_d0x = my_d0x :* (1:-d)
	my_d1x = my_d1x :* d

	te  = (d :* (y :- my_d1x) :/ md_x) :-  ((1 :- d) :* (y :- my_d0x) :/ (1 :- md_x)) :+ my_d1x :- my_d0x  
	ate = mean(te)
	ate_V =  variance(te)/n

	st_matrix(outate,ate)
	st_matrix(outatese,ate_V)
}

void LATE(  string scalar yvar,		 // Y
			string scalar dvar,	  // D
			string scalar zvar,	  // Z
			string scalar y0tilde,   // E[Y|X,Z=0]
			string scalar y1tilde,   // E[Y|X,Z=1]
			string scalar d0tilde,   // E[D|X,Z=0]
			string scalar d1tilde,   // E[D|X,Z=1]
			string scalar ztilde,	// E[Z|X]
			string scalar sample,	// sample
			string scalar outlate,   // output: name of matrix to store b
			string scalar outlatese  // output: name of matrix to store V
			)
{
	st_view(my_z0x,.,y0tilde,sample)
	st_view(my_z1x,.,y1tilde,sample)
	st_view(md_z0x,.,d0tilde,sample)
	st_view(md_z1x,.,d1tilde,sample)
	st_view(mz_x,.,ztilde,sample)
	st_view(d,.,dvar,sample)
	st_view(y,.,yvar,sample)
	st_view(z,.,zvar,sample)

	n = rows(y)

	late =  mean( z :* (y :- my_z1x) :/ mz_x :-  ((1 :- z) :* (y :- my_z0x) :/ (1 :- mz_x)) :+ my_z1x :- my_z0x ) /
				mean( z :* (d :- md_z1x) :/ mz_x :-  ((1 :- z) :* (d :- md_z0x) :/ (1 :- mz_x)) :+ md_z1x :- md_z0x ) 
	late_se = variance(( z :* (y :- my_z1x) :/ mz_x :-  ((1 :- z) :* (y :- my_z0x) :/ (1 :- mz_x)) :+ my_z1x :- my_z0x ) :/ 
			   mean( z :* (d :- md_z1x) :/ mz_x :-  ((1 :- z) :* (d :- md_z0x) :/ (1 :- mz_x)) :+ md_z1x :- md_z0x )) / n

	st_matrix(outlate,late)
	st_matrix(outlatese,late_se)
}

end 

/*

ATE <- function(y, d, my_d1x, my_d0x, md_x)
{
  return( mean( (d * (y - my_d1x) / md_x) -  ((1 - d) * (y - my_d0x) / (1 - md_x)) + my_d1x - my_d0x ) );
}

SE.ATE <- function(y, d, my_d1x, my_d0x, md_x)
{
  return( sd( (d * (y - my_d1x) / md_x) -  ((1 - d) * (y - my_d0x) / (1 - md_x)) + my_d1x - my_d0x )/sqrt(length(y)) );
}

LATE <- function(y, d, z, my_z1x, my_z0x, mz_x, md_z1x, md_z0x)
{
  return( mean( z * (y - my_z1x) / mz_x -  ((1 - z) * (y - my_z0x) / (1 - mz_x)) + my_z1x - my_z0x ) / 
			mean( z * (d - md_z1x) / mz_x -  ((1 - z) * (d - md_z0x) / (1 - mz_x)) + md_z1x - md_z0x ) );
}

SE.LATE <- function(y, d, z, my_z1x, my_z0x, mz_x, md_z1x, md_z0x)
{
  return( sd(( z * (y - my_z1x) / mz_x -  ((1 - z) * (y - my_z0x) / (1 - mz_x)) + my_z1x - my_z0x ) / 
			   mean( z * (d - md_z1x) / mz_x -  ((1 - z) * (d - md_z0x) / (1 - mz_x)) + md_z1x - md_z0x )) / sqrt(length(y)) );
}
