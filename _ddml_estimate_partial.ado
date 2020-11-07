*** ddml estimation: partial linear model
program _ddml_estimate_partial, eclass sortpreserve

	syntax namelist(name=mname) [if] [in] , /// 
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

	// loop through all Ytildes and all possible combinations of Dtildes
	mata: st_local("ylist",invtokens(`mname'.nameYtilde))
	// allcombos program put combinations in r(oplists) separated by commas
	allcombos `mname'.eqnlistD
	local xlists `r(oplists)'

	foreach yvar of varlist `ylist' {
		tokenize "`xlists'", parse(",")
		local i 1
		while "``i''" ~= "" {
			di
			di as res "DML with Y~ `yvar' and D~ ``i'':"
			reg `yvar' ``i'' , nocons `robust' noheader
			// since commas are in local2, increment by 2
			local i = `i'+2
		}
	}
	
   	mata: st_local("Yopt",`mname'.nameYopt)
   	mata: st_local("Dopt",`mname'.nameDopt)

	*** estimate best model
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
	ereturn post `b' `V', depname(`Yopt')
	if "`robust'"~="" {
		ereturn local vcetype	robust
	}
	di
	di as res "Optimal model: DML with optimal Y~ `Yopt' and optimal D~ `Dopt':"
	ereturn display

end
