*** ddml estimation: partial linear model
program _ddml_estimate_partial, eclass sortpreserve

	syntax namelist(name=mname) [if] [in] , /// 
								[  ///
								ROBust ///
								show(string) /// dertermines which to post
								clear /// deletes all tilde-variables (to be implemented)
								avplot ///
								* ]


	// loop through all Ytildes and all possible combinations of Dtildes
	mata: st_local("ylist",invtokens(`mname'.nameYtilde))
	// allcombos program put combinations in r(oplists) separated by commas
	allcombos `mname'.eqnlistD
	local xlists `r(oplists)'
	di "`xlists'"

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
	   twoway (scatter `Yopt' `Dopt') (lfit `Yopt' `Dopt')
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
