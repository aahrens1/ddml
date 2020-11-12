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
	
	// ylist is not separated by commas so just count the words
	local ycount : word count `ylist'
	// xlist is separated by commas unless there's only one list
	local xlists : subinstr local xlists "," ",", all count(local xcount)
	// multi=0 if there's only one estimation
	local multi = `ycount' - 1 + `xcount'
	if `multi' {
		foreach yvar of varlist `ylist' {
			tokenize "`xlists'", parse(",")
			local i 1
			while "``i''" ~= "" {
				di
				qui reg `yvar' ``i'' , nocons `robust' noheader
				di as res "DML with Y~ `yvar' and D~ ``i'' (N=`e(N)'):"
				// replay
				reg, noheader
				// since commas are in local2, increment by 2
				local i = `i'+2
			}
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
	local N = e(N)
	ereturn clear
	ereturn post `b' `V', depname(`Yopt') obs(`N')
	if "`robust'"~="" {
		ereturn local vcetype	robust
	}
	di
	if `multi' > 0 {
		di as res "Optimal model: DML with optimal Y~ `Yopt' and optimal D~ `Dopt' (N=`N'):"
	}
	else {
		di as res "DML with Y~ `Yopt' and D~ `Dopt' (N=`N'):"
	}
	ereturn display

end
