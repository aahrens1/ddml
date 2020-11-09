*** ddml estimation: partial linear model
program _ddml_estimate_iv, eclass sortpreserve

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
	local Dlists `r(oplists)'
	di "`Dlists'"
    allcombos `mname'.eqnlistZ
	local Zlists `r(oplists)'
	di "`Zlists'"

	mata: st_local("Yopt",`mname'.nameYopt)
   	mata: st_local("Dopt",`mname'.nameDopt)
    mata: st_local("Zopt",`mname'.nameZopt)

	foreach yvar of varlist `ylist' {
		local i = 1
		local d `Dlists'
		while "`d'" ~= "" {
			tokenize "`Dlists'", parse(",")
			local d ``i''
			local j = 1
			local z `Zlists'
            while "`z'" ~= "" {
				tokenize "`Zlists'", parse(",")
				local z ``j''
				if ("`yvar'"=="`Yopt'"&"`d'"=="`Dopt'"&"`z'"=="`Zopt'") {
					// do nothing; optimal model comes last
					di "" _c
				}
				else {
					di
					di as res "DML with Y=`yvar' and D=`d', Z=`z':"
					ivreg2 `yvar' (`d'=`z') , nocons `robust' noheader nofooter
				}
                // since commas are in local2, increment by 2
                local j = `j' + 2
				local z ``j''
            }
            // since commas are in local2, increment by 2
            local i = `i' + 2
			local d ``i''
		}
	}
	
	*** estimate best model
    di as res "DML with Y=`Yopt' and D=`Dopt', Z=`Zopt':"
    qui ivreg2 `Yopt' (`Dopt'=`Zopt') , nocons `robust' noheader nofooter
	
    // plot
	//if ("`avplot'"!="") {
    //   // only works with one Dopt
	//   twoway (scatter `Yopt' `Dopt') (lfit `Yopt' `Dopt')
	//}

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
	di as res "Optimal model: DML with optimal Y=`Yopt' and optimal D= `Dopt':"
	ereturn display

end
