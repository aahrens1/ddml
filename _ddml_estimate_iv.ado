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
	//allcombos `mname'.eqnlistD
	//local Dlists `r(oplists)'
	//di "`Dlists'"
    //allcombos `mname'.eqnlistZ
	//local Zlists `r(oplists)'
	//di "`Zlists'"

	mata: st_local("Ztilde",mat_to_varlist(`mname'.nameZtilde))
	mata: st_local("Dtilde",mat_to_varlist(`mname'.nameDtilde))
	mata: st_local("Ytilde",invtokens(`mname'.nameYtilde))
	mata: st_local("Yopt",`mname'.nameYopt)
   	mata: st_local("Dopt",`mname'.nameDopt)
    mata: st_local("Zopt",`mname'.nameZopt)
    di "`Ztilde'"
    di "`Dtilde'"
    di "`Ytilde'"

    _ddml_allcombos `Ytilde' | `Dtilde' | `Ztilde' , putlast(`Yopt' `Dopt' `Zopt') ///
    													debug  ///
    													dpos_start(2) dpos_end(2) ///
    													zpos_start(3) zpos_end(3)
	return list
	local ncombos = r(ncombos)
	local tokenlen = `ncombos'*2 -1
	local ylist `r(ystr)'
	local Dlist `r(dstr)'
	local Zlist `r(zstr)' 

	if ("`show'"=="all") {
	    local j = 1
	    di `tokenlen'
	    forvalues i = 1(2)`tokenlen' {
	    	tokenize `ylist'
	    	local y ``i''
	    	tokenize `Dlist' , parse("|")
	    	local d ``i''
	    	tokenize `Zlist' , parse("|")
	    	local z ``i''
	    	if (`j'==`ncombos') {
	    		di as res "Optimal model: DML with Y=`y' and D=`d', Z=`z':"
	        	qui ivreg2 `y' (`d'=`z') , nocons `robust' noheader nofooter
	    	}
	    	else {
	    		di as res "DML with Y=`y' and D=`d', Z=`z':"
	       		ivreg2 `y' (`d'=`z') , nocons `robust' noheader nofooter
	    	}

	        local j= `j'+1
	     }
	}

	if ("`show'"=="opt") {
		*** estimate best model
    	di as res "DML with Y=`Yopt' and D=`Dopt', Z=`Zopt':"
    	qui ivreg2 `Yopt' (`Dopt'=`Zopt') , nocons `robust' noheader nofooter
	}

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

mata:

string scalar mat_to_varlist(string matrix inmat)
{
	r = rows(inmat)
	for (i=1;i<=r;i++) {

		if (i==1) {
			str = invtokens(inmat[i,]) 
		}
		else {
			str = str + " | " + invtokens(inmat[i,]) 
		}
	} 
	return(str)
}
end
