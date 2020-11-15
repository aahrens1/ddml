*** ddml estimation: partial linear IV model

program _ddml_estimate_iv, eclass sortpreserve

	syntax namelist(name=mname) [if] [in] , /// 
								[  ///
								ROBust ///
								show(string) /// dertermines which to post
								clear /// deletes all tilde-variables (to be implemented)
								avplot ///
								debug ///
								* ]

	if ("`show'"=="") {
		local show all 
	}

    //mata: `mname'.nameDtilde
    mata: st_local("Ztilde",invtokens(`mname'.nameZtilde))
    mata: st_local("Dtilde",invtokens(`mname'.nameDtilde))
    mata: st_local("Ytilde",invtokens(`mname'.nameYtilde))
    mata: st_local("Yopt",`mname'.nameYopt)
    mata: st_local("Dopt",`mname'.nameDopt)
    mata: st_local("Zopt",`mname'.nameZopt)

    if ("`debug'"!="") {
        di "`Ytilde'"
        di "`Ztilde'"
        di "`Dtilde'"
    }

    _ddml_allcombos `Ytilde' - `Dtilde' - `Ztilde' , putlast(`Yopt' `Dopt' `Zopt') ///
                                                        `debug' ///
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
	    	tokenize `ylist' , parse("-")
	    	local y ``i''
	    	tokenize `Dlist' , parse("-")
	    	local d ``i''
	    	tokenize `Zlist' , parse("-")
	    	local z ``i''
	    	if (`j'==`ncombos') {
	    		di as res "Optimal model: " _c
	    	}
	    	di as res "DML with Y=`y' and D=`d', Z=`z':"
	       	ivreg2 `y' (`d'=`z') , nocons `robust' noheader nofooter

	        local j= `j'+1
	     }
	}

	if ("`show'"=="opt") {
		*** estimate best model
    	di as res "Optimal model: DML with Y=`Yopt' and D=`Dopt', Z=`Zopt':"
    	qui ivreg2 `Yopt' (`Dopt'=`Zopt') , nocons `robust' noheader nofooter
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
	ereturn display

end

/*
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
