*** ddml estimation: partial linear IV model

program _ddml_estimate_optimaliv, eclass sortpreserve

	syntax namelist(name=mname) [if] [in] , /// 
								[  ///
								ROBust ///
								show(string) /// dertermines which to post
								clear /// deletes all tilde-variables (to be implemented)
								avplot ///
								debug ///
								* ]

	if ("`show'"=="") {
		local show opt 
	}

	// base sample for estimation - determined by if/in
	marksample touse
	// also exclude obs already excluded by ddml sample
	qui replace `touse' = 0 if `mname'_sample==0

    //mata: `mname'.nameDtilde
    mata: st_local("Yopt",`mname'.nameYopt)
    mata: st_local("Dopt",invtokens(`mname'.nameDopt))
    mata: st_local("DHopt",invtokens(`mname'.nameDHopt))
    mata: st_local("nameD",invtokens(`mname'.nameD))
    mata: st_local("nameY",invtokens(`mname'.nameY))

    if ("`debug'"!="") {
        di "`Ytilde'"
        di "`DHtilde'"
        di "`Dtilde'"
    }

    _ddml_make_varlists, mname(`mname')
    _ddml_allcombos `r(eq)' , putlast(`Yopt' `Dopt' `DHopt') ///
                                                `debug' ///
                                                dpos_end(`r(dpos_end)') ///
                                                zpos_start(`r(zpos_start)') zpos_end(`r(zpos_end)') ///
                                                addprefix("`mname'_")

	local ncombos = r(ncombos)
	local tokenlen = `ncombos'*2 -1
	local ylist `r(ystr)'
	local Dlist `r(dstr)'
	local DHlist `r(zstr)' 

	local j = 1
	forvalues i = 1(2)`tokenlen' {
		if ("`show'"=="all"|`i'==`tokenlen') {
		   	tokenize `ylist' , parse("-")
		    local y ``i''
		    tokenize `Dlist' , parse("-")
		    local d ``i''
		    tokenize `DHlist' , parse("-")
		    local dh ``i''
		    if (`j'==`ncombos') {
		        if "`show'"=="all" di as res "Optimal model: " _c
		    	local qui qui
		    } 
		    else {
		    	local qui
		    }
		    di as res "DML with E[Y|X]=`y' and E[D|X]=`d', E[D|X,Y]=`dh':"
		     `qui' _ddml_optiv, yvar(`nameY') dvar(`nameD') dhtilde(`dh') ytilde(`y') dtilde(`d') touse(`touse')  `debug'
		}
	    local j= `j'+1
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
	ereturn post `b' `V', depname(`Yopt') obs(`N') esample(`touse')
	if "`robust'"~="" {
		ereturn local vcetype	robust
	}
    ereturn display

end