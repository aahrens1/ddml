*** ddml estimation: interactive model

program _ddml_estimate_interactive, eclass sortpreserve

    syntax namelist(name=mname) [if] [in] , /// 
                                [  ///
                                ROBust ///
                                show(string) /// dertermines which to post
                                clear /// deletes all tilde-variables (to be implemented)
                                avplot ///
                                debug ///
                                * ]

    // base sample for estimation - determined by if/in
    marksample touse
    // also exclude obs already excluded by ddml sample
    qui replace `touse' = 0 if `mname'_sample==0

    if ("`show'"=="") {
        local show opt 
    }
    //mata: `mname'.nameDtilde
    //mata: st_local("Ztilde",invtokens(`mname'.nameZtilde))
    mata: st_local("Dtilde",invtokens(`mname'.nameDtilde))
    mata: st_local("Ytilde",invtokens(`mname'.nameYtilde))
    mata: st_local("nameD",invtokens(`mname'.nameD))
    mata: st_local("nameY",invtokens(`mname'.nameY))
    mata: st_local("Y0opt",`mname'.nameY0opt)
    mata: st_local("Y1opt",`mname'.nameY1opt)
    mata: st_local("Dopt",`mname'.nameDopt)
    //mata: st_local("Zopt",`mname'.nameZopt)

    if ("`debug'"!="") {
        di "`Ytilde'"
        di "`Ztilde'"
        di "`Dtilde'"
    }

    _ddml_allcombos `Ytilde' - `Ytilde' - `Dtilde', ///
                                                        putlast(`Y0opt' `Y1opt' `Dopt') ///
                                                        addprefix("")
                                                        `debug'  
    local ncombos = r(ncombos)
    local tokenlen = `ncombos'*2 -1
    local y0list `r(colstr1)'
    local y1list `r(colstr2)'
    local Dlist `r(colstr3)'

    local j = 1
    forvalues i = 1(2)`tokenlen' {
        if ("`show'"=="all"|`i'==`tokenlen') {
            tokenize `y0list' , parse("-")
            local y0 ``i''
            tokenize `y1list' , parse("-")
            local y1 ``i''
            tokenize `Dlist' , parse("-")
            local d ``i''
            if (`j'==`ncombos') {
                if "`show'"=="all" di as res "Optimal model: " _c
                local qui qui
            } 
            else {
                local qui
            }
            di as res "DML with Y0=`y0', Y1=`y1' and D=`d':"
            `qui' _ddml_ate, yvar(`nameY') dvar(`nameD') y0tilde(`y0') y1tilde(`y1') dtilde(`d') touse(`touse')
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
        ereturn local vcetype   robust
    }
    ereturn display

end
