*** ddml estimation: LATE model

program _ddml_estimate_late, eclass sortpreserve

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
    mata: st_local("nameD",invtokens(`mname'.nameD))
    mata: st_local("nameY",invtokens(`mname'.nameY))
    mata: st_local("nameZ",invtokens(`mname'.nameZ))
    mata: st_local("Y0opt",`mname'.nameY0opt)
    mata: st_local("Y1opt",`mname'.nameY1opt)
    mata: st_local("D0opt",`mname'.nameD0opt)
    mata: st_local("D1opt",`mname'.nameD1opt)
    mata: st_local("Zopt",`mname'.nameZopt)

    if ("`debug'"!="") {
        di "`Ytilde'"
        di "`Ztilde'"
        di "`Dtilde'"
    }

    _ddml_allcombos `Ytilde' - `Ytilde' - `Dtilde' - `Dtilde' - `Ztilde' , ///
                                                        putlast(`Y0opt' `Y1opt' `D0opt' `D1opt' `Zopt') ///
                                                        `debug'  ///
                                                        addprefix("`mname'_")
    return list
    local ncombos = r(ncombos)
    local tokenlen = `ncombos'*2 -1
    local y0list `r(colstr1)'
    local y1list `r(colstr2)'
    local d0list `r(colstr3)'
    local d1list `r(colstr4)'
    local Zlist `r(colstr5)' 

    if ("`show'"=="all") {
        local j = 1
        forvalues i = 1(2)`tokenlen' {
            tokenize `y0list' , parse("-")
            local y0 ``i''
            tokenize `y1list' , parse("-")
            local y1 ``i''
            tokenize `d0list' , parse("-")
            local d0 ``i''
            tokenize `d1list' , parse("-")
            local d1 ``i''
            tokenize `Zlist' , parse("-")
            local z ``i''
            if (`j'==`ncombos') {
                di as res "Optimal model: " _c
            }
            di as res "DML with Y0=`y0', Y1=`y1' and D=`d':"
            _ddml_late, yvar(`nameY') y0tilde(`y0') y1tilde(`y1') ///
                        dvar(`nameD') d0tilde(`d0') d1tilde(`d1') ///
                        zvar(`nameZ') ztilde(`z')  

            local j= `j'+1
         }
    }

    if ("`show'"=="opt") {
        *** estimate best model
        di as res "Optimal model: DML with Y=`Yopt' and D=`Dopt'"
        qui _ddml_late, yvar(`nameY') y0tilde(`Y0opt') y1tilde(`Y1opt') ///
                        dvar(`nameD') d0tilde(`D0opt') d1tilde(`D1opt') ///
                        zvar(`nameZ') ztilde(`Zopt')   
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
        ereturn local vcetype   robust
    }
    ereturn display

end
 