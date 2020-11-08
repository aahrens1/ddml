*** ddml estimation: partial linear model
program _ddml_estimate_iv, eclass sortpreserve

    syntax [anything] [if] [in] , /// 
								[  ///
                                robust ///
                                show(string) /// dertermines which to post
                                clear /// deletes all tilde-variables (to be implemented)
                                * ]

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
    forvalues i = 1(1)`zestn' {
        local ztilde`i' `e(z`i')'
        local zname`i' `e(z`i')'
        local zcmd`i' `e(zcmd`i')'
    }

    *** do estimation
    if ("`show'"=="all") {
        forvalues i = 1(1)`yestn' {
            forvalues j = 1(1)`destn' {
                    forvalues l = 1(1)`zestn' {
                        if (`i'==`yoptid' & `j'==`doptid' & `l'==`zoptid') {
                            // do nothing: optimal model always comes last 
                            // and is estimated below
                            di "" _c
                        }
                        else {
                            di as text "DML with `ycmd`i'' (`yname`i''), `dcmd`j'' (`dname`j'') and instrument `zcmd`j'' (`zname`j''):"
                            qui ivreg2 `ytilde`i'' (`dtilde`j''=`ztilde`l''), nocons `robust' noheader
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
    }
    *** estimate best model
    di as res "Optimal model: DML with `ycmd`yoptid'' (`yname`yoptid''), `dcmd`doptid'' (`dname`doptid'') and instrument `zcmd`zoptid'' (`zname`zoptid''):"
    qui ivreg2 `ytilde`yoptid'' (`dtilde`doptid''=`ztilde`zoptid''), nocons `robust' noheader

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

end
