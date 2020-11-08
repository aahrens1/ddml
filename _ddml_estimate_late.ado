*** ddml estimation: LATE model
program _ddml_estimate_late, eclass sortpreserve

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
    local y0optid = e(y0optid)
    local y1optid = e(y1optid)
    local d0optid = e(d0optid)
    local d1optid = e(d1optid)
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
        forvalues i0 = 1(1)`yestn' {
            forvalues i1 = 1(1)`yestn' {
                forvalues j0 = 1(1)`destn' {
                    forvalues j1 = 1(1)`destn' {
                        forvalues l = 1(1)`zestn' {
                            if (`i0'==`y0optid' & `i1'==`y1optid' & `j0'==`d0optid' & `j1'==`d1optid' & `l'==`zoptid') {
                                // do nothing: optimal model always comes last 
                                // and is estimated below
                                di "" _c
                            }
                            else {
                                di as text "DML with `ycmd`i'' (`yname`i''), `dcmd`j'' (`dname`j'') and instrument `zcmd`j'' (`zname`j''):"
                                tempname b
                                tempname V 
                                mata: LATE("`yvar'","`dvar'","`zvar'","`ytilde`i0''", "`ytilde`i1''", "`dtilde`j0''","`dtilde`j1''","`ztilde`l''","`touse'","`b'","`V'")
                                matrix colnames `b' = "`dvar'"
                                matrix rownames `b' = "`yvar'"
                                matrix colnames `V' = "`dvar'"
                                matrix rownames `V' = "`dvar'"
                                // display
                                ereturn clear
                                ereturn post `b' `V' 
                                ereturn display
                            }
                        }
                    }
                }
            }
        }
    }
    *** estimate best model
    di as res "Optimal model: DML with `ycmd`yoptid'' (`yname`yoptid''), `dcmd`doptid'' (`dname`doptid'') and instrument `zcmd`zoptid'' (`zname`zoptid''):"
	tempname b
	tempname V 
    mata: LATE("`yvar'","`dvar'","`zvar'","`ytilde`y0optid''", "`ytilde`y1optid''", "`dtilde`d0optid''","`dtilde`d1optid''","`ztilde`zoptid''","`touse'","`b'","`V'")
	matrix colnames `b' = "`dvar'"
	matrix rownames `b' = "`yvar'"
 	matrix colnames `V' = "`dvar'"
	matrix rownames `V' = "`dvar'"
    // display
	ereturn clear
	ereturn post `b' `V' 
	ereturn display

end

mata: 

void LATE(  string scalar yvar,	     // Y
            string scalar dvar,      // D
            string scalar zvar,      // Z
            string scalar y0tilde,   // E[Y|X,Z=0]
            string scalar y1tilde,   // E[Y|X,Z=1]
            string scalar d0tilde,   // E[D|X,Z=0]
            string scalar d1tilde,   // E[D|X,Z=1]
            string scalar ztilde,    // E[Z|X]
            string scalar sample,    // sample
            string scalar outlate,   // output: name of matrix to store b
            string scalar outlatese  // output: name of matrix to store V
            )
{
    st_view(my_z0x,.,y0tilde,sample)
    st_view(my_z1x,.,y1tilde,sample)
    st_view(md_z0x,.,d0tilde,sample)
    st_view(md_z1x,.,d1tilde,sample)
    st_view(mz_x,.,ztilde,sample)
    st_view(d,.,dvar,sample)
    st_view(y,.,yvar,sample)
    st_view(z,.,zvar,sample)

    n = rows(y)

    late =  mean( z :* (y :- my_z1x) :/ mz_x :-  ((1 :- z) :* (y :- my_z0x) :/ (1 :- mz_x)) :+ my_z1x :- my_z0x ) /
                mean( z :* (d :- md_z1x) :/ mz_x :-  ((1 :- z) :* (d :- md_z0x) :/ (1 :- mz_x)) :+ md_z1x :- md_z0x ) 
    late_se = variance(( z :* (y :- my_z1x) :/ mz_x :-  ((1 :- z) :* (y :- my_z0x) :/ (1 :- mz_x)) :+ my_z1x :- my_z0x ) :/ 
               mean( z :* (d :- md_z1x) :/ mz_x :-  ((1 :- z) :* (d :- md_z0x) :/ (1 :- mz_x)) :+ md_z1x :- md_z0x )) / n

    st_matrix(outlate,late)
    st_matrix(outlatese,late_se)
}

end 

