*** ddml estimation: interactive model
program _ddml_estimate_interactive, eclass sortpreserve

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
    local y0optid = e(y0optid)
    local y1optid = e(y1optid)
    local doptid = e(doptid)
    local yvar = e(depvar)
    local dvar = e(dvar)
    
    *** retrieve variable names
    local alltilde
    forvalues i = 1(1)`yestn' {
        local ytilde`i' `e(y`i')'
        local yname`i' `e(y`i')'
        local ycmd`i' `e(ycmd`i')'
        local alltilde `alltilde' `ytilde`i''
        di "`ytilde`i''"
    }
    forvalues i = 1(1)`destn' {
        local dtilde`i' `e(d`i')'
        local dname`i' `e(d`i')'
        local dcmd`i' `e(dcmd`i')'
        local alltilde `alltilde' `dtilde`i''
    }

    *** mark sample
    marksample touse 
    markout `touse' `yvar' `dvar' `alltilde'

    *** do estimation
    if ("`show'"=="all") {
        forvalues i0 = 1(1)`yestn' {
            forvalues i1 = 1(1)`yestn' {
                forvalues j = 1(1)`destn' {
                    if (`i0'==`y0optid' & `i1'==`y1optid' & `j'==`doptid') {
                        // do nothing: optimal model always comes last 
                        // and is estimated below
                        di "" _c
                    }
                    else {
                        di as res "Optimal model: DML with `ycmd`i0'' (`yname`i0''), `ycmd`i1'' (`yname`i1'') and `dcmd`j'' (`dname`j''):"
                        tempname b
                        tempname V 
                        mata: ATE("`yvar'","`dvar'","`ytilde`i0''", "`ytilde`i1''", "`dtilde`j''","`touse'","`b'","`V'")

                        // display
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
    di as res "Optimal model: DML with `ycmd`y0optid'' (`yname`y0optid''), `ycmd`y1optid'' (`yname`y0optid'') and `dcmd`doptid'' (`dname`doptid''):"
    tempname b
	tempname V 
    mata: ATE("`yvar'","`dvar'","`ytilde`y0optid''", "`ytilde`y1optid''", "`dtilde`doptid''","`touse'","`b'","`V'")

    // display
	matrix colnames `b' = "`dvar'"
	matrix rownames `b' = "`yvar'"
 	matrix colnames `V' = "`dvar'"
	matrix rownames `V' = "`dvar'"
	ereturn clear
	ereturn post `b' `V' 
	ereturn display

end



********************************************************************************
*** Mata section															 ***
********************************************************************************

mata:

void ATE(   string scalar yvar,	    // Y
            string scalar dvar,     // D
            string scalar y0tilde,  // E[Y|X,D=0]
            string scalar y1tilde,  // E[Y|X,D=1]
            string scalar dtilde,   // E[D|X]
            string scalar sample,   // sample
            string scalar outate,   // output: name of matrix to store b
            string scalar outatese  // output: name of matrix to store V
            )
{
    st_view(my_d0x,.,y0tilde,sample)
    st_view(my_d1x,.,y1tilde,sample)
    st_view(md_x,.,dtilde,sample)
    st_view(d,.,dvar,sample)
    st_view(y,.,yvar,sample)

    n = rows(y)

    my_d0x = my_d0x :* (1:-d)
    my_d1x = my_d1x :* d

    te  = (d :* (y :- my_d1x) :/ md_x) :-  ((1 :- d) :* (y :- my_d0x) :/ (1 :- md_x)) :+ my_d1x :- my_d0x  
    ate = mean(te)
    ate_V =  variance(te)/n

    st_matrix(outate,ate)
    st_matrix(outatese,ate_V)
}

end

/*

ATE <- function(y, d, my_d1x, my_d0x, md_x)
{
  return( mean( (d * (y - my_d1x) / md_x) -  ((1 - d) * (y - my_d0x) / (1 - md_x)) + my_d1x - my_d0x ) );
}

SE.ATE <- function(y, d, my_d1x, my_d0x, md_x)
{
  return( sd( (d * (y - my_d1x) / md_x) -  ((1 - d) * (y - my_d0x) / (1 - md_x)) + my_d1x - my_d0x )/sqrt(length(y)) );
}

