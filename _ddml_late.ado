
program _ddml_late, eclass

	version 13

	syntax [anything] , yvar(varname) ///
				dvar(varname) ///
				y0tilde(varname) ///
				y1tilde(varname) ///
				d0tilde(varname) ///
                d1tilde(varname) ///
                zvar(varname) ///
                ztilde(varname) ///
				[touse(varname)]

    if "`touse'"=="" {
        tempvar touse
        mark `touse'
        markout `touse' `yvar' `dvar' `y0tilde' `y1tilde' `d0tilde' `d1tilde' `ztilde'
    }

    tempname b
    tempname V 
    mata: LATE("`yvar'","`dvar'","`zvar'","`y0tilde'", "`y1tilde'", "`d0tilde'","`d1tilde'","`ztilde'","`touse'","`b'","`V'")
    matrix colnames `b' = "`dvar'"
    matrix rownames `b' = "`yvar'"
    matrix colnames `V' = "`dvar'"
    matrix rownames `V' = "`dvar'"


    // display
    local N = `r(N)'
    matrix colnames `b' = "`dvar'"
    matrix rownames `b' = "`yvar'"
    matrix colnames `V' = "`dvar'"
    matrix rownames `V' = "`dvar'"
    ereturn clear
    ereturn post `b' `V', depname(`yvar') obs(`N')
    ereturn display

end


********************************************************************************
*** Mata section                                                             ***
********************************************************************************

mata: 
void LATE(  string scalar yvar,      // Y
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
    st_numscalar("r(N)",n)
    st_matrix(outlate,late)
    st_matrix(outlatese,late_se)
}
end 

/*
LATE <- function(y, d, z, my_z1x, my_z0x, mz_x, md_z1x, md_z0x)
{
  return( mean( z * (y - my_z1x) / mz_x -  ((1 - z) * (y - my_z0x) / (1 - mz_x)) + my_z1x - my_z0x ) / 
            mean( z * (d - md_z1x) / mz_x -  ((1 - z) * (d - md_z0x) / (1 - mz_x)) + md_z1x - md_z0x ) );
}
SE.LATE <- function(y, d, z, my_z1x, my_z0x, mz_x, md_z1x, md_z0x)
{
  return( sd(( z * (y - my_z1x) / mz_x -  ((1 - z) * (y - my_z0x) / (1 - mz_x)) + my_z1x - my_z0x ) / 
               mean( z * (d - md_z1x) / mz_x -  ((1 - z) * (d - md_z0x) / (1 - mz_x)) + md_z1x - md_z0x )) / sqrt(length(y)) );
}