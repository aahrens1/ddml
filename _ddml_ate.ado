
program _ddml_ate, eclass

	version 13

	syntax [anything] , yvar(varname) ///
				dvar(varname) ///
				y0tilde(varname) ///
				y1tilde(varname) ///
				dtilde(varname) ///
				[touse(varname)]

	tempname b
    tempname V 

    if "`touse'"=="" {
    	tempvar touse
    	mark `touse'
    	markout `touse' `yvar' `dvar' `y0tilde' `y1tilde' `dtilde'
    }

    mata: ATE("`yvar'","`dvar'","`y0tilde'", "`y1tilde'", "`dtilde'","`touse'","`b'","`V'")

    // display
    local N = `r(N)'
    matrix colnames `b' = "`dvar'"
    matrix rownames `b' = "`yvar'"
    matrix colnames `V' = "`dvar'"
    matrix rownames `V' = "`dvar'"
    ereturn clear
    ereturn post `b' `V', depname(`yvar') obs(`N')

end


********************************************************************************
*** Mata section                                                             ***
********************************************************************************

mata:

void ATE(   string scalar yvar,     // Y
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

    //my_d0x = my_d0x :* (1:-d)
    //my_d1x = my_d1x :* d

    te  = (d :* (y :- my_d1x) :/ md_x) :-  ((1 :- d) :* (y :- my_d0x) :/ (1 :- md_x)) :+ my_d1x :- my_d0x  
    ate = mean(te)
    ate_V =  variance(te)/n

    st_numscalar("r(N)",n)
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

