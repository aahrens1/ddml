*! ddml v1.1
*! last edited: 28 dec 2022
*! authors: aa/ms

program define _ddml_nnls_p
    version 14
 
    syntax newvarname [if] [in] , [ xb r]
 
    marksample touse, novarlist
 
    local nopts : word count `xb' `r'
    if `nopts' >1 {
        display "{err}only one statistic may be specified"
        exit 498
    }
 
    if `nopts' == 0 {
        local xb xb
        display "expected xb"
    }   
 
    if "`xb'" != "" {
        _predict `typlist' `varlist' if `touse' , xb
    }
    else {
        tempvar xbv
        quietly _predict double `xbv' if `touse' , xb
        generate `typlist' `varlist' = `e(depvar)' - `xbv' if `touse'
    }
    
end
