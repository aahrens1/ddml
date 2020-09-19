*! ddml2 v0.1.2 (19 sep 2020)

program ddml, eclass

    version 13

    syntax [anything] , [*]

    local subcmd: word 1 of `anything'
    
    *** initialize new estimation
    if "`subcmd'"=="init" {
        local model: word 2 of `anything'
        if ("`model'" != "partial") {
            di as err "no or wrong model specified: currently only 'ddml init partial' supported"
            exit 1
        }
        ereturn clear
        ereturn local cmd "ddml_init"
        ereturn local model "`model'"
        ereturn scalar yest = 0 // number of y-equation estimations
        ereturn scalar dest = 0 // number of d-equation estimations
    } 
    else {
        if ("`e(cmd)'"!="ddml_init") {
            di as err "no ddml object found. first initialize using 'ddml init'"
            exit 1
        }
    }

    *** add y-equation estimation
    if "`subcmd'"=="yeq" {
        local estname: word 2 of `anything'
        gettoken left theeq: anything, parse(":")
        gettoken left theeq: theeq, parse(":") 
        gettoken theeq: theeq, parse("(") match(paren)
        ereturn local y`=`e(yest)'+1' "y_`estname'"
        ereturn local ycmd`=`e(yest)'+1' "`theeq'"
        ereturn scalar yest = `e(yest)'+1
    }

    *** add d-equation estimation
    if "`subcmd'"=="deq" {
        local estname: word 2 of `anything'
        gettoken left theeq: anything, parse(":")
        gettoken left theeq: theeq, parse(":") 
        gettoken theeq: theeq, parse("(") match(paren)
        ereturn local d`=`e(dest)'+1' "y_`estname'"
        ereturn local dcmd`=`e(dest)'+1' "`theeq'"
        ereturn scalar dest = `e(dest)'+1
    }

    *** estimate
    if "`subcmd'" =="estimate" {
        _ddml_estimate, `options'
    }
end


*** ddml estimation command
program ddml2, eclass sortpreserve


end