clear all
 
if ("`c(username)'"=="kahrens") {
	adopath + "/Users/kahrens/MyProjects/ddml"
	adopath + "/Users/kahrens/MyProjects/pystacked"
	cd "/Users/kahrens/MyProjects/ddml/cert"
}

//webuse cattaneo2, clear
local cattaneo2 cattaneo2.dta
//save `cattaneo2'

//use "http://fmwww.bc.edu/repec/bocode/j/jtpa.dta",clear   
local jtpa jtpa.dta
//save `jtpa'

//use "https://github.com/VC2015/DMLonGitHub/raw/master/sipp1991.dta"
local sipp1991 sipp1991.dta
//save `sipp1991'


frame create Rresults
frame change Rresults
insheet using DoubleML_results.csv, comma  
list 
rename coef Rcoef
rename se Rse
gen Scoef=.
gen Sse=.
frame change default 

which ddml

global tol = 10e-4

local i = 0

********************************************************************************
*** catteneo2 data: interactive model										 ***
********************************************************************************

use `cattaneo2', clear
 
global Y bweight
global D mbsmoke
global X mage prenatal1 mmarried fbaby mage medu

gen fid=_n>2321

ddml init interactive, foldvar(fid)
ddml E[Y|X,D], gen(regy): reg $Y $X
ddml E[D|X], gen(regd): logit $D $X
ddml crossfit

ddml estimate  ,  trim(0)  

local b = _b[mbsmoke]
local se = _se[mbsmoke] 
 
frame change Rresults
replace Scoef = `b' if application == "cattaneo2" & model =="interactive ATE"
replace Sse = `se' if application == "cattaneo2" & model =="interactive ATE"
frame change default

ddml estimate   , atet trim(0)

local b = _b[mbsmoke]
local se = _se[mbsmoke] 
 
frame change Rresults
replace Scoef = `b' if application == "cattaneo2" & model =="interactive ATTE"
replace Sse = `se' if application == "cattaneo2" & model =="interactive ATTE"
frame change default

********************************************************************************
*** catteneo2 data: partial linear model									 ***
********************************************************************************

use `cattaneo2', clear
 
global Y bweight
global D mbsmoke
global X mage prenatal1 mmarried fbaby mage medu

gen fid=_n>2321

ddml init partial, foldvar(fid)
ddml E[Y|X], gen(regy): reg $Y $X
ddml E[D|X], gen(regd): reg $D $X
ddml crossfit

ddml estimate, vce(hc2)

local b = _b[mbsmoke]
local se = _se[mbsmoke] 

frame change Rresults
replace Scoef = `b' if application == "cattaneo2" & model =="partial"
replace Sse = `se' if application == "cattaneo2" & model =="partial"
frame change default

********************************************************************************
*** jtpa: LATE 			 													 ***
********************************************************************************
 
use `jtpa',clear   

gen fid=_n>5602

global Y earnings
global D training
global Z assignmt 
global X sex age married black hispanic 

ddml init late, foldvar(fid)
ddml E[Y|X,Z]: reg $Y $X 
ddml E[D|X,Z]: logit $D $X
ddml E[Z|X]: logit $Z $X
ddml crossfit 
ddml estimate ,  trim(0)

local b = _b[$D]
local se = _se[$D] 
 
frame change Rresults
replace Scoef = `b' if application == "jtpa" & model =="interactive IV"
replace Sse = `se' if application == "jtpa" & model =="interactive IV"
frame change default

********************************************************************************
*** jtpa: partial IV.  	 													 ***
********************************************************************************
 
use `jtpa',clear   

gen fid=_n>5602

global Y earnings
global D training
global Z assignmt 
global X sex age married black hispanic 

ddml init iv, foldvar(fid)
ddml E[Y|X]: reg $Y $X 
ddml E[D|X]: reg $D $X
ddml E[Z|X]: reg $Z $X
ddml crossfit 
ddml estimate , trim(0)

local b = _b[$D]
local se = _se[$D] 
 
frame change Rresults
replace Scoef = `b' if application == "jtpa" & model =="partial IV"
replace Sse = `se' if application == "jtpa" & model =="partial IV"
frame change default

********************************************************************************
*** jtpa: interactive	 													 ***
********************************************************************************
 
use `jtpa',clear   

gen fid=_n>5602

global Y earnings
global D assignmt 
global X sex age married black hispanic 

ddml init interactive, foldvar(fid)
ddml E[Y|X,D]: reg $Y $X 
ddml E[D|X]: logit $D $X
ddml crossfit 
ddml estimate ,  trim(0)

local b = _b[$D]
local se = _se[$D] 
 
frame change Rresults
replace Scoef = `b' if application == "jtpa" & model =="interactive ATE"
replace Sse = `se' if application == "jtpa" & model =="interactive ATE"
frame change default

ddml estimate , atet trim(0)

local b = _b[$D]
local se = _se[$D] 
 
frame change Rresults
replace Scoef = `b' if application == "jtpa" & model =="interactive ATTE"
replace Sse = `se' if application == "jtpa" & model =="interactive ATTE"
frame change default

********************************************************************************
*** jtpa: partial		 													 ***
********************************************************************************
 
use `jtpa',clear   

gen fid=_n>5602

global Y earnings
global D assignmt 
global X sex age married black hispanic 

ddml init partial, foldvar(fid)
ddml E[Y|X]: reg $Y $X 
ddml E[D|X]: reg $D $X
ddml crossfit 
ddml estimate ,  vce(hc3)

local b = _b[$D]
local se = _se[$D] 
 
frame change Rresults
replace Scoef = `b' if application == "jtpa" & model =="partial"
replace Sse = `se' if application == "jtpa" & model =="partial"
frame change default


********************************************************************************
*** 401k: partial		 													 ***
********************************************************************************
 
clear
use `sipp1991'

global Y net_tfa
global D e401
global X tw   age   inc fsize  educ    db  marr twoearn  pira  hown

gen fid = mod(_n,3)+1

ddml init partial, foldvar(fid)
ddml E[Y|X]: reg $Y $X 
ddml E[D|X]: reg $D $X
ddml crossfit 
ddml estimate , vce(hc3) 

local b = _b[$D]
local se = _se[$D] 
 
frame change Rresults
replace Scoef = `b' if application == "401k" & model =="partial"
replace Sse = `se' if application == "401k" & model =="partial"
frame change default

********************************************************************************
*** 401k: interactive	 													 ***
********************************************************************************
 
clear
use `sipp1991'

global Y net_tfa
global D e401
global X tw   age   inc fsize  educ    db  marr twoearn  pira  hown

gen fid = mod(_n,3)+1

ddml init interactive, foldvar(fid)
ddml E[Y|X,D]: reg $Y $X 
ddml E[D|X]: logit $D $X
ddml crossfit 
ddml estimate  , trim(0)

local b = _b[$D]
local se = _se[$D] 
 
frame change Rresults
replace Scoef = `b' if application == "401k" & model =="interactive ATE"
replace Sse = `se' if application == "401k" & model =="interactive ATE"
frame change default

ddml estimate , atet trim(0)

local b = _b[$D]
local se = _se[$D] 
 
frame change Rresults
replace Scoef = `b' if application == "401k" & model =="interactive ATTE"
replace Sse = `se' if application == "401k" & model =="interactive ATTE"
frame change default

********************************************************************************
*** 401k: iv			 													 ***
********************************************************************************
 
clear
use `sipp1991'

global Y net_tfa
global D p401
global Z e401
global X tw   age   inc fsize  educ    db  marr twoearn  pira  hown

gen fid = mod(_n,3)+1

ddml init iv, foldvar(fid)
ddml E[Y|X]: reg $Y $X 
ddml E[Z|X]: reg $Z $X
ddml E[D|X]: reg $D $X
ddml crossfit 
ddml estimate , robust

local b = _b[$D]
local se = _se[$D] 
 
frame change Rresults
replace Scoef = `b' if application == "401k" & model =="partial IV"
replace Sse = `se' if application == "401k" & model =="partial IV"
frame change default


********************************************************************************
*** 401k: late			 													 ***
********************************************************************************
 
clear
use `sipp1991'

global Y net_tfa
global D p401
global Z e401
global X tw   age   inc fsize  educ    db  marr twoearn  pira  hown

gen fid = mod(_n,3)+1

ddml init late, foldvar(fid)  
ddml E[Y|X,Z]: reg $Y $X 
ddml E[D|X,Z]: logit $D $X
ddml E[Z|X]: logit $Z $X

ddml crossfit 
ddml estimate ,  trim(0)

local b = _b[$D]
local se = _se[$D] 
 
frame change Rresults
replace Scoef = `b' if application == "401k" & model =="interactive IV"
replace Sse = `se' if application == "401k" & model =="interactive IV"
frame change default


********************************************************************************
*** compare			 													     ***
********************************************************************************
 
frame change Rresults

gen Dcoef = reldif(Rcoef,Scoef) 
gen Dse = reldif(Rse,Sse) 
list

assert reldif(Rcoef,Scoef) <$tol
assert  reldif(Rse,Sse) <$tol


 