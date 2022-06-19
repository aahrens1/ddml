clear all
 
if ("`c(username)'"=="kahrens") {
	adopath + "/Users/kahrens/MyProjects/ddml"
	adopath + "/Users/kahrens/MyProjects/pystacked"
	cd "/Users/kahrens/Dropbox (PP)/ddml/Examples"
}

which ddml
which pystacked

********************************************************************************
*** interactive model     													 ***
********************************************************************************

webuse cattaneo2, clear
 
global Y bweight
global D mbsmoke
global X mage

gen fid=_n>2321

ddml init interactive, foldvar(fid)
ddml yeq, gen(regy): reg $Y $X
ddml deq, gen(regd): logit $D $X
ddml crossfit

ddml estimate m0   

local b = _b[mbsmoke]
local se = _se[mbsmoke] 
 
assert reldif(`b',-275.981 )<10e-4
assert reldif(`se',22.63938  ) <10e-4


ddml estimate m0  , atet

local b = _b[mbsmoke]
local se = _se[mbsmoke] 
 
 
di reldif(`b',-254.9799 ) 
di reldif(`se',21.62284  )  

assert reldif(`b',-254.9799 )<10e-4
assert reldif(`se',21.62284  ) <10e-4




********************************************************************************
*** LATE				 													 ***
********************************************************************************
 
use "http://fmwww.bc.edu/repec/bocode/j/jtpa.dta",clear   

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
ddml estimate 

local b = _b[training]
local se = _se[training] 
 
assert reldif(`b',-1767.331 )<10e-4
assert reldif(`se',513.3541 ) <10e-4




********************************************************************************
*** partial model	      													 ***
********************************************************************************

webuse cattaneo2, clear
 
global Y bweight
global D mbsmoke
global X mage

gen fid=_n>2321

ddml init partial, foldvar(fid)
ddml yeq, gen(regy): reg $Y $X
ddml deq, gen(regd): logit $D $X
ddml crossfit
ddml estimate m0 
 