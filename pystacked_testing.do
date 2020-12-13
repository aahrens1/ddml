
clear all
 
if ("`c(username)'"=="kahrens") {
	adopath + "/Users/kahrens/MyProjects/ddml"
	cd "/Users/kahrens/Dropbox (PP)/ddml"
}

insheet using https://web.stanford.edu/~hastie/ElemStatLearn/datasets/prostate.data,  tab

//replace lcavol = . if _n ==30
 
pystacked lpsa lcavol lweight age lbph svi lcp gleason pgg45, type(regress) lasso rf grad

predict xb if _n<95, xb  
list xb*
predict xb if _n<95, transform
list xb* 
