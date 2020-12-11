insheet using https://web.stanford.edu/~hastie/ElemStatLearn/datasets/prostate.data, clear tab

timer clear
cap profiler clear

timer on 1
profiler on
forvalues i=1/100 {
	qui pylasso2 lpsa lcavol lweight age lbph svi lcp gleason pgg45, lambda(0.2)
}
profiler off
timer off 1
di "pylasso2"
profiler report
profiler clear

timer on 2
profiler on
di
forvalues i=1/100 {
	qui lasso2 lpsa lcavol lweight age lbph svi lcp gleason pgg45, lglmnet lambda(0.2) prestd
}
profiler off
timer off 2
di "lasso2"
profiler report
profiler clear
di
timer list