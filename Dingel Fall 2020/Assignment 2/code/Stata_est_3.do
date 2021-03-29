* Start log
capture log close
log using Stata_est_3.log, replace

* Setup
clear all
set maxvar 32000

* Load data
use col_regfile09
* Only keep what we will use
keep flow distw iso_o year iso_d contig comlang_off
keep if (flow > 0)


* Log the data and 1+data
gen logflow = log(flow)
gen logdist = log(distw)


* Make  dummies
egen year_iso_o = group(year iso_o)
egen year_iso_d = group(year iso_d)


* label variables for table
label var contig "Contiguous countries (binary)"
label var comlang_off "Common language (binary)"
label var logdist "log(Distance)"
label var logflow "log(Flows)"


* Clear estimate storage and timer
eststo clear
timer clear


* Table 3

timer on 1
eststo: quietly reghdfe logflow logdist contig comlang_off, a(year_iso_o year_iso_d)
timer off 1
timer list 1
estadd local time = r(t1), replace

esttab using "Table_3_Stata.tex", tex label se title(Table 3) ///
	mgroups("Stata" "R" "Julia", pattern(1 1 1) prefix(\multicolumn{@span}{c}{) suffix(})) ///
	keep(logdist contig comlang_off) ///
	s(N r2 r2_a time, label("N" "R^2" "Adj. R2" "Time (sec)"))  replace
















