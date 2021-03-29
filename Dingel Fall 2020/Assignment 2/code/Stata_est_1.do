* Start log
capture log close
log using Stata_est_1.log, replace

* Setup
clear all
set maxvar 32000

* Load data
use col_regfile09
* Only keep what we will use
keep flow distw iso_o year iso_d contig comlang_off
keep if year >= 2000 & year <= 2006
keep if flow > 0


* Log the data
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

* Table 1

* reg
* Start timer
timer on 1
* Quietly store regression
* Here we use FE dummies
eststo: quietly reg logflow logdist contig comlang_off i.year_iso_o i.year_iso_d
* Stop timer
timer off 1
timer list 1
* Add to list of estimate outputs for "time"
estadd local time = r(t1), replace

* xtreg
timer on 2
// xtset origin
xtset year_iso_o
// eststo: quietly xtreg logflow logdist contig comlang_off i.destination if year >= 2000 & year <= 2006, fe
eststo: quietly xtreg logflow logdist contig comlang_off i.year_iso_d, fe
timer off 2
timer list 2
estadd local time = r(t2), replace

* areg
timer on 3
// eststo: quietly areg logflow logdist contig comlang_off i.origin i.destination if year >= 2000 & year <= 2006, a(year)
eststo: quietly areg logflow logdist contig comlang_off i.year_iso_o, a(year_iso_d)
timer off 3
timer list 3
estadd local time = r(t3), replace

* reghdfe
timer on 4
eststo: quietly reghdfe logflow logdist contig comlang_off, a(year_iso_o year_iso_d)
timer off 4
timer list 4
estadd local time = r(t4), replace


// esttab
esttab using "Table_1.tex", tex label se title(Table 1) ///
	mgroups("reg" "xtreg" "areg" "reghdfe", pattern(1 1 1 1) prefix(\multicolumn{@span}{c}{) suffix(})) ///
	keep(logdist contig comlang_off) ///
	s(N r2 r2_a time, label("N" "R^2" "Adj. R2" "Time (sec)"))  replace














