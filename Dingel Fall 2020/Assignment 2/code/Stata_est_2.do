* Start log
capture log close
log using Stata_est_2.log, replace

* Setup
clear all
set maxvar 32000

* Load data
use col_regfile09
* Only keep what we will use
keep flow distw iso_o year iso_d contig comlang_off
keep if year >= 2004 & year <= 2006


* Log the data and 1+data
gen logflow = log(flow)
gen logdist = log(distw)
gen logflow_1 = log(flow + 1)

* Make dummies
egen year_iso_o = group(year iso_o)
egen year_iso_d = group(year iso_d)

* label variables for table
label var contig "Contiguous countries (binary)"
label var comlang_off "Common language (binary)"
label var logdist "log(Distance)"
label var logflow "log(Flows)"
label var logflow_1 "log(1+Flows)"

* Clear estimate storage and timer
eststo clear
timer clear

* Table 2

* 1. A log-linear regression that omits observations in which flow equals zero.
timer on 1
eststo: quietly reghdfe logflow logdist contig comlang_off if (flow > 0), a(year_iso_o year_iso_d) resid
predict fitted
predict resid, resid
timer off 1
timer list 1
estadd local time = r(t1), replace
estadd local zeros "No"
estadd local log_plus "No"
estadd local regtype "Log-Linear"

* 2. A log-linear regression that omits observations in which flow equals zero. In addition, set the dependent variable to log of flow plus one.
timer on 2
eststo: quietly reghdfe logflow_1 logdist contig comlang_off if (flow > 0), a(year_iso_o year_iso_d)
timer off 2
timer list 2
estadd loca time = r(t2), replace
estadd local zeros "No"
estadd local log_plus "Yes"
estadd local regtype "Log-Linear"

* 3. A log-linear regression in which the dependent variable is log of flow plus one.
timer on 3
eststo: quietly reghdfe logflow_1 logdist contig comlang_off, a(year_iso_o year_iso_d)
timer off 3
timer list 3
estadd loca time = r(t3), replace
estadd local zeros "Yes"
estadd local log_plus "Yes"
estadd local regtype "Log-Linear"

* 4. An estimate of the same constant-elasticity specification that uses the ppml command to implement the PPML estimator of Silva and Tenreyro (REStat 2006). Use all observations, including zeros, in this and following three columns.
timer on 4
qui tab year_iso_o, gen(yo_)
qui tab year_iso_d, gen(yd_)
eststo: ppml flow logdist contig comlang_off yo_* yd_*
timer off 4
timer list 4
estadd loca time = r(t4), replace
estadd local zeros "Yes"
estadd local log_plus "No"
estadd local regtype "PPML"


* 5. An estimate of the same constant-elasticity specification that uses the poi2hdfe command to implement the PPML estimator.
timer on 5
eststo: quietly poi2hdfe flow logdist contig comlang_off, id1(year_iso_o) id2(year_iso_d)
timer off 5
timer list 5
estadd loca time = r(t5), replace
estadd local zeros "Yes"
estadd local log_plus "No"
estadd local regtype "PPML"

* 6. An estimate of the same constant-elasticity specification that uses the ppml_panel_sg command to implement the PPML estimator.
timer on 6
eststo: ppml_panel_sg flow logdist contig comlang_off, exporter(iso_o) importer(iso_d) year(year) nopair olsguess nocheck
timer off 6
timer list 6
estadd loca time = r(t6), replace
estadd local zeros "Yes"
estadd local log_plus "No"
estadd local regtype "PPML"

* 7. An estimate of the same constant-elasticity specification that uses the ppmlhdfe command to implement the PPML estimator.
timer on 7
eststo: quietly ppmlhdfe flow logdist contig comlang_off, a(year_iso_o year_iso_d)
timer off 7
timer list 7
estadd loca time = r(t7), replace
estadd local zeros "Yes"
estadd local log_plus "No"
estadd local regtype "PPML"

* 8. An estimate of the same constant-elasticity specification that uses the ppmlhdfe command to implement the PPML estimator. Omit observations in which flow equals zero.
timer on 8
eststo: quietly ppmlhdfe flow logdist contig comlang_off if (flow > 0), a(year_iso_o year_iso_d)
timer off 8
timer list 8
estadd loca time = r(t8), replace
estadd local zeros "No"
estadd local log_plus "No"
estadd local regtype "PPML"


esttab using "Table_2.tex", tex label se title(Table 2) ///
	mgroups("reghdfe" "reghdfe" "reghdfe" "ppml" "poi2hdfe" "ppml_panel_sg" "ppmlhdfe" "ppmlhdfe", ///
	pattern(1 1 1 1 1 1 1 1) prefix(\multicolumn{@span}{c}{) suffix(})) ///
	keep(logdist contig comlang_off) ///
	s(N r2 r2_a zeros log_plus regtype time, label("N" "R^2" "Adj. R2" "Include zero flows?" "log(flows+1)?" "Approach" "Time (sec)")) replace

	
* Breusch-Pagan test
quietly reg logflow logdist contig comlang_off i.year_iso_o i.year_iso_d if (flow > 0)
hettest

* Residual scatterplot
scatter resid fitted, graphregion(color(white)) msize(.3) mcolor(blue)
graph export residuals.png, replace
















