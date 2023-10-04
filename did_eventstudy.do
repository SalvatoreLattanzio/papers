********************************************************************************
* Difference-in-differences and event studies
* Application using data from the paper 
* "Behind the Child Penalty: What Contributes to the Labour Market Costs of Motherhood"
* Alessandra Casarico and Salvatore Lattanzio (CESifo working paper, 2021)

* Salvatore Lattanzio
* Public Economics, Bocconi University, A.Y. 2022/23
********************************************************************************

clear
graph drop _all
set more off
set scheme s2color

if `"`c(os)'"' == "MacOSX"   global   stem   `"/users/`c(username)'/Dropbox"'
if `"`c(os)'"' == "Windows"  global   stem   `"C:/Users/`c(username)'/Dropbox"'

glo data  "${stem}/PubEconBocconi/lezioni/data"
glo log   "${stem}/PubEconBocconi/lezioni/log"
glo figs  "${stem}/PubEconBocconi/lezioni/figures"
glo table "${stem}/PubEconBocconi/lezioni/tables"

cap log close
log using ${log}/did_eventstudy, text replace

* Load dataset - social security matched employer employee data
use ${data}/losai_women, clear

* summarize data
describe
summarize

* treatment status
tab treated

* generate logs
gen learn = log(real_earnings)
gen lwage = log(wwage)
gen lweeks = log(weeks)
gen lweeksraw = log(weeksraw)

* descriptive graphs on the evolution of child penalty over time
preserve // we "save" the data
	collapse (mean) learn lweeks lwage part, by(k treated)
	
	foreach var in learn lweeks lwage part { // for loop over outcomes 
		
		if "`var'" == "learn" 	local lbe "Log annual earnings"
		if "`var'" == "lweeks" 	local lbe "Log weeks worked"
		if "`var'" == "lwage" 	local lbe "Log weekly wage"
		if "`var'" == "part" 	local lbe "Part-time = 1"
		
		two (connected `var' k if treated == 0, msymbol(O)) ///
			(connected `var' k if treated == 1, msymbol(S)), ///
			graphregion(color(white)) legend(order(1 "Women without children" 2 "Women with children")) ///
			xtitle("Years from first childbirth") ytitle("`lbe'") ///
			xline(-0.5) name(`var', replace)
	}

restore
grc1leg learn lweeks lwage part, graphregion(color(white))
graph export ${figs}/means.png, replace
graph drop _all


********************************************************************************
* Static DiD
********************************************************************************
gen post = (k >= 0)

reg learn i.post##i.treated i.year i.age, cluster(worker_id)
reg learn i.post##i.treated year age, cluster(worker_id)

* (almost) same specification as Kleven et al 2019
foreach var in learn lweeks lwage part { 
	reg `var' i.post##i.treated i.year i.age, cluster(worker_id)
	estimates store `var'
}
estout learn lweeks lwage part, keep(1.post#1.treated) cells(b se(par))

* invidual fixed effects
foreach var in learn lweeks lwage part { 
	areg `var' i.post##i.treated i.year, absorb(worker_id) cluster(worker_id)
	estimates store `var'
}
estout learn lweeks lwage part, keep(1.post#1.treated) cells(b se(par))
estimates clear


********************************************************************************
* Event study
********************************************************************************

* Event time dummies pre
forvalues i = 5 (-1) 2 {
	gen m_`i' = (k == -`i')
}

* auxiliary variable, just for graph
gen m_1 = 1

* Event time dummies post
forvalues i = 0/15 {
	gen p_`i' = (k == `i')
}

* Year and age dummies
xi i.year i.age

* Example with log earnings
* Separate regressions for women w/ and w/o children
reg learn m_* p_* _Iy* _Ia* if treated==1, cluster(worker_id)
estimates store treat
global treat15 = _b[p_15]

reg learn m_* p_* _Iy* _Ia* if treated==0, cluster(worker_id)
estimates store contr
global contr15 = _b[p_15]

* store long-run child penalty in a macro
global childpenalty : display %3.2f ${treat15} - ${contr15}

/* ssc install coefplot */ // uncomment if not already installed
coefplot (treat, ms(O)) (contr, ms(S)), vertical omitted baselevels ///
	xlabel(1 "-5" 6 "0" 11 "5" 16 "10" 21 "15")  ///
	recast(connected) ciopts(recast(rarea) color(%40) lwidth(none)) ///
	keep(m_5 m_4 m_3 m_2 m_1 p_0 p_1 p_2 p_3 p_4 p_5 p_6 p_7 p_8 p_9 p_10 p_11 p_12 p_13 p_14 p_15) ///
	xline(5.5, lc(black) lp(dash)) yline(0, lc(black)) ///
	legend(subtitle("Long-run penalty = $childpenalty", tstyle(body) pos(11)) ///
	order(2 "Women with children" 4 "Women without children") cols(1) pos(5) ring(0)) ///
	offset(0) xti("Years from maternity") yti("Coefficients relative to t-1") ///
	graphregion(color(white))
graph export ${figs}/childpenalty.png, replace


* Compare estimates with and without age fixed effects
reg learn m_* p_* _Iy* if treated==1, cluster(worker_id)
estimates store treat_noage
global treat15 = _b[p_15]

reg learn m_* p_* _Iy* if treated==0, cluster(worker_id)
estimates store contr_noage
global contr15 = _b[p_15]

global childpenalty_noage : display %3.2f ${treat15} - ${contr15}

coefplot (treat, ms(O)) (contr, ms(S)) (treat_noage, ms(O) lp(dash)) (contr_noage, ms(S) lp(dash)), ///
	vertical omitted baselevels ///
	xlabel(1 "-5" 6 "0" 11 "5" 16 "10" 21 "15")  ///
	recast(connected) ciopts(recast(rarea) color(%40) lwidth(none)) ///
	keep(m_5 m_4 m_3 m_2 m_1 p_0 p_1 p_2 p_3 p_4 p_5 p_6 p_7 p_8 p_9 p_10 p_11 p_12 p_13 p_14 p_15) ///
	xline(5.5, lc(black) lp(dash)) yline(0, lc(black)) ///
	legend(subtitle("Long-run penalty, age FE = $childpenalty" "Long-run penalty, no age FE = $childpenalty_noage", tstyle(body) pos(11) justification(left)) ///
	order(2 "Women with children, age FE" 4 "Women without children, age FE" ///
	6 "Women with children, no age FE" 8 "Women without children, no age FE") cols(1) pos(5) ring(0)) ///
	offset(0) xti("Years from maternity") yti("Coefficients relative to t-1") ///
	graphregion(color(white)) ylabel(-1.2(.2).4)

	
**************
* Decomposition graph
**************
* estimate coefficients for earnings, wage rates and weeks worked and store them in vectors
foreach var in learn lwage lweeksraw {
	reg `var' m_* p_* _Iy* _Ia* if treated==1, cluster(worker_id)
	mat `var'1 = (_b[m_5],_b[m_4],_b[m_3],_b[m_2],_b[m_1],_b[p_0],_b[p_1],_b[p_2],_b[p_3],_b[p_4],_b[p_5],_b[p_6],_b[p_7],_b[p_8],_b[p_9],_b[p_10],_b[p_11],_b[p_12],_b[p_13],_b[p_14],_b[p_15])'

	reg `var' m_* p_* _Iy* _Ia* if treated==0, cluster(worker_id)
	mat `var'0 = (_b[m_5],_b[m_4],_b[m_3],_b[m_2],_b[m_1],_b[p_0],_b[p_1],_b[p_2],_b[p_3],_b[p_4],_b[p_5],_b[p_6],_b[p_7],_b[p_8],_b[p_9],_b[p_10],_b[p_11],_b[p_12],_b[p_13],_b[p_14],_b[p_15])'
}

* graph
preserve

	* clear 
	drop _all
	
	* convert matrix into data
	svmat learn0
	svmat learn1
	svmat lwage0 
	svmat lwage1 
	svmat lweeksraw0
	svmat lweeksraw1
	
	* compute child penalty
	foreach var in learn lwage lweeksraw {
		gen p_`var' = `var'1 - `var'0 
	}
	
	* generate values for x-axis
	gen t=_n-6
	
	* auxiliary variable to draw "area" graphs
	gen x=0
	
	two (rarea p_learn x t, lwidth(none) fcolor(navy)) ///
		(rarea p_lweeks x t, lwidth(none) fcolor(dkorange)) ///
		(rarea p_lwage x t, lwidth(none) fcolor(gs10)) ///
		(line p_learn t, lwidth(1) lc(black)), ///
		legend(order(1 "Part-time related penalty" 2 "Weeks related penalty" ///
		3 "Wage related penalty" 4 "Total child penalty in earnings") symxsize(*0.5)) ///
		xti("Years from maternity")	yti("Child penalty, coefficients relative to t-1") ///
		xline(-0.5, lc(black) lp(dash)) yline(0, lc(black)) xlabel(-5(5)15) ///
		graphregion(color(white))
		
	graph export ${figs}/childpenalty_deco.png, replace
restore



capture log close

