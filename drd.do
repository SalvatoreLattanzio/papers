********************************************************************************
* Difference-in-discontinuities 
* Application using data from the paper 
* "Do Fiscal Rules Matter?"
* Veronica Grembi, Tommaso Nannici, Ugo Troiano (2016). AEJ: Applied 8(3): 1-30

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
log using ${log}/drd, text replace

* Load dataset - municipalities balance sheet and socio-demographic info 1997-2004
use ${data}/fiscal_aej.dta, clear

********************************************************************************
* Data cleaning
********************************************************************************
* diff-in-disc sample selection 
* - they keep only municipalities between 3500 and 7000 residents
* - and only years 1999-2004
drop if popcens<3500|popcens>7000
drop if anno<1999|anno>2004
sort id_comune anno 

* unbalanced panel
tab year

* create balanced sample with non-missing information on outcomes:
egen n1  = count(deficit_pc), by(id_codente)
egen n2  = count(saldofinanziario_pc), by(id_codente)
egen n3  = count(spesecor_pc), by(id_codente)
egen n4  = count(spesecocap_pc), by(id_codente)
egen n5  = count(expend_interest_pc), by(id_codente)
egen n6  = count(imposte), by(id_codente)
egen n7  = count(tasse), by(id_codente)
egen n8  = count(state_transfers), by(id_codente)
egen n9  = count(entrate_altre_pc), by(id_codente)
egen n10 = count(aliquota_ordinaria), by(id_codente)
egen n11 = count(other_transfers), by(id_codente)
egen n12 = count(trafgrants_pc), by(id_codente)
gen sample = (n1==6 & n2==6 & n3==6 & n4==6&n5==6&n6==6&n7==6&n8==6&n9==6&n10==6&n11==6&n12==6)
keep if sample==1
drop n1 n2 n3 n4 n5 n6 n7 n8 n9 n10 n11 n12 sample

* now we have a balanced panel
tab year

* variable labels
label variable deficit_pc "Deficit"       
label variable saldofinanziario_pc "Fiscal Gap"


********************************************************************************
* Generating variables for regressions
********************************************************************************
* normalized running variable
gen pop5000 = popcens - 5000

* generate TREATMENT STATUS:
gen treatment = pop5000<0

* post dummy for years 2001-2004
gen post = anno>=2001

* interaction treatment X running variable
gen treatment_pop5000 = treatment * pop5000

* interaction treatment X post
gen treatment_post = treatment * post

* interaction post X running variable
gen post_pop5000 = post * pop5000

* interaction treatment X post X running variable 
gen treatment_post_pop5000 = treatment_post * pop5000



********************************************************************************
* Tables
********************************************************************************

* Cols 1-2 - Fiscal Discipline (only Calonico et al. 2014 optimal bandwidth)
foreach var in deficit_pc saldofinanziario_pc {
	
	* compute optimal bandwidth
	rdrobust `var' pop5000 if post == 0
	local band1=e(h_l)
	rdrobust `var' pop5000 if post == 1
	local band2=e(h_l)
	local band=(`band1'+`band2')/2
	display "`band'"
	
	* parametric regression
	reg `var' treatment_post treatment post pop5000 ///
		post_pop5000 treatment_pop5000 ///
		treatment_post_pop5000 ///
		if abs(pop5000)<`band', cluster(id_codente)
	estimates store `var'
}
estout deficit_pc saldofinanziario_pc, cells(b(star fmt(2)) se(fmt(2) par)) keep(treatment_post) 
/* Note: since the publication of the paper, Calonico et al have updated their
	code for the computation of the optimal bw, which is different here than 
	that reported in the paper. If you use same bw as those in the paper, results
	are just the same */
	

********************************************************************************
* Sensitivity to bandwidth
********************************************************************************
foreach var in deficit_pc saldofinanziario_pc { //loop over dep. vars.
	
	* save variable label in a local
	local lbe : variable label `var'

	* empty matrix to store coefficients
	mat coef = J(13,3,.)
	local i = 1
	
	forvalues band=300(100)1500 { //loop over bandwidths
		
		* regressions
		reg `var' treatment_post treatment post pop5000 ///
			post_pop5000 treatment_pop5000 ///
			treatment_post_pop5000 ///
			if abs(pop5000)<`band', cluster(id_codente)
		
		* store beta and confidence intervals
		mat coef[`i',1] = _b[treatment_post]
		mat coef[`i',2] = _b[treatment_post] - 1.96*_se[treatment_post]
		mat coef[`i',3] = _b[treatment_post] + 1.96*_se[treatment_post]
	
	local ++i
		
	}
	coefplot matrix(coef[,1]), ci((coef[,2] coef[,3])) vert ///
		ciopts(recast(rline) lpattern(dash)) recast(connected) ///
		xlabel(1 "300" 4 "600" 7 "900" 10 "1200" 13 "1500") ///
		yline(0) graphregion(color(white)) xtitle("Bandwidth") ytitle("`lbe'") name(`var', replace)
}
graph combine deficit_pc saldofinanziario_pc, graphregion(color(white))
graph export ${figs}/bwsens.png, replace



********************************************************************************
* Parallel trends
********************************************************************************
foreach var in deficit_pc saldofinanziario_pc { //loop over dep. vars.
	
	* save variable label in a local
	local lbe : variable label `var'

	* empty matrix to store coefficients
	mat coef = J(6,3,.)
	local i = 1
	
	forvalues y=1999/2004 { //loop over years
		
		* bandwidth
		rdrobust `var' pop5000 if anno == `y'
		local band=e(h_l)

		* regressions
		reg `var' treatment pop5000 treatment_pop5000 ///
			 if anno == `y' & abs(pop5000)<`band', robust
		
		* store beta and confidence intervals
		mat coef[`i',1] = _b[treatment]
		mat coef[`i',2] = _b[treatment] - 1.96*_se[treatment]
		mat coef[`i',3] = _b[treatment] + 1.96*_se[treatment]
	
	local ++i
		
	}
	coefplot matrix(coef[,1]), ci((coef[,2] coef[,3])) vert ///
		ciopts(recast(rline) lpattern(dash)) recast(connected) ///
		xlabel(1 "1999" 2 "2000" 3 "2001" 4 "2002" 5 "2003" 6 "2004") ///
		yline(0) graphregion(color(white)) xtitle("Year") ytitle("`lbe'") name(`var', replace)
}
graph combine deficit_pc saldofinanziario_pc, graphregion(color(white))
graph export ${figs}/rdyear.png, replace
/* Note: the paper uses rdrobust (non-parametric estimates) to produce the 
	yearly RD estimates. Here I use the same approach used in the rest of the paper,
	i.e. rdrobust to compute optimal bandwidth and then parametric estimates.
	One weird thing is that rdrobust, by default, consider as treated units those 
	above the threshold (while here those below the threshold are treated). 
	Hence, figure 2, panel A and B, reports coefficients of opposite sign
	w.r.t. those reported in Table 4 */


********************************************************************************
* Density of the forcing variable
********************************************************************************
* density of population before 2001
cap drop Xj Yj r0 *fhat*
DCdensity popcens if anno<2001, breakpoint(5000) generate(Xj Yj r0 fhat se_fhat)
graph display, name(g1, replace)

* density of population after 2001
cap drop Xj Yj r0 *fhat*
DCdensity popcens if anno>=2001, breakpoint(5000) generate(Xj Yj r0 fhat se_fhat)
graph display, name(g2, replace)

graph combine g1 g2
graph export ${figs}/mccrary.png, replace


********************************************************************************
* Graphs
********************************************************************************
* generate first differences wrt 1999 & 2000
rename (saldofinanziario_pc deficit_pc) (saldo deficit)
keep id_codente anno post saldo deficit popcens 

* we do this in a for loop for fiscal gap (saldo) and deficit
foreach var in saldo deficit {
	forvalues y = 1999/2000 {
		gen temp = `var' if anno == `y'
		bysort id_codente (anno): egen temp2 = max(temp)
		gen diff_`var'`y' = `var' - temp2
		drop temp*
	}
}

* keep only post-policy period
keep if post == 1

* reshape data to long format to have differences with respect to both pre-periods
drop saldo deficit post
reshape long diff_saldo diff_deficit, i(anno id_codente) j(baseyear)
sort id_codente

* normalized running variable & polynomial
gen pop5000 = popcens - 5000
gen pop5000_2 = pop5000 * pop5000

* treatment status
gen treatment = pop5000<0

* interactions
gen treatment_int1 = treatment * pop5000
gen treatment_int2 = treatment * pop5000 * pop5000

* population bins
xtile bin=pop5000, nq(47) // using same bins as paper (it's arbitrary)

* average population in each bin (just for graphs)
egen newx=mean(pop5000), by(bin)

* labels
label variable diff_saldo "Fiscal gap"
label variable diff_deficit "Deficit"

foreach var in diff_saldo diff_deficit {
	
	* save variable label in a local
	local lbe : variable label `var'
	
	* mean outcomes in each population bin
	egen newy=mean(`var'), by(bin)

	* regression 
	reg `var' treatment pop5000 pop5000_2 treatment_int1 treatment_int2
	predict yhat
	predict SE, stdp 
	gen low = yhat - 1.96*(SE)
	gen high = yhat + 1.96*(SE)	

	twoway ///
		(scatter newy newx, msymbol(circle_hollow) mcolor(gray)) ///
		(line yhat low high pop5000 if treatment == 1, pstyle(p1 p3 p3) sort) ///
		(line yhat low high pop5000 if treatment == 0, pstyle(p1 p3 p3) sort) if inrange(pop5000,-1500,1500), ///
		ytitle("`lbe'") xtitle("Normalized population") graphregion(color(white)) legend(off) ///
		xline(0) xlabel(-1500 0 1500) ylabel(#3) name(`var', replace)
	drop yhat low high SE newy
}
graph combine diff_saldo diff_deficit, graphregion(color(white))
graph export ${figs}/drd_graphs.png, replace

/* Note: 
	the authors use a 3rd order polynomial on both sides of the threshold,
	but Gelman and Imbens (2019) advise against going above 2nd order, 
	which is what is used here. Therefore, the picture does not look exactly 
	the same (Moreover, they do an asymmetric winsorization of the fiscal gap 
	variable, which further explains differences with the paper */

log close