********************************************************************************
* Replication of the paper "Gender Interactions within Hierarchies"
* Stefano Gagliarducci and Daniele Paserman (Review of Economic Studies, 2012)

* Salvatore Lattanzio
* Public Economics, Bocconi University, A.Y. 2022/23
********************************************************************************

clear
graph drop _all
set more off
set scheme s2color

if `"`c(os)'"' == "MacOSX"   global   stem   `"/users/`c(username)'/Dropbox"'
if `"`c(os)'"' == "Windows"  global   stem   `"C:/Users/`c(username)'/Dropbox"'

global data  "${stem}/PubEconBocconi/lezioni/data"
glo log   "${stem}/PubEconBocconi/lezioni/log"
glo figs  "${stem}/PubEconBocconi/lezioni/figures"
glo table "${stem}/PubEconBocconi/lezioni/tables"

cap log close
log using ${log}/rdd, text replace

* Load dataset - municipality panel on mayors and candidates to mayor position
use "${data}/SIN19932003.dta", clear
describe

* Sharp design
tabstat votes_gender, by(gender) stats(mean min max)

twoway scatter gender votes_gender, msize(vsmall) graphregion(c(white)) xline(0)
graph export ${figs}/sharp.png, replace

********************************************************************************
* Figure 1 - Early termination probability
********************************************************************************
* 2% Bins of margin of victory
gen votes_gender_grp = 2*int((votes_gender)/2)
egen p_dimiss1 = mean(early_termination1), by(votes_gender_grp)
*bysort votes_gender_grp: egen p_dimiss1 = mean(early_termination1)
keep if p_dimiss1<1
egen  mean_votes_gender = mean(votes_gender), by(votes_gender_grp)
two (lowess p_dimiss1 votes_gender if votes_gender<0 & p_dimiss1<1,mean lcolor(blue)) ///
	(lowess p_dimiss1 votes_gender if votes_gender>=0,mean lcolor(red) lpattern(dash))  ///
	(scatter p_dimiss1 mean_votes_gender if p_dimiss1<1, ///
	msize(small) msymb(circle) mc(black) mlab(votes_gender_grp) mlabcolor(black) mlabsize(vsmall)), ///
	xtitle("Relative % votes female candidate") ytitle("Early termination (smooth average)") ///
	xscale(range(-50,50)) xlabel(-50(10)50,labsize(vsmall)) yscale(range(0.02,0.25)) ///
	ylabel(0.00(0.05)0.25,labsize(vsmall)) xline(0,lwidth(small)) ///
	legend(label(1 "winning man") label(2 "winning woman") label(3 "observed values") colgap(*3) width(100)) ///
	xscale(titlegap(3)) yscale(titlegap(3)) graphregion(color(white))
graph export ${figs}/fig1.png, replace

* alternative commands
* findit rdrobust //and then click and install
rdplot early_termination1 votes_gender , p(2) nbins(20 20) 
rdplot early_termination1 votes_gender  if abs(votes_gender)<20, p(1) nbins(20 20) 
cmogram early_termination1 votes_gender, lfitci scatter cut(0) lineat(0) graphopts(xtitle("Relative % votes female candidate") ytitle("Early termination (smooth average)"))
cmogram early_termination1 votes_gender if abs(votes_gender)<20, lfitci scatter cut(0) lineat(0) graphopts(xtitle("Relative % votes female candidate") ytitle("Early termination (smooth average)"))
cmogram early_termination1 votes_gender if abs(votes_gender)<20, qfitci scatter cut(0) lineat(0) graphopts(xtitle("Relative % votes female candidate") ytitle("Early termination (smooth average)"))

* "by hand"
reg early_termination1 gender votes_gender genderXvotes_gender, robust
predict fit
predict fitsd, stdp
gen upfit=fit+1.96*fitsd
gen downfit=fit-1.96*fitsd
twoway ///
	(line fit upfit downfit votes_gender if votes_gender<0, sort pstyle(p1 p3 p3))  ///
	(line fit upfit downfit votes_gender if votes_gender>0, sort pstyle(p1 p3 p3))  ///
	(scatter p_dimiss1 votes_gender_grp, mcolor(gray) msymbol(circle_hollow)) if p_dimiss1<1, ///
	ytitle("Early termination (smooth average)") xtitle("Relative % votes female candidate") ///
	xlabel(-50(10)50) xline(0) graphregion(c(white)) ///
	legend(off)
	

********************************************************************************
* Figure 2: Frequency of margin of victory, mixed-gender races
********************************************************************************
use "${data}/SIN19932003.dta", clear
twoway (histogram votes_gender, color(sienna) lcolor(olive) freq width(2) lwidth(vvthin)) ///
	(histogram votes_gender, color(sand) lcolor(olive) freq width(1) lwidth(vvthin)) ///
	(histogram votes_gender, color(brown) lcolor(olive) freq width(0.5) lwidth(vvthin)), ///
	xline(0) xtitle("Relative % votes female candidate") ytitle("Absolute frequency") ///
	xlabel(-40(10)40, labsize(small)) ylabel(0(50)200, labsize(small)) xscale(titlegap(3)) yscale(titlegap(3)) ///
	legend(label(1 "bandwidth=2%") label(2 "bandwidth=1%") label(3 "bandwidth=0.5%") colgap(*3) width(100)) graphregion(color(white))
graph export ${figs}/density.png, replace

* more formal McCrary test (McCrary, 2008)
/* To install DCdensity see README here: https://eml.berkeley.edu/~jmccrary/DCdensity/ */
capture drop Xj Yj r0 *fhat*
qui DCdensity votes_gender, breakpoint(0) generate(Xj Yj r0 fhat se_fhat) nograph
local bmc : display %4.3f `r(theta)'
local smc : display %4.3f `r(se)'
gen upfhat=fhat+1.96*se_fhat
gen lofhat=fhat-1.96*se_fhat
twoway (scatter Yj Xj if abs(Xj)<40, ms(circle_hollow) mc(gray)) ///
	(line fhat upfhat lofhat r0 if r0<0 & abs(r0)<40, pstyle(p1 p3 p3)) ///
	(line fhat upfhat lofhat r0 if r0>0 & abs(r0)<40, pstyle(p1 p3 p3)), ///
	xline(0) xti(Margin of victory) yti(Density) ///
	legend(order(- "McCrary test: `bmc' (`smc')") ///
	cols(1) pos(2) ring(0) region(lc(black))) name(mccrary, replace)
graph export ${figs}/mccrary.png, replace

* another test [Cattaneo, Jansson, and Ma, 2017]
rddensity votes_gender, plot


********************************************************************************
* Figure 3: Balance tests, mixed-gender races: city characteristics
********************************************************************************
use "${data}/SIN19932003.dta", clear
preserve
	keep if votes_gender>=-30 & votes_gender<=30 //sample restriction (arbitrary)
	gen votes_gender_grp = 2*int((votes_gender)/2) //mv bins
	egen  mean_votes_gender = mean(votes_gender), by(votes_gender_grp)
	* compute means of covariates in each bin
	foreach var of varlist pop_res elderly_index income_pc active_pop n_firms_pc NW NE CE SOU {
	  egen `var'_mean = mean(`var'), by(votes_gender_grp)
	}
	
	replace pop_res_mean=pop_res_mean/1000 // divide by 1000 just for better visualization
	replace income_pc_mean=income_pc_mean/1000
	* labels
	la var pop_res_mean "Population"
	la var NW_mean "North-west"
	la var NE_mean "North-east"
	la var CE_mean "Center"
	la var SOU_mean "South"
	la var active_pop_mean "Active population"
	la var n_firms_pc_mean "Firms per capita"
	la var elderly_index_mean "Old-age index"
	la var income_pc_mean "Income pc"
	
	* graphs
	foreach var of varlist pop_res elderly_index income_pc active_pop n_firms_pc NW NE CE SOU {
		local lab : variable label `var'_mean
		
		two ///
			(lowess `var'_mean votes_gender if votes_gender>=-30 & votes_gender<0, mean lcolor(gs4)) ///
			(lowess `var'_mean votes_gender if votes_gender>=0 & votes_gender<=30, mean lcolor(gs8) lpattern(dash))  ///
			(scatter `var'_mean mean_votes_gender if votes_gender>=-30 & votes_gender<=30, msize(small) msymb(circle) mcolor(black)), ///
			xtitle("Relative % votes female") name(RDDgraph_`var', replace) graphregion(color(white)) ///
			ytitle(" `lab' ") xlabel(-30(10)30) ylabel(#3) xline(0, lcolor(black)) ///
			legend(label(1 "winning man") label(2 "winning woman") label(3 "observed values") colgap(*3) width(100))
	}
restore

* combine graphs
grc1leg RDDgraph_pop_res RDDgraph_NW RDDgraph_NE RDDgraph_CE RDDgraph_SOU RDDgraph_active_pop RDDgraph_n_firms_pc RDDgraph_elderly_index RDDgraph_income_pc, graphregion(color(white))
graph export ${figs}/balancetest.png, replace


*****************************************************************************
* Table 3: The Effect of Mayor Gender on the Probability of Early Termination
*****************************************************************************
use ${data}/SIN19932003.dta, clear


#delimit ;
global controls "NE NW SOU IS 
	a1993 a1994 a1995 a1996 a1997 a1998 a1999 a2000 a2001 a2002
	age term_limit n_terms exp born_prov born_prov_miss lowsecondary upsecondary college educ_miss  
	empl_not empl_profess empl_entrepeneur empl_white empl_miss
	logpop_res lnincome_pc elderly_index active_pop n_firms_pc demographic_miss  
	electors_female_per electors_female_per_miss seats_counc seats_maycoal_per seats_maycoal_per_miss 
	number_parties_maycoal pct_councilors_mparty pct_councilors_mparty_miss 
	centerright2 centerleft2 separatist2 center2" ;
#delimit cr

* Add percentage of councilors in the mayor's *party*
gen pct_councilors_mparty = 1-polar_maycoal
replace pct_councilors_mparty = 0 if polar_maycoal==.
gen pct_councilors_mparty_miss = (polar_maycoal==.)

*** Top panel - Full sample

* Determine sample
qui reg early_termination1 gender ${controls}, cluster(id_city_istat)
keep if e(sample)

* column 1 - unadjusted estimates
reg early_termination1 gender, cluster(id_city_istat)		
xtreg early_termination1 gender, cluster(id_city_istat) fe i(id_city_istat)
		
* column 2 - adjusted estimates
reg early_termination1 gender ${controls}, cluster(id_city_istat)
xtreg early_termination1 gender ${controls}, i(id_city_istat) fe cluster(id_city_istat)


*** Bottom panel - mixed gender elections only
use ${data}/SIN19932003.dta, clear
glo controls "age lowsecondary upsecondary college educ_miss empl_not empl_profess empl_entrepeneur empl_white empl_miss born_prov born_prov_miss n_terms exp electors_female_per electors_female_per_miss centerleft2 separatist2 center2 centerright2 a1993 a1994 a1995 a1996 a1997 a1998 a1999 a2000 a2001 a2002 NW NE CE SOU IS logpop_res elderly_index active_pop income_pc n_firms_pc demographic_miss"

* Determine sample - drop non-mixed gender elections
qui reg early_termination1 gender /*
*/ age lowsecondary upsecondary college educ_miss empl_not empl_profess empl_entrepeneur empl_white empl_miss born_prov born_prov_miss n_terms exp electors_female_per electors_female_per_miss centerleft2 separatist2 center2 centerright2 a1993 a1994 a1995 a1996 a1997 a1998 a1999 a2000 a2001 a2002 NW NE CE SOU IS logpop_res elderly_index active_pop income_pc n_firms_pc demographic_miss /*
*/ if votes_gender!=., cluster(id_city_istat)
keep if e(sample)


* linear regression
reg early_termination1 gender, cluster(id_city_istat)
reg early_termination1 gender ${controls}, cluster(id_city_istat)
			
* local linear regression
reg early_termination1 gender votes_gender genderXvotes_gender, cluster(id_city_istat)
reg early_termination1 gender votes_gender genderXvotes_gender ${controls}, cluster(id_city_istat)

* local linear regression, two candidates
reg early_termination1 gender votes_gender genderXvotes_gender if number_rivals==1, cluster(id_city_istat)
reg early_termination1 gender votes_gender genderXvotes_gender ${controls} if number_rivals==1, cluster(id_city_istat)

* local linear regression, optimal bw
reg early_termination1 gender votes_gender genderXvotes_gender if abs(votes_gender)<=25, cluster(id_city_istat)
reg early_termination1 gender votes_gender genderXvotes_gender ${controls} if abs(votes_gender)<=25, cluster(id_city_istat)

* local linear regression, two candidates, optimal bw
reg early_termination1 gender votes_gender genderXvotes_gender if abs(votes_gender)<=29 & number_rivals==1, cluster(id_city_istat)
reg early_termination1 gender votes_gender genderXvotes_gender ${controls} if abs(votes_gender)<=29 & number_rivals==1, cluster(id_city_istat)

* local linear regression, half optimal bw
reg early_termination1 gender votes_gender genderXvotes_gender if abs(votes_gender)<=13, cluster(id_city_istat)
reg early_termination1 gender votes_gender genderXvotes_gender ${controls} if abs(votes_gender)<=13, cluster(id_city_istat)

* second-order polynomial on both sides of the cut-off
gen votes_gender2 = votes_gender^2
gen genderXvotes_gender2 = gender*votes_gender2
reg early_termination1 gender votes_gender votes_gender2 genderXvotes_gender genderXvotes_gender2 , cluster(id_city_istat)
*reg early_termination1 i.gender##c.votes_gender##c.votes_gender
reg early_termination1 gender votes_gender votes_gender2 genderXvotes_gender genderXvotes_gender2 ${controls}, cluster(id_city_istat)

* second-order polynomial on both sides of the cut-off, two candidates
reg early_termination1 gender votes_gender votes_gender2 genderXvotes_gender genderXvotes_gender2 if number_rivals==1, cluster(id_city_istat)
reg early_termination1 gender votes_gender votes_gender2 genderXvotes_gender genderXvotes_gender2 ${controls} if number_rivals==1, cluster(id_city_istat)


********************************************************************************
* Extras
********************************************************************************
	
* extra 1a - non-parametric estimation and optimal bandwidth from Calonico, Cattaneo and Titiunik
rdrobust early_termination1 votes_gender, all

* extra 1b - parametric estimation and optimal bandwidth from CCT
quietly rdrobust early_termination1 votes_gender
reg early_termination1 gender votes_gender genderXvotes_gender if abs(votes_gender)<=e(h_l), cluster(id_city_istat)


* extra 2 - sensitivity to bandwidth
matrix coef=J(20,4,.)
local i=1
* estimate regressions for different bandwidths (between 2% and 40%)
forvalues h=2(2)40 {
	reg early_termination1 gender votes_gender genderXvotes_gender if abs(votes_gender)<`h', cluster(id_city_istat)
	mat coef[`i',1]=`h'
	mat coef[`i',2]=_b[gender]
	mat coef[`i',3]=_b[gender]+1.96*_se[gender]
	mat coef[`i',4]=_b[gender]-1.96*_se[gender]
	local ++i
}
mat colnames coef=bw b ul ll

* report estimates in a graph
preserve
	drop _all
	svmat coef, names(col)
	two ///
		(scatter b bw, mc(black)) (rcap ul ll bw, lc(black)), ///
		legend(off) xtitle("Bandwidth (%)") graphregion(colo(white)) ///
		xlabel(2(2)40) yline(0, lp(dash)) ytitle("Coefficient")
	graph export ${figs}/bwsens.png, replace
restore

* extra 3 - sensitivity to polynomial order (now less "popular" in the RD literature
* see Andrew Gelman & Guido Imbens (2019) Why High-Order Polynomials Should Not Be Used in Regression Discontinuity Designs, Journal of Business & Economic Statistics, 37:3, 447-456)
mat coef=J(6,3,.)
gen votes_gender3=votes_gender^3
gen votes_gender4=votes_gender^4
gen votes_gender5=votes_gender^5
gen votes_gender6=votes_gender^6
gen genderXvotes_gender3=gender*votes_gender3
gen genderXvotes_gender4=gender*votes_gender4
gen genderXvotes_gender5=gender*votes_gender5
gen genderXvotes_gender6=gender*votes_gender6

* poly order 1
reg early_termination1 gender votes_gender genderXvotes_gender , cluster(id_city_istat)
	mat coef[1,1]=_b[gender]
	mat coef[1,2]=_b[gender]+1.96*_se[gender]
	mat coef[1,3]=_b[gender]-1.96*_se[gender]
	
* poly order 2
reg early_termination1 gender votes_gender votes_gender2 genderXvotes_gender genderXvotes_gender2 , cluster(id_city_istat)
	mat coef[2,1]=_b[gender]
	mat coef[2,2]=_b[gender]+1.96*_se[gender]
	mat coef[2,3]=_b[gender]-1.96*_se[gender]
	
* poly order 3
reg early_termination1 gender votes_gender votes_gender2 votes_gender3 genderXvotes_gender genderXvotes_gender2 genderXvotes_gender3, cluster(id_city_istat)
	mat coef[3,1]=_b[gender]
	mat coef[3,2]=_b[gender]+1.96*_se[gender]
	mat coef[3,3]=_b[gender]-1.96*_se[gender]
	
* poly order 4
reg early_termination1 gender votes_gender votes_gender2 votes_gender3 votes_gender4 genderXvotes_gender genderXvotes_gender2 genderXvotes_gender3 genderXvotes_gender4, cluster(id_city_istat)
	mat coef[4,1]=_b[gender]
	mat coef[4,2]=_b[gender]+1.96*_se[gender]
	mat coef[4,3]=_b[gender]-1.96*_se[gender]
	
* poly order 5
reg early_termination1 gender votes_gender votes_gender2 votes_gender3 votes_gender4 votes_gender5 genderXvotes_gender genderXvotes_gender2 genderXvotes_gender3 genderXvotes_gender4 genderXvotes_gender5, cluster(id_city_istat)
	mat coef[5,1]=_b[gender]
	mat coef[5,2]=_b[gender]+1.96*_se[gender]
	mat coef[5,3]=_b[gender]-1.96*_se[gender]
	
* poly order 6
reg early_termination1 gender votes_gender votes_gender2 votes_gender3 votes_gender4 votes_gender5 votes_gender6 genderXvotes_gender genderXvotes_gender2 genderXvotes_gender3 genderXvotes_gender4 genderXvotes_gender5 genderXvotes_gender6, cluster(id_city_istat)
	mat coef[6,1]=_b[gender]
	mat coef[6,2]=_b[gender]+1.96*_se[gender]
	mat coef[6,3]=_b[gender]-1.96*_se[gender]

mat colnames coef=b ul ll

* report estimates in a graph
preserve
	drop _all
	svmat coef, names(col)
	gen X=_n
	two ///
		(scatter b X, mc(black)) (rcap ul ll X, lc(black)), ///
		legend(off) xtitle("Polynomial order") graphregion(color(white)) ///
		xlabel(1(1)6) yline(0, lp(dash)) ytitle("Coefficient")
	graph export ${figs}/polordersens.png, replace
restore

log close





