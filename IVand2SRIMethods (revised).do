/* Stata code for IV and 2SRI for multilevel models*/

/* Code by Eliot R. Weinstein 
eliotw@uchicago.edu and erweinstein.com */
* Examples based on my own work, where I'm studying "authors" and "physicians"
* And I have outcome variables that are highly skewed but not suitable for zero-inflated poisson, hence my use of negative binomial


* First stage scatter and regression
use predicted_LOS_by_author_menbreg.dta, clear
capture confirm variable endog
if _rc {
    display "ERROR: replace 'endog' with your endogenous regressor variable name"
    exit 1
}
* Replace endog and instr below if your variable names differ
local endog = "mean_LOS_MJ1"   // example endogenous regressor; change if needed
local instr = "instr_var"      // replace with your instrument name

twoway (scatter `endog' `instr', msize(medium) mcolor(navy)) (lfit `endog' `instr', lcolor(maroon) lwidth(medium)), ///
    title("First stage: `endog' vs `instr'") xtitle("Instrument (`instr')") ytitle("Endogenous regressor (`endog')")

regress `endog' `instr'
display "Check the t and F for `instr' in the regress output (partial F should be large)"


* First-stage residual diagnostics
use predicted_LOS_by_author_menbreg.dta, clear
local endog = "mean_LOS_MJ1"
local instr = "instr_var"
regress `endog' `instr'
predict double endog_hat, xb
predict double endog_resid, resid

histogram endog_resid, normal width(0.5) title("Histogram of first-stage residuals")
twoway (scatter endog_resid endog_hat, msize(small) mcolor(black)) (lfit endog_resid endog_hat, lcolor(gs12)), ///
    title("First-stage residuals vs predicted") xtitle("Predicted `endog'") ytitle("Residuals")


* 2SLS vs OLS comparison (linear outcome)
use predicted_LOS_by_author_menbreg.dta, clear
local outcome = "outcome_var"   // replace with your outcome
local endog   = "mean_LOS_MJ1"  // replace if different
local instr   = "instr_var"     // replace with your instrument
local controls = "cov1 cov2"    // replace with your covariates or leave blank

* Run 2SLS
ivregress 2sls `outcome' (`endog' = `instr') `controls', robust
estimates store ivmodel

* Run OLS
regress `outcome' `endog' `controls', robust
estimates store olsmodel

* Coefficient plot (coefplot recommended)
capture which coefplot
if _rc {
    display "coefplot not installed; install with: ssc install coefplot"
    coefplot olsmodel ivmodel, keep(`endog') xline(0) ciopts(recast(rcap)) title("OLS vs IV: `endog'")
}
else {
    coefplot olsmodel ivmodel, keep(`endog') xline(0) ciopts(recast(rcap)) title("OLS vs IV: `endog'")
}


* 2SRI for nonlinear outcome (Poisson example)
use predicted_LOS_by_author_menbreg.dta, clear
local outcome = "outcome_var"   // replace with your outcome
local endog   = "mean_LOS_MJ1"  // replace if different
local instr   = "instr_var"     // replace with your instrument
local controls = "cov1 cov2"    // replace with your covariates or leave blank

* First stage (linear)
regress `endog' `instr' `controls'
predict double endog_resid, resid
predict double endog_hat, xb

* Second stage (Poisson)
poisson `outcome' `endog' endog_resid `controls', vce(robust)
estimates store 2sri_poisson
display "Coefficient on endog_resid tests for endogeneity (significant => endogeneity present)"

* If outcome is binary, replace poisson with logit:
* logit `outcome' `endog' endog_resid `controls', vce(robust)






* First stage scatter and regression
use predicted_LOS_by_author_menbreg.dta, clear
* Replace endog and instr with your variable names if different
local endog "mean_LOS_MJ1"
local instr "instr_var"
twoway scatter `endog' `instr', msize(medium) mcolor(navy)
twoway lfit `endog' `instr', lcolor(maroon) lwidth(medium)
regress `endog' `instr'
display "Check the t and F for `instr' in the regress output (partial F should be large)"



* First-stage residual diagnostics
use predicted_LOS_by_author_menbreg.dta, clear
local endog "mean_LOS_MJ1"
local instr "instr_var"
regress `endog' `instr'
predict double endog_hat, xb
predict double endog_resid, resid
histogram endog_resid, normal width(0.5) title("Histogram of first-stage residuals")
twoway scatter endog_resid endog_hat, msize(small) mcolor(black)
twoway lfit endog_resid endog_hat, lcolor(gs12)



* 2SLS vs OLS comparison (linear outcome)
use predicted_LOS_by_author_menbreg.dta, clear
local outcome "outcome_var"
local endog "mean_LOS_MJ1"
local instr "instr_var"
local controls "cov1 cov2"
ivregress 2sls `outcome' (`endog' = `instr') `controls', robust
estimates store ivmodel
regress `outcome' `endog' `controls', robust
estimates store olsmodel
capture which coefplot
if _rc {
    display "coefplot not installed; install with: ssc install coefplot"
    coefplot olsmodel ivmodel, keep(`endog') xline(0) ciopts(recast(rcap)) title("OLS vs IV: `endog'")
}
else {
    coefplot olsmodel ivmodel, keep(`endog') xline(0) ciopts(recast(rcap)) title("OLS vs
	
	
	
* ---------------------------
* A. Save current Top-20 subset (so you can return to it)
* ---------------------------
tempfile top20_save
save `top20_save', replace
display "Top-20 subset saved to temp file: `top20_save'"

* ---------------------------
* B. Reload full dataset (replace path if your file is named differently)
* ---------------------------
use "path/to/your_full_dataset.dta", clear
display "Full dataset loaded."

* ---------------------------
* C. Sanity checks (warn if any required var missing)
* ---------------------------
capture confirm variable LOSCalc
if _rc display as error "Missing variable: LOSCalc"
capture confirm variable ICU
if _rc display as error "Missing variable: ICU"
capture confirm variable MYHighJeo
if _rc display as error "Missing variable: MYHighJeo"
capture confirm variable AuthorCode
if _rc display as error "Missing variable: AuthorCode"

* ---------------------------
* D. First stage: multilevel logit for ICU (melogit)
*    Instrument: MYHighJeo; random intercept by AuthorCode
* ---------------------------
melogit ICU MYHighJeo ib12.DisD Age sex i.zipCat ib12.Service i.DayofWeek ib4.MDC ib3.FClass_c ///
        HMProviderCount HMServiceCount NonHM1 NonICUNonHM || AuthorCode:
display "First-stage melogit complete. Inspect coefficient on MYHighJeo and its z-stat."

* Predicted conditional mean (subject-specific predicted probability)
predict pICU, mu

* First-stage residual for 2SRI
gen double rICU = ICU - pICU

* ---------------------------
* E. Second stage: menbreg with ICU and residual rICU (2SRI)
*    Include same controls and random intercept by AuthorCode; report incidence-rate ratios
* ---------------------------
menbreg LOSCalc ICU rICU MYHighJeo ib12.DisD Age sex i.zipCat ib12.Service i.DayofWeek ///
       ib4.MDC ib3.FClass_c HMProviderCount HMServiceCount NonHM1 NonICUNonHM || AuthorCode:, irr
estimates store menbreg_2sri
display "2SRI menbreg complete; estimates stored as menbreg_2sri."

* ---------------------------
* F. Diagnostics and quick checks
*    - Check relevance: coefficient on MYHighJeo in first stage
*    - Check endogeneity: significance of rICU in second stage
* ---------------------------
display "---- First-stage (melogit) coefficient on MYHighJeo ----"
melogit ICU MYHighJeo ib12.DisD Age sex i.zipCat ib12.Service i.DayofWeek ib4.MDC ib3.FClass_c ///
        HMProviderCount HMServiceCount NonHM1 NonICUNonHM || AuthorCode:
* (re-run to show output for inspection)

display "---- Test: is rICU significant in the menbreg (evidence of endogeneity)? ----"
test rICU

* Optional: ICC for melogit (between-author clustering)
estat icc

* ---------------------------
* G. If you want to compare Top-20 vs Full-sample estimates
* ---------------------------
* (Full-sample estimates already stored as menbreg_2sri)
use `top20_save', clear
display "Top-20 subset reloaded."
* Re-run the same menbreg on Top-20 (to store for comparison)
melogit ICU MYHighJeo ib12.DisD Age sex i.zipCat ib12.Service i.DayofWeek ib4.MDC ib3.FClass_c ///
        HMProviderCount HMServiceCount NonHM1 NonICUNonHM || AuthorCode:
predict pICU_top20, mu
gen double rICU_top20 = ICU - pICU_top20
menbreg LOSCalc ICU rICU_top20 MYHighJeo ib12.DisD Age sex i.zipCat ib12.Service i.DayofWeek ///
       ib4.MDC ib3.FClass_c HMProviderCount HMServiceCount NonHM1 NonICUNonHM || AuthorCode:, irr
estimates store menbreg_top20

* Compare stored estimates (if esttab/estout installed you can use that; otherwise use estimates table)
estimates table menbreg_2sri menbreg_top20, b se p

display "Done. Full-sample estimates stored as menbreg_2sri; Top-20 stored as menbreg_top20."




* Install coefplot if not already installed
ssc install coefplot, replace

* Replace reg_model and tsls_model with your stored estimate names
* This draws only the ICU coefficient; remove keep() to show more terms
coefplot (Regular = menbreg) (2SRI = menbreg_2sri),  keep(ICU)  drop(_cons)  xline(0, lpattern(dash) lcolor(black))  horizontal ciopts(recast(rcap) lwidth(medium))  msymbol(circle) msize(medium)  legend(order(1 "Regular" 2 "2SRI") ring(0) pos(6))    title("ICU treatment effect: Regular vs 2SRI") xlabel(, grid)

* Export publication PNG (1200x800 px)
graph export "forest_pub.png", width(1200) height(800) replace

* Export poster PNG (36x48 in at 300 dpi = 10800x14400 px)
* Note: very large exports may be slow or hit system limits; if Stata errors, try a smaller size or export vector from another tool.
graph export "forest_poster.png", width(10800) height(14400) replace

* Save the numeric table used for plotting
* coefplot can create a matrix of estimates; here we create a small table manually for reproducibility
tempfile forest_table

tempfile forest_table
postfile handle str20 model double estimate double se double lower double upper double p ///
    double irr double irr_lo double irr_hi using `forest_table', replace
foreach m in menbreg menbreg_2sri {
    estimates restore `m'
    scalar bICU  = _b[ICU]
    scalar seICU = _se[ICU]
    scalar lo    = bICU - invnormal(0.975)*seICU
    scalar hi    = bICU + invnormal(0.975)*seICU
    scalar pval  = 2*normal(-abs(bICU/seICU))      // z-based; menbreg has no e(df_r)
    post handle ("`m'") (bICU) (seICU) (lo) (hi) (pval) ///
        (exp(bICU)) (exp(lo)) (exp(hi))
}
postclose handle
use `forest_table', clear
list, noobs
* NOTE: these are the NAIVE model SEs/CIs. For the 2SRI model use the
* cluster-bootstrap CI from the Appendix instead of lower/upper here.
save "forest_table.dta", replace
export delimited using "forest_table.csv", replace




* -------------------------
* 1. Replace these placeholders with your actual variable names
* -------------------------
* outcome: LOS
* treatment: ICU            (binary indicator for treatment)
* instruments: Z1 Z2        (replace with your instrument variable(s))
* covariates: age male comorbidity_score othercov1 othercov2  (replace with your covariates)
* cluster variable: physician_id
* exposure (if used): exposure_var   (remove exposure() if not applicable)
* -------------------------

* 0. Confirm dataset loaded and working directory
display "Working directory: " c(pwd)
describe in 1/1

* 1. Naive model: multilevel negative binomial (store as menbreg)
menbreg LOS ICU age male comorbidity_score othercov1 othercov2 || physician_id:, cov(unstructured)
estimates store menbreg

* 2. First stage for 2SRI: linear probability model for ICU (recommended for residuals)
regress ICU Z1 Z2 age male comorbidity_score othercov1 othercov2, vce(cluster physician_id)
predict double ICU_hat, xb
gen double ICU_resid = ICU - ICU_hat

* 2b. (Optional alternative) If you prefer probit first stage, uncomment below:
* probit ICU Z1 Z2 age male comorbidity_score othercov1 othercov2, vce(cluster physician_id)
* predict double ICU_phat, pr
* gen double ICU_resid = ICU - ICU_phat

* 3. Second stage: include residual from first stage (2SRI)
menbreg LOS ICU ICU_resid age male comorbidity_score othercov1 othercov2 || physician_id: , vce(unstructured)
estimates store menbreg_2sri

* 4. Save model diagnostics (N and number of clusters)
scalar N_obs = e(N)
display "Observations in last model: " N_obs
* number of clusters (unique physician_id)
preserve
keep physician_id
bysort physician_id: keep if _n==1
count
scalar N_clusters = r(N)
restore
display "Physician clusters: " N_clusters

* 5. Extract all non-constant coefficients from both stored estimates into a dataset and CSV
tempfile forest_table_recreate
postfile ph str40 model str80 term double estimate double lower double upper double p using `forest_table_recreate', replace

* Extract from menbreg (naive)
estimates restore menbreg
matrix b = e(b)'
matrix V = e(V)
local coefnames : colnames b
foreach cname of local coefnames {
    if "`cname'" != "_cons" {
        scalar be = b[1, "`cname'"]
        scalar se = sqrt(V[`cname', "`cname'"])
        scalar lo = be - invnormal(0.975)*se
        scalar hi = be + invnormal(0.975)*se
        scalar pval = 2*ttail(e(df_r), abs(be/se))
        post ph ("Naive_menbreg") ("`cname'") (be) (lo) (hi) (pval)
    }
}

* Extract from menbreg_2sri
estimates restore menbreg_2sri

tempfile forest_table_recreate
postfile ph str40 model str80 term double estimate double lower double upper double p ///
    using `forest_table_recreate', replace
foreach m in menbreg menbreg_2sri {
    estimates restore `m'
    local eq = e(depvar)                 // main fixed-effects equation name
    matrix b = e(b)                      // no transpose
    matrix V = e(V)
    local fnames : colfullnames b        // "eq:term" for each column
    local j = 0
    foreach fn of local fnames {
        local ++j
        if !regexm("`fn'", "^`eq':") continue   // keep main equation only
        if  regexm("`fn'", ":_cons$") continue  // drop the intercept
        scalar se = sqrt(V[`j', `j'])
        if se == 0 | missing(se) continue       // drop omitted/base levels
        scalar be   = b[1, `j']
        scalar lo   = be - invnormal(0.975)*se
        scalar hi   = be + invnormal(0.975)*se
        scalar pval = 2*normal(-abs(be/se))
        local term  = subinstr("`fn'", "`eq':", "", 1)
        post ph ("`m'") ("`term'") (be) (lo) (hi) (pval)
    }
}
postclose ph
use `forest_table_recreate', clear

* 6. Tidy and export numeric table for PowerPoint
* Create a nicely formatted string for estimate (estimate (lower, upper))
gen str est_ci = string(estimate, "%9.3f") + " (" + string(lower, "%9.3f") + ", " + string(upper, "%9.3f") + ")"
order model term estimate lower upper p est_ci
sort term model

* Save outputs with unique filenames (change path if needed)
save "forest_table_recreated.dta", replace
export delimited using "forest_table_recreated.csv", replace

* 7. Display the table in Results for quick copy/paste
list model term est_ci p, noobs


* -------------------------
* Updated Steps 6 and 7
* -------------------------
local outdir "C:\Users\YourUser\Documents"    // <- change to a writable folder
cd "`outdir'"

* --- Step 6: create or update est_ci safely ---
capture confirm variable est_ci
if !_rc {
    drop est_ci
}
gen strL est_ci = string(estimate, "%9.3f") + " (" + string(lower, "%9.3f") + ", " + string(upper, "%9.3f") + ")"

* Ensure ordering and sorting for presentation
order model term estimate lower upper p est_ci
sort term model

* Create a timestamp token for unique filenames (use current date)
local datestr = c(current_date)

* Save numeric table and export CSV with unique name
save "`outdir'\forest_table_recreated_`datestr'.dta", replace
export delimited using "`outdir'\forest_table_recreated_`datestr'.csv", replace

* Display the table for quick copy/paste into PowerPoint
list model term est_ci p, noobs

* --- Also save a small summary with N and clusters (if scalars exist) ---
capture confirm scalar N_obs
if _rc {
    scalar N_obs = . 
}
capture confirm scalar N_clusters
if _rc {
    scalar N_clusters = .
}
clear
set obs 1
gen N = N_obs
gen clusters = N_clusters
save "`outdir'\forest_summary_recreated_`datestr'.dta", replace
export delimited using "`outdir'\forest_summary_recreated_`datestr'.csv", replace

* --- Step 7: build plotting positions and y-axis label pairs safely ---
use "`outdir'\forest_table_recreated_`datestr'.dta", clear

* Create a numeric id for each unique term (one center per term)
egen termid = group(term), label

* For stacking models within each term, create a sequence
bysort termid model: gen seq = _n

* Compute y position so models for the same term are adjacent and terms are separated
gen y = termid * 3 - seq

* Build ylab pairs (position "label") for twoway ylab()
preserve
keep term termid
bysort termid: keep if _n==1
gen ymid = termid * 3 - 1
local pairs ""
forvalues i = 1/`=_N' {
    local ypos = ymid[`i']
    local tname = term[`i']
    local pairs `pairs' `ypos' "`tname'"
}
restore

* Optional: choose colors for models (adjust as desired)
levelsof model, local(models)
local color1 "blue"
local color2 "red"
* Map model names to colors (assumes two models; extend if more)
local m1 : word 1 of `models'
local m2 : word 2 of `models'

* --- Draw the forest plot: horizontal CIs and points, color by model ---
twoway (rcap lower upper y if model=="`m1'", horizontal lwidth(medium) lcolor(black)) ///
       (scatter estimate y if model=="`m1'", msymbol(circle) msize(medium) mcolor(`color1')) ///
       (rcap lower upper y if model=="`m2'", horizontal lwidth(medium) lcolor(black)) ///
       (scatter estimate y if model=="`m2'", msymbol(circle) msize(medium) mcolor(`color2')), ///
       ylab(`pairs') xline(0, lpattern(dash) lcolor(black)) xlabel(, grid) ///
       title("Forest plot: menbreg vs menbreg_2sri (all coefficients)") ///
       legend(order(1 "`m1'" 2 "`m2'") ring(0) pos(6))

* Save graph object and export images with unique filenames
graph save "`outdir'\forest_graph_recreated_`datestr'.gph", replace
graph use "`outdir'\forest_graph_recreated_`datestr'.gph"
graph export "`outdir'\forest_pub_recreated_`datestr'.png", replace width(1200) height(800
)
graph use "`outdir'\forest_graph_recreated_`datestr'.gph"
graph export "`outdir'\forest_poster_recreated_`datestr'.png", replace width(3600) height(4800)

display "Step 6 and 7 complete. CSVs and graphs saved to `outdir'."




* 8. Also save a small summary with N and clusters
clear
set obs 1
gen N = N_obs
gen clusters = N_clusters
save "forest_summary_recreated.dta", replace
export delimited using "forest_summary_recreated.csv", replace

display "Done. CSV files: forest_table_recreated.csv and forest_summary_recreated.csv"


* Appendix A: CLUSTER BOOTSTRAP for 2SRI standard errors  (better for inference)

capture program drop boot2sri
program define boot2sri, rclass
    melogit ICU MYHighJeo <controls> || newid:
    predict double pICU, mu
    gen double rICU = ICU - pICU
    menbreg LOSCalc ICU rICU MYHighJeo <controls> || newid:
    return scalar bICU = _b[ICU]
    drop pICU rICU
end

bootstrap bICU = r(bICU), reps(500) cluster(AuthorCode) idcluster(newid) seed(1): boot2sri

* Appendix B: ESTTAB/ESTOUT to make tables

* ============================================================================
* esttab / estout alternative to the manual postfile forest tables
* (requires: ssc install estout, replace)
* ============================================================================
capture which esttab
if _rc ssc install estout, replace

* ---------------------------------------------------------------------------
* 1. FOCUSED TABLE  -- replaces the ICU-only postfile block.
*    Reports incidence-rate ratios (exp(b)) for the treatment and, for the
*    2SRI model, the first-stage residual (the endogeneity-test term).
* ---------------------------------------------------------------------------
* eform  -> exponentiate coefficients to IRR (menbreg is on the log scale)
* keep() -> show only the terms of interest; change rICU to ICU_resid if that
*           is your residual's name in this model
* ci(3)  -> 95% CI (also exponentiated) printed beneath each IRR, 3 decimals
* a CI that excludes 1 for rICU is the evidence-of-endogeneity check
esttab menbreg menbreg_2sri,                               ///
    eform b(3) ci(3)                                       ///
    keep(ICU rICU)                                         ///
    star(* 0.05 ** 0.01 *** 0.001)                         ///
    mtitles("Naive" "2SRI")                                ///
    title("ICU effect (IRR): naive vs 2SRI")               ///
    note("rICU = first-stage residual; IRR CI excluding 1 => endogeneity.")

* Same table to file (CSV for Excel, or swap to .rtf for Word):
esttab menbreg menbreg_2sri using "icu_irr_esttab.csv", replace ///
    eform b(3) ci(3) keep(ICU rICU) mtitles("Naive" "2SRI")

* ---------------------------------------------------------------------------
* 2. FULL COEFFICIENT TABLE  -- replaces the matrix-scraping extraction block.
*    esttab pulls every fixed effect automatically and, crucially, handles the
*    z-based p-values and CIs correctly for an ML estimator -- no hand-rolled
*    invnormal()/normal() arithmetic, and no risk of the e(b) transpose bug.
* ---------------------------------------------------------------------------
* keep() the model regressors so the dispersion (/lnalpha) and the random-
*   effect variance terms are excluded. Listing them explicitly is the robust
*   route because the variance-component's internal name is version-specific.
*   (Alternative: drop(_cons lnalpha:_cons <variance term>) -- run esttab once
*    with no keep/drop to see the exact ancillary names, then add them.)
* scalars() prints N and the number of clusters/groups under the table.
esttab menbreg menbreg_2sri,                               ///
    eform b(3) ci(3)                                       ///
    keep(ICU ICU_resid age male comorbidity_score othercov1 othercov2) ///
    star(* 0.05 ** 0.01 *** 0.001)                         ///
    label nonumbers mtitles("Naive" "2SRI")                ///
    scalars("N Observations" "N_g Clusters")               ///
    title("All fixed effects (IRR): naive vs 2SRI")

esttab menbreg menbreg_2sri using "forest_table_esttab.rtf", replace ///
    eform b(3) ci(3)                                       ///
    keep(ICU ICU_resid age male comorbidity_score othercov1 othercov2) ///
    label mtitles("Naive" "2SRI") scalars("N Observations" "N_g Clusters")

* ---------------------------------------------------------------------------
* 3A. (CAPSTONE) Report the VALID 2SRI interval from the bootstrap appendix.
*    The CIs in tables 1-2 are the naive model CIs -- too narrow for 2SRI.
*    The bootstrap command leaves an estimation result, so store it and let
*    esttab format it. eform turns the log coefficient into an IRR; the CI is
*    now the cluster-bootstrap one.
* ---------------------------------------------------------------------------
bootstrap bICU = r(bICU), reps(500) cluster(AuthorCode) idcluster(newid) seed(1): boot2sri
estimates store icu_boot
esttab icu_boot, eform b(3) ci(3) ///
    mtitles("2SRI (cluster bootstrap)") title("ICU IRR with valid CI")
* For a percentile (rather than normal-approx) interval, follow the bootstrap
* with:  estat bootstrap, percentile

* If you are instead running this after Appendix A above, rather than stand-alone
* the bootstrap already ran, so we just store + format it

* 3B. (ALTERNATE CAPSTONE)
estimates store icu_boot
esttab icu_boot, eform b(3) ci(3) ///
    mtitles("2SRI (cluster bootstrap)") title("ICU IRR with valid CI")
* estat bootstrap, percentile   // for a percentile interval