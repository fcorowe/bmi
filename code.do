************************************************** 
* Code to replicate the results reported in the paper: Explaining the widening distribution of Body Mass Index. Evidence from England, 2004-2013.* By: F. Rowe
**************************************************
** See paper for details about data; 
	**variable time = 1 is 2012/4, time = 0 is 2002/4.

clear all
set more off
set matsize 2000
set memory 2g
cd "../Data"
use "<data_file_name>"
compress

svyset [pweight= weight]

***********************************************************************
****** Dummy for ethnicity and physically active & labelling
**********************************************************************
	*Ethnicity
gen ethnicity2=0 if ethnicity!=.
replace ethnicity2=1 if ethnicity==1
replace ethnicity2=1 if ethnicity==2
replace ethnicity2=0 if ethnicity==3
replace ethnicity2=0 if ethnicity==4
	label var ethnicity2 `"Ethnicity"'
	label define ethnicity2 0 `"0 Nonwhite"', modify
	label define ethnicity2 1 `"White"', modify
	label values ethnicity2 ethnicity2	
	
	*65 and over
gen age65plus=1 if age>64
replace age65plus=0 if age<64

	*Three age band classification: 20-39, 40-59, 60+
gen age2039=0 if age!=.
replace age2039=1 if age>19 & age<40
gen age4059=0 if age!=.
replace age4059=1 if age>39 & age<60
gen age60plus=0 if age!=.
replace age60plus=1 if age>59
gen agebands=0
replace agebands=1 if age2039==1
replace agebands=2 if age4059==1
replace agebands=3 if age60plus==1
	
	*Physical activity
gen phyact=0 if physact1==0
replace phyact=1 if physact1==1

	*Labelling for Time
	label var time `"Ref. Year"'
	label define time 0 `"2002-4"', modify
	label define time 1 `"2012-4"', modify
	label values time time	

*Clean data
save "<new_data_file_name>", replace

*******************************
*** Table 1: Raw Differences in BMI
*******************************
clear all
set memory 3g
use "<new_data_file_name>"
compress

	* "raw" BMI difference at percentiles 10, 25, 50, 75 and 90;
	**s: start; f:final
	sort time
	_pctile bmi [pweight=weight] if time==0, p(10 25 50 75 90)
	return list
	gen s10=r(r1)
	gen s25=r(r2)
	gen s50=r(r3)
	gen s75=r(r4)
	gen s90=r(r5)
	_pctile bmi [pweight=weight] if time==1, p(10 25 50 75 90)
	return list
	gen f10=r(r1)
	gen f25=r(r2)
	gen f50=r(r3)
	gen f75=r(r4)
	gen f90=r(r5)

	*Difference between percentiles
	gen dif90=f90-s90
	gen dif75=f75-s75
	gen dif50=f50-s50
	gen dif25=f25-s25
	gen dif10=f10-s10
	di "  f13-s04 q10 =" dif10 "  f13-s04 q25 =" dif25 " f13-s04 q50 =" dif50 "  f13-s04 q75 ="  dif75 "  f13-s04 q90 ="  dif90
	
	format %7.6g bmi
    cendif bmi [pweight=weight], by(time) centile(10 25 50 75 90)
	
	*xtile quan = bmi [pweight=weight] if time==0, c(22 24 27 29 34)
	*xtile bmi quan [pweight=weight] if time==1, p(10 25 50 75 90)
	*ttest bmi [pweight=weight] if quan==, by(year) level(95)

	    
*******************************
*** Table 2. Variable weighted sample means at selected percentiles of the BMI distribution.
*******************************	

*< 25p

qui svy: mean male if bmi<s25, over(time)
matrix define meanest25 = e(b)
qui svy: mean male if bmi<s25, over(time)
lincom _b[male:_subpop_2]- _b[male:_subpop_1]
matrix define ci_lb = r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub = r(estimate) + invnorm(0.975)*r(se)

qui svy: mean age2039 if bmi<s25, over(time)
matrix define meanest25 = meanest25\e(b)
qui svy: mean age2039 if bmi<s25, over(time)
lincom _b[age2039:_subpop_2] - _b[age2039:_subpop_1]
matrix define ci_lb = ci_lb\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub = ci_ub\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean age4059 if bmi<s25, over(time)
matrix define meanest25 = meanest25\e(b)
qui svy: mean age4059 if bmi<s25, over(time)
lincom _b[age4059:_subpop_2] - _b[age4059:_subpop_1]
matrix define ci_lb = ci_lb\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub = ci_ub\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean age60plus if bmi<s25, over(time)
matrix define meanest25 = meanest25\e(b)
qui svy: mean age60plus if bmi<s25, over(time)
lincom _b[age60plus:_subpop_2] - _b[age60plus:_subpop_1]
matrix define ci_lb = ci_lb\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub = ci_ub\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean ethnicity2 if bmi<s25, over(time)
matrix define meanest25 = meanest25\e(b)
qui svy: mean ethnicity2 if bmi<s25, over(time)
lincom _b[ethnicity2:_subpop_2] - _b[ethnicity2:_subpop_1]
matrix define ci_lb = ci_lb\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub = ci_ub\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean phyact if bmi<s25, over(time)
matrix define meanest25 = meanest25\e(b)
qui svy: mean phyact if bmi<s25, over(time)
lincom _b[phyact:_subpop_2] - _b[phyact:_subpop_1]
matrix define ci_lb = ci_lb\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub = ci_ub\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean noqual if bmi<s25, over(time)
matrix define meanest25 = meanest25\e(b)
qui svy: mean noqual if bmi<s25, over(time)
lincom _b[noqual:_subpop_2] - _b[noqual:_subpop_1]
matrix define ci_lb = ci_lb\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub = ci_ub\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean secondary if bmi<s25, over(time)
matrix define meanest25 = meanest25\e(b)
qui svy: mean secondary if bmi<s25, over(time)
lincom _b[secondary:_subpop_2] - _b[secondary:_subpop_1]
matrix define ci_lb = ci_lb\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub = ci_ub\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean alevel if bmi<s25, over(time)
matrix define meanest25 = meanest25\e(b)
qui svy: mean alevel if bmi<s25, over(time)
lincom _b[alevel:_subpop_2] - _b[alevel:_subpop_1]
matrix define ci_lb = ci_lb\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub = ci_ub\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean degree if bmi<s25, over(time)
matrix define meanest25 = meanest25\e(b)
qui svy: mean degree if bmi<s25, over(time)
lincom _b[degree:_subpop_2] - _b[degree:_subpop_1]
matrix define ci_lb = ci_lb\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub = ci_ub\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean occlow if bmi<s25, over(time)
matrix define meanest25 = meanest25\e(b)
qui svy: mean occlow if bmi<s25, over(time)
lincom _b[occlow:_subpop_2] - _b[occlow:_subpop_1]
matrix define ci_lb = ci_lb\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub = ci_ub\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean occmed if bmi<s25, over(time)
matrix define meanest25 = meanest25\e(b)
qui svy: mean occmed if bmi<s25, over(time)
lincom _b[occmed:_subpop_2] - _b[occmed:_subpop_1]
matrix define ci_lb = ci_lb\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub = ci_ub\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean occhigh if bmi<s25, over(time)
matrix define meanest25 = meanest25\e(b)
qui svy: mean occhigh if bmi<s25, over(time)
lincom _b[occhigh:_subpop_2] - _b[occhigh:_subpop_1]
matrix define ci_lb = ci_lb\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub = ci_ub\r(estimate) + invnorm(0.975)*r(se)

matlist meanest25
matlist ci_lb
matlist ci_ub

	*>25p <75p
qui svy: mean male if bmi> s25 | bmi< s75, over(time)	
matrix define meanest50 = e(b)
qui svy: mean male if bmi> s25 | bmi< s75, over(time)
lincom _b[male:_subpop_2] - _b[male:_subpop_1]
matrix define ci_lb2 = r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub2 = r(estimate) + invnorm(0.975)*r(se)

qui svy: mean age2039 if bmi> s25 | bmi< s75, over(time)
matrix define meanest50 = meanest50\e(b)
qui svy: mean age2039 if bmi> s25 | bmi< s75, over(time)
lincom _b[age2039:_subpop_2] - _b[age2039:_subpop_1]
matrix define ci_lb2 = ci_lb2\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub2 = ci_ub2\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean age4059 if bmi> s25 | bmi< s75, over(time)
matrix define meanest50 = meanest50\e(b)
qui svy: mean age4059 if bmi> s25 | bmi< s75, over(time)
lincom _b[age4059:_subpop_2] - _b[age4059:_subpop_1]
matrix define ci_lb2 = ci_lb2\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub2 = ci_ub2\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean age60plus if bmi> s25 | bmi< s75, over(time)
matrix define meanest50 = meanest50\e(b)
qui svy: mean age60plus if bmi> s25 | bmi< s75, over(time)
lincom _b[age60plus:_subpop_2] - _b[age60plus:_subpop_1]
matrix define ci_lb2 = ci_lb2\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub2 = ci_ub2\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean ethnicity2 if bmi> s25 | bmi< s75, over(time)
matrix define meanest50 = meanest50\e(b)
qui svy: mean ethnicity2 if bmi> s25 | bmi< s75, over(time)
lincom _b[ethnicity2:_subpop_2] - _b[ethnicity2:_subpop_1]
matrix define ci_lb2 = ci_lb2\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub2 = ci_ub2\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean phyact if bmi> s25 | bmi< s75, over(time)
matrix define meanest50 = meanest50\e(b)
qui svy: mean phyact if bmi> s25 | bmi< s75, over(time)
lincom _b[phyact:_subpop_2] - _b[phyact:_subpop_1]
matrix define ci_lb2 = ci_lb2\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub2 = ci_ub2\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean noqual if bmi> s25 | bmi< s75, over(time)
matrix define meanest50 = meanest50\e(b)
qui svy: mean noqual if bmi> s25 | bmi< s75, over(time)
lincom _b[noqual:_subpop_2] - _b[noqual:_subpop_1]
matrix define ci_lb2 = ci_lb2\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub2 = ci_ub2\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean secondary if bmi> s25 | bmi< s75, over(time)
matrix define meanest50 = meanest50\e(b)
qui svy: mean secondary if bmi> s25 | bmi< s75, over(time)
lincom _b[secondary:_subpop_2] - _b[secondary:_subpop_1]
matrix define ci_lb2 = ci_lb2\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub2 = ci_ub2\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean alevel if bmi> s25 | bmi< s75, over(time)
matrix define meanest50 = meanest50\e(b)
qui svy: mean alevel if bmi> s25 | bmi< s75, over(time)
lincom _b[alevel:_subpop_2] - _b[alevel:_subpop_1]
matrix define ci_lb2 = ci_lb2\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub2 = ci_ub2\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean degree if bmi> s25 | bmi< s75, over(time)
matrix define meanest50 = meanest50\e(b)
qui svy: mean degree if bmi> s25 | bmi< s75, over(time)
lincom _b[degree:_subpop_2] - _b[degree:_subpop_1]
matrix define ci_lb2 = ci_lb2\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub2 = ci_ub2\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean occlow if bmi> s25 | bmi< s75, over(time)
matrix define meanest50 = meanest50\e(b)
qui svy: mean occlow if bmi> s25 | bmi< s75, over(time)
lincom _b[occlow:_subpop_2] - _b[occlow:_subpop_1]
matrix define ci_lb2 = ci_lb2\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub2 = ci_ub2\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean occmed if bmi> s25 | bmi< s75, over(time)
matrix define meanest50 = meanest50\e(b)
qui svy: mean occmed if bmi> s25 | bmi< s75, over(time)
lincom _b[occmed:_subpop_2] - _b[occmed:_subpop_1]
matrix define ci_lb2 = ci_lb2\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub2 = ci_ub2\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean occhigh if bmi> s25 | bmi< s75, over(time)
matrix define meanest50 = meanest50\e(b)
qui svy: mean occhigh if bmi> s25 | bmi< s75, over(time)
lincom _b[occhigh:_subpop_2] - _b[occhigh:_subpop_1]
matrix define ci_lb2 = ci_lb2\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub2 = ci_ub2\r(estimate) + invnorm(0.975)*r(se)
 
matlist meanest50
matlist ci_lb2
matlist ci_ub2


	*>75p
qui svy: mean male if bmi> s75, over(time)	
matrix define meanest75 = e(b)
qui svy: mean male if bmi> s75, over(time)
lincom _b[male:_subpop_2] - _b[male:_subpop_1]
matrix define ci_lb3 = r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub3 = r(estimate) + invnorm(0.975)*r(se)

qui svy: mean age2039 if bmi > s75, over(time)
matrix define meanest75 = meanest75\e(b)
qui svy: mean age2039 if bmi>s75, over(time)
lincom _b[age2039:_subpop_2] - _b[age2039:_subpop_1]
matrix define ci_lb3 = ci_lb3\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub3 = ci_ub3\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean age4059 if bmi>s75, over(time)
matrix define meanest75 = meanest75\e(b)
qui svy: mean age4059 if bmi>s75, over(time)
lincom _b[age4059:_subpop_2] - _b[age4059:_subpop_1]
matrix define ci_lb3 = ci_lb3\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub3 = ci_ub3\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean age60plus if bmi> s75, over(time)
matrix define meanest75 = meanest75\e(b)
qui svy: mean age60plus if bmi> s75, over(time)
lincom _b[age60plus:_subpop_2] - _b[age60plus:_subpop_1]
matrix define ci_lb3 = ci_lb3\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub3 = ci_ub3\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean ethnicity2 if bmi> s75, over(time)
matrix define meanest75 = meanest75\e(b)
qui svy: mean ethnicity2 if bmi> s75, over(time)
lincom _b[ethnicity2:_subpop_2] - _b[ethnicity2:_subpop_1]
matrix define ci_lb3 = ci_lb3\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub3 = ci_ub3\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean phyact if bmi>s75, over(time)
matrix define meanest75 = meanest75\e(b)
qui svy: mean phyact if bmi>s75, over(time)
lincom _b[phyact:_subpop_2] - _b[phyact:_subpop_1]
matrix define ci_lb3 = ci_lb3\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub3 = ci_ub3\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean noqual if bmi>s75, over(time)
matrix define meanest75 = meanest75\e(b)
qui svy: mean noqual if bmi>s75, over(time)
lincom _b[noqual:_subpop_2] - _b[noqual:_subpop_1]
matrix define ci_lb3 = ci_lb3\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub3 = ci_ub3\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean secondary if bmi>s75, over(time)
matrix define meanest75 = meanest75\e(b)
qui svy: mean secondary if bmi>s75, over(time)
lincom _b[secondary:_subpop_2] - _b[secondary:_subpop_1]
matrix define ci_lb3 = ci_lb3\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub3 = ci_ub3\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean alevel if bmi>s75, over(time)
matrix define meanest75 = meanest75\e(b)
qui svy: mean alevel if bmi>s75, over(time)
lincom _b[alevel:_subpop_2] - _b[alevel:_subpop_1]
matrix define ci_lb3 = ci_lb3\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub3 = ci_ub3\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean degree if bmi>s75, over(time)
matrix define meanest75 = meanest75\e(b)
qui svy: mean degree if bmi>s75, over(time)
lincom _b[degree:_subpop_2] - _b[degree:_subpop_1]
matrix define ci_lb3 = ci_lb3\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub3 = ci_ub3\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean occlow if bmi>s75, over(time)
matrix define meanest75 = meanest75\e(b)
qui svy: mean occlow if bmi>s75, over(time)
lincom _b[occlow:_subpop_2] - _b[occlow:_subpop_1]
matrix define ci_lb3 = ci_lb3\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub3 = ci_ub3\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean occmed if bmi>s75, over(time)
matrix define meanest75 = meanest75\e(b)
qui svy: mean occmed if bmi>s75, over(time)
lincom _b[occmed:_subpop_2] - _b[occmed:_subpop_1]
matrix define ci_lb3 = ci_lb3\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub3 = ci_ub3\r(estimate) + invnorm(0.975)*r(se)

qui svy: mean occhigh if bmi>s75, over(time)
matrix define meanest75 = meanest75\e(b)
qui svy: mean occhigh if bmi>s75, over(time)
lincom _b[occhigh:_subpop_2] - _b[occhigh:_subpop_1]
matrix define ci_lb3 = ci_lb3\r(estimate) + invnorm(0.025)*r(se)
matrix define ci_ub3 = ci_ub3\r(estimate) + invnorm(0.975)*r(se)

matlist meanest75
matlist ci_lb3
matlist ci_ub3


*****************************************************************************
******* Table 3. Unconditional Quantile Regressions exploring Body Mass Index, 2002-4 and 2012-4
******************************************************************************
clear all
set more off
set matsize 2000
cd "../Data"
use "<new_data_file_name>"
compress
				
forvalues qt = 10(40)90 {	
   gen rif_`qt'=.
}
	*Density for time=0
pctile eval1 = bmi if time == 1 , nq(100) 
kdensity bmi if time == 1, at(eval1) gen(evalf densf) width(0.2) nograph 
forvalues qt = 10(40)90 {	
 local qc = `qt'/100.0
 replace rif_`qt'= evalf[`qt'] + `qc'/ densf[`qt'] if bmi >= evalf[`qt'] & time==1
 replace rif_`qt'= evalf[`qt'] - (1-`qc') / densf[`qt'] if bmi < evalf[`qt'] & time==1
}

pctile eval2=bmi if time==0, nq(100) 
kdensity bmi if time==0, at(eval2) gen(evals denss) width(0.2) nograph 
forvalues qt = 10(40)90 {	
 local qc = `qt'/100.0
 replace rif_`qt'=evals[`qt'] + `qc'/denss[`qt'] if bmi>=evals[`qt'] & time==0
 replace rif_`qt'=evals[`qt'] - (1-`qc')/denss[`qt'] if bmi<evals[`qt'] & time==0
}

gen baseyear=0
replace baseyear=1 if time==0
sort baseyear
by baseyear: sum rif_10 rif_50 rif_90
svy: mean rif_10 rif_50 rif_90, over(baseyear)

 * RIF regressions
				reg rif_10 male age4059 age60plus ethnicity2 phyact secondary alevel degree occmed occhigh [pw = weight] if baseyear==1, robust
				estimates store UQR10_t1	
				reg rif_50 male age4059 age60plus ethnicity2 phyact secondary alevel degree occmed occhigh [pw = weight] if baseyear==1, robust
				estimates store UQR50_t1	
				reg rif_90 male age4059 age60plus ethnicity2 phyact secondary alevel degree occmed occhigh [pw = weight] if baseyear==1, robust
				estimates store UQR90_t1	
				reg rif_10 male age4059 age60plus ethnicity2 phyact secondary alevel degree occmed occhigh [pw = weight] if baseyear==0, robust
				estimates store UQR10_t2	
				reg rif_50 male age4059 age60plus ethnicity2 phyact secondary alevel degree occmed occhigh [pw = weight] if baseyear==0, robust
				estimates store UQR50_t2	
				reg rif_90 male age4059 age60plus ethnicity2 phyact secondary alevel degree occmed occhigh [pw = weight] if baseyear==0, robust
				estimates store UQR90_t2
				esttab UQR10_t1 UQR50_t1 UQR90_t1 UQR10_t2 UQR50_t2 UQR90_t2 /*
					*/ using "../Results/UQR_V1.rtf", b(3) se(3) nogaps star(* 0.10 ** 0.05 *** 0.01)/*
					*/ title("Unconditional Quantile Regressions, 2002-4 and 2012-4") /*
					*/ mtitle("QR=0.1" "QR=0.5" "QR=0.9" "QR=0.1" "QR=0.5" "QR=0.9") /*
					*/ refcat(male "{\i Gender}" age4059 "{\i Age}" ethnicity2 "{\i Ethnicity}" fv12 "{\i Diet}" /*
					*/ phyact "{\i Physical activity}" secondary "{\i Education}" occmed "{\i Occupation}", label(" ")) /*
					*/ coeflabels(male "Male" age4059 "Age 40-59" age60plus "Age 60+" ethnicity2 "White" /*
					*/ phyact "Moderate exercise" secondary "Secondary" alevel "A level" degree "Degree level" /*
					*/ occmed "Medium" occhigh "High" _cons "Constant") addnotes("Reference category: Age 20-39, Non-white, No physical activity, No qualification, Low-skilled occupation") /*
					*/ compress label nodepvar replace


**************************************************
******* Table 4. Table 4. Quantile decomposition of changes in Body Mass Index (BMI) between 2002-4 and 2012-4.
**************************************************
clear all
set more off
set matsize 2000
set memory 2g
cd "../Data"
use "<new_data_file_name>"
compress
				
forvalues qt = 10(40)90 {	
   gen rif_`qt'=.
}
	*Density for time=0
pctile eval1 = bmi if time == 1 , nq(100) 
kdensity bmi if time == 1, at(eval1) gen(evalf densf) width(0.2) nograph 
forvalues qt = 10(40)90 {	
 local qc = `qt'/100.0
 replace rif_`qt'= evalf[`qt'] + `qc'/ densf[`qt'] if bmi >= evalf[`qt'] & time==1
 replace rif_`qt'= evalf[`qt'] - (1-`qc') / densf[`qt'] if bmi < evalf[`qt'] & time==1
}

pctile eval2=bmi if time==0, nq(100) 
kdensity bmi if time==0, at(eval2) gen(evals denss) width(0.2) nograph 
forvalues qt = 10(40)90 {	
 local qc = `qt'/100.0
 replace rif_`qt'=evals[`qt'] + `qc'/denss[`qt'] if bmi>=evals[`qt'] & time==0
 replace rif_`qt'=evals[`qt'] - (1-`qc')/denss[`qt'] if bmi<evals[`qt'] & time==0
}

gen baseyear=0
replace baseyear=1 if time==0
sort baseyear
by baseyear: sum rif_10 rif_50 rif_90
svy: mean rif_10 rif_50 rif_90, over(baseyear)

		*Raw Differences
			qreg bmi time [aweight = weight], quantile(10)
			qreg bmi time [aweight = weight], quantile(50)
			qreg bmi time [aweight = weight], quantile(90)
				
		* Machado-MataÐMelly decomposition
			*net install counterfactual, from("https://sites.google.com/site/blaisemelly/home")
			rqdeco bmi male age4059 age60plus ethnicity2 phyact secondary alevel degree occmed occhigh [pw = weight], by(time) quantile(0.1 0.5 0.9) vce(boot) reps(100)
				
		* RIF decomposition		
		oaxaca rif_10 male normalize(age2039-age60plus) ethnicity2 phyact normalize(degree-noqual) /*
			*/ normalize(occhigh-occlow), by(baseyear) weight(1) detail(groupgen: male, groupage: age2039-age60plus, /*
			*/ groupeth: ethnicity2, grouppact: phyact, groupses: degree-noqual occhigh-occlow)
		estimates store RIF10
		oaxaca rif_50 male normalize(age2039-age60plus) ethnicity2 phyact normalize(degree-noqual) /*
			*/ normalize(occhigh-occlow), by(baseyear) weight(1) detail(groupgen: male, groupage: age2039-age60plus, /*
			*/ groupeth: ethnicity2, grouppact: phyact, groupses: degree-noqual occhigh-occlow) 
		estimates store RIF50
		oaxaca rif_90 male normalize(age2039-age60plus) ethnicity2 phyact normalize(degree-noqual) /*
			*/ normalize(occhigh-occlow), by(baseyear) weight(1) detail(groupgen: male, groupage: age2039-age60plus, /*
			*/ groupeth: ethnicity2, grouppact: phyact, groupses: degree-noqual occhigh-occlow)
		estimates store RIF90

		esttab RIF10 RIF50 RIF90 /*
		*/ using "../Results/RIF_V1.rtf", b(3) se(3) nogaps star(* 0.10 ** 0.05 *** 0.01)/*
		*/ title("Decomposition method: RIF regression, 2002-4 vs. 2012-4") /*
		*/ mtitle("10th Percentile" "50th Percentile" "90th Percentile") /*
		*/ coeflabels(groupgen "Gender" groupage "Age" groupeth "Ethnicity" /*
		*/ grouppact "Physical activity" groupses "Socio-economic" /*
		*/ unexplained "Exp. by coefficients" explained "Exp. by characteristics" /*
		*/ group_1 "2012-4" group_2 "2002-4" difference "Difference" _cons "Constant") compress label nodepvar replace


**************************************************
******* Appendix B. Unconditional Quantile Regressions by Sex
**************************************************

clear all
set more off
set matsize 2000
cd "../Data"
use "<new_data_file_name>"
compress
				
forvalues qt = 10(40)90 {	
   gen rif_`qt'=.
}
	*Density for time=0
pctile eval1 = bmi if time == 1 , nq(100) 
kdensity bmi if time == 1, at(eval1) gen(evalf densf) width(0.2) nograph 
forvalues qt = 10(40)90 {	
 local qc = `qt'/100.0
 replace rif_`qt'= evalf[`qt'] + `qc'/ densf[`qt'] if bmi >= evalf[`qt'] & time==1
 replace rif_`qt'= evalf[`qt'] - (1-`qc') / densf[`qt'] if bmi < evalf[`qt'] & time==1
}

pctile eval2=bmi if time==0, nq(100) 
kdensity bmi if time==0, at(eval2) gen(evals denss) width(0.2) nograph 
forvalues qt = 10(40)90 {	
 local qc = `qt'/100.0
 replace rif_`qt'=evals[`qt'] + `qc'/denss[`qt'] if bmi>=evals[`qt'] & time==0
 replace rif_`qt'=evals[`qt'] - (1-`qc')/denss[`qt'] if bmi<evals[`qt'] & time==0
}

gen baseyear=0
replace baseyear=1 if time==0
sort baseyear
by baseyear: sum rif_10 rif_50 rif_90
svy: mean rif_10 rif_50 rif_90, over(baseyear)

			* Table B.1. Unconditional Quantile Regressions, Male, 2002-4 and 2012-4
				reg rif_10 age4059 age60plus ethnicity2 phyact secondary alevel degree occmed occhigh [pw = weight] if baseyear==1 & male==1, robust
				estimates store M_UQR10_04	
				reg rif_50 age4059 age60plus ethnicity2 phyact secondary alevel degree occmed occhigh [pw = weight] if baseyear==1 & male==1, robust
				estimates store M_UQR50_04	
				reg rif_90 age4059 age60plus ethnicity2 phyact secondary alevel degree occmed occhigh [pw = weight] if baseyear==1 & male==1, robust
				estimates store M_UQR90_04	
				reg rif_10 age4059 age60plus ethnicity2 phyact secondary alevel degree occmed occhigh [pw = weight] if baseyear==0 & male==1, robust
				estimates store M_UQR10_13	
				reg rif_50 age4059 age60plus ethnicity2 phyact secondary alevel degree occmed occhigh [pw = weight] if baseyear==0 & male==1, robust
				estimates store M_UQR50_13	
				reg rif_90 age4059 age60plus ethnicity2 phyact secondary alevel degree occmed occhigh [pw = weight] if baseyear==0 & male==1, robust
				estimates store M_UQR90_13
				esttab M_UQR10_04  M_UQR50_04 M_UQR90_04 M_UQR10_13 M_UQR50_13 M_UQR90_13 /*
					*/ using "../Results/Male_UQR_V1.rtf", b(3) se(3) nogaps star(* 0.10 ** 0.05 *** 0.01)/*
					*/ title(Unconditional Quantile Regressions, Male, 2004 and 2013) /*
					*/ mtitle("QR=0.1" "QR=0.5" "QR=0.9" "QR=0.1" "QR=0.5" "QR=0.9") /*
					*/ refcat(male "{\i Gender}" age4059 "{\i Age}" ethnicity2 "{\i Ethnicity}" /*
					*/ phyact "{\i Physical activity}" secondary "{\i Education}" occmed "{\i Occupation}", label(" ")) /*
					*/ coeflabels(male "Male" age4059 "Age 40-59" age60plus "Age 60+" ethnicity2 "White" /*
					*/ fv34 "3-4 Fruit & Veg." fv5p "5+ Fruit & Veg." phyact "Moderate exercise" secondary "Secondary" alevel "A level" degree "Degree level" /*
					*/ occmed "Medium" occhigh "High" _cons "Constant") addnotes("Reference category: Age 20-39, Non-white, No physical activity, No qualification, Low-skilled occupation") /*
					*/ compress label nodepvar replace

			* Table B.2. Unconditional Quantile Regressions, Female, 2002-4 and 2012-4
				
				reg rif_10 age4059 age60plus ethnicity2 phyact secondary alevel degree occmed occhigh [pw = weight] if baseyear==1 & male==0, robust
				estimates store F_UQR10_04	
				reg rif_50 age4059 age60plus ethnicity2 phyact secondary alevel degree occmed occhigh [pw = weight] if baseyear==1 & male==0, robust
				estimates store F_UQR50_04	
				reg rif_90 age4059 age60plus ethnicity2 phyact secondary alevel degree occmed occhigh [pw = weight] if baseyear==1 & male==0, robust
				estimates store F_UQR90_04	
				reg rif_10 age4059 age60plus ethnicity2 phyact secondary alevel degree occmed occhigh [pw = weight] if baseyear==0 & male==0, robust
				estimates store F_UQR10_13	
				reg rif_50 age4059 age60plus ethnicity2 phyact secondary alevel degree occmed occhigh [pw = weight] if baseyear==0 & male==0, robust
				estimates store F_UQR50_13	
				reg rif_90 age4059 age60plus ethnicity2 phyact secondary alevel degree occmed occhigh [pw = weight] if baseyear==0 & male==0, robust
				estimates store F_UQR90_13
				esttab F_UQR10_04 F_UQR50_04 F_UQR90_04 F_UQR10_13 F_UQR50_13  F_UQR90_13 /*
					*/ using "../Results/Female_UQR_V1.rtf", b(3) se(3) nogaps star(* 0.10 ** 0.05 *** 0.01)/*
					*/ title("Unconditional Quantile Regressions, Female, 2004 and 2013") /*
					*/ mtitle("QR=0.1" "QR=0.5" "QR=0.9" "QR=0.1" "QR=0.5" "QR=0.9") /*
					*/ refcat(male "{\i Gender}" age4059 "{\i Age}" ethnicity2 "{\i Ethnicity}" /*
					*/ phyact "{\i Physical activity}" secondary "{\i Education}" occmed "{\i Occupation}", label(" ")) /*
					*/ coeflabels(male "Male" age4059 "Age 40-59" age60plus "Age 60+" ethnicity2 "White" /*
					*/ fv34 "3-4 Fruit & Veg." fv5p "5+ Fruit & Veg." phyact "Moderate exercise" secondary "Secondary" alevel "A level" degree "Degree level" /*
					*/ occmed "Medium" occhigh "High" _cons "Constant") addnotes("Reference category: Age 20-39, Non-white, No physical activity, No qualification, Low-skilled occupation") /*
					*/ compress label nodepvar replace



**************************************************
******* Appendix C. Decomposition by Sex
**************************************************
* Table C.1. Decomposition method: RIF regression by Sex, 2002-4 vs. 2012-4
		
			* Male
		oaxaca rif_10 age4059 age60plus ethnicity2 phyact secondary alevel degree occmed occhigh if male==1, by(baseyear) weight(1) detail(groupage: age4059 age60plus, groupeth: ethnicity2, grouppact: phyact, groupses: secondary alevel degree occmed occhigh)
		estimates store M_RIF10
		oaxaca rif_50 age4059 age60plus ethnicity2 phyact secondary alevel degree occmed occhigh if male==1, by(baseyear) weight(1) detail(groupage: age4059 age60plus, groupeth: ethnicity2, grouppact: phyact, groupses: secondary alevel degree occmed occhigh)
		estimates store M_RIF50
		oaxaca rif_90 age4059 age60plus ethnicity2 phyact secondary alevel degree occmed occhigh if male==1, by(baseyear) weight(1) detail(groupage: age4059 age60plus, groupeth: ethnicity2, grouppact: phyact, groupses: secondary alevel degree occmed occhigh)
		estimates store M_RIF90	

			*Female
		
		oaxaca rif_10 age4059 age60plus ethnicity2 phyact secondary alevel degree occmed occhigh if male==0, by(baseyear) weight(1) detail(groupage: age4059 age60plus, groupeth: ethnicity2, grouppact: phyact, groupses: secondary alevel degree occmed occhigh)
		estimates store F_RIF10
		oaxaca rif_50 age4059 age60plus ethnicity2 phyact secondary alevel degree occmed occhigh if male==0, by(baseyear) weight(1) detail(groupage: age4059 age60plus, groupeth: ethnicity2, grouppact: phyact, groupses: secondary alevel degree occmed occhigh)
		estimates store F_RIF50
		oaxaca rif_90 age4059 age60plus ethnicity2 phyact secondary alevel degree occmed occhigh if male==0, by(baseyear) weight(1) detail(groupage: age4059 age60plus, groupeth: ethnicity2, grouppact: phyact, groupses: secondary alevel degree occmed occhigh)
		estimates store F_RIF90
		
		esttab M_RIF10 M_RIF50 M_RIF90 F_RIF10 F_RIF50 F_RIF90/*
		*/ using "../Results/Sex_RIF_V1.rtf", b(3) se(3) nogaps star(* 0.10 ** 0.05 *** 0.01)/*
		*/ title("Decomposition method: RIF regression by Sex, 2002-4 vs. 2012-4") /*
		*/ mtitle("10th Percentile" "50th Percentile" "90th Percentile" "10th Percentile" "50th Percentile" "90th Percentile") /*
		*/ coeflabels(groupgen "Gender" groupage "Age" groupeth "Ethnicity" /*
		*/ grouppact "Physical activity" groupses "Socio-economic" /*
		*/ unexplained "Exp. by coefficients" explained "Exp. by characteristics" /*
		*/ group_1 "2012-4" group_2 "2002-4" difference "Difference" _cons "Constant") compress label nodepvar replace
		
				
