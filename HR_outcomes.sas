proc datasets library=work memtype=data kill nolist ; run ; quit ;

proc import datafile = "" out=stroke replace ;
run ;

* Strata: cerebrovasculardisease ;
data stroke ;
	set stroke ;
	if cerebrovasculardisease = 0 then strata = "Strata_1" ;
	else if cerebrovasculardisease = 1 then strata = "Strata_2" ;
run ;

* derive varibal;
data cathlab_stoke;
	set cathlab_stoke ;
	if calcification in (1, 2) then calc = 1 ;
	if calcification = 0 then calc = 0 ;
run ;

* derive varibal comp2_1y and tt_comp2_1y ;
data stroke;
	set stroke ;
	if death_1y = 1 or mi_1y = 1 or cva_1y = 1 then comp2_1y = 1 ;
	else comp2_1y = 0 ;
	tt_comp2_1y = MIN(tt_death_1y, tt_mi_1y, tt_cva_1y) ;
run ;

/* validation ;
proc freq data= stroke;
	table cerebrovasculardisease*comp2_1y/ norow nocol ;
run ;

proc means data =stroke ;
	var tt_comp2_1y ;
run ;
*/

proc sql noprint ;
SELECT COUNT(unitno)INTO: cat1
FROM cathlab_stoke
WHERE strata = "Strata_1" ;

SELECT COUNT(unitno)INTO: cat2
FROM cathlab_stoke
WHERE strata = "Strata_2" ;
quit ;

proc sort data=cathlab_stoke ; by strata ; run ;


%macro surv (outds=, vart=, vare=, label=) ;

* count number for specific survival outcome ;
proc freq data= NOPRINT ;
table strata*&vare / out=m_&outds ;
run ;

* derive failure rates ;
proc lifetest data = method=km intervals=(0 to 366) OUTSURV = rate;
	time &vart*&vare(0);
	strata strata ;
run; 

proc sort data= rate; where survival NE . ; by strata &vart ; run ;

data rate(keep = strata survival failure) ;
	set rate ;
	by strata &vart ;
	if last.strata ;
	failure = 1-survival ;
run ;

data m_&outds ;
	merge m_&outds rate ;
	by strata ;
	n = cats(put(count, 4.0))||" ("||cats(put(failure*100, 4.1))||")" ;
run ;

proc transpose data= m_&outds out= m_&outds(keep = &vare strata_1 strata_2) ;
	where &vare = 1 ;
	id strata ;
	var n ;
	by &vare ;
run ;

* unadjusted cox model ;
ods output ParameterEstimates = p_&outds(keep = parameter ProbChiSq HazardRatio HRLowerCL HRUpperCL label) ;
proc phreg data= cathlab_stoke ;
	class cerebrovasculardisease(ref = '0') ;
	model &vart*&vare(0) = cerebrovasculardisease / risklimits;
	*hazardratio cerebrovasculardisease / CL=WALD ;
run ;

* adjusted cox mode ;
ods output ParameterEstimates = p_&outds._adjust(keep = parameter ProbChiSq HazardRatio HRLowerCL HRUpperCL label) ;
proc phreg data= cathlab_stoke ;
	class cerebrovasculardisease(ref = '0') sex(ref = '0') race_eth ;
	model &vart*&vare(0) = cerebrovasculardisease age sex race_eth 
	peripheralarterialdisease atrialfibrillation anemia multivessel calc / risklimits;
	*hazardratio cerebrovasculardisease / CL=WALD ;
run ;
ods output close ;

data p_&outds._adjust ;
	set p_&outds._adjust ;
	where parameter = "cerebrovasculardisea" ;
	rename ProbChiSq= Probadjust HazardRatio= HRadjust HRLowerCL= lciadjust HRUpperCL= uciadjust ;
run ;

%macro add_lab(inds=, var=, label_=) ;
data &inds(drop = &var) ;
length parameter $200. label $200. ;
	set &inds ;
	parameter = "&var." ;
	label = &label_ ;
run ;
%mend add_lab ;

%add_lab(inds= m_&outds, var= &vare, label_= &label)  ;
%add_lab(inds= p_&outds, var= &vare, label_= &label)  ;
%add_lab(inds= p_&outds._adjust, var= &vare, label_= &label)  ;

proc sql noprint ;
create table &outds as 
	select *
	from m_&outds

	left join p_&outds
	on m_&outds..parameter = p_&outds..parameter

	left join p_&outds._adjust
	on m_&outds..parameter = p_&outds._adjust.parameter ;

quit ;

data &outds ;
length hr_ci ci_adjust $20. ;
	set &outds ;
	hr_ci = put(HazardRatio, 4.2)||" ("||put(HRLowerCL, 4.2)||"-"||put(HRUpperCL, 4.2)||")" ;
	hrci_adjust = put(HRadjust, 4.2)||" ("||put(lciadjust, 4.2)||"-"||put(uciadjust, 4.2)||")" ;
run ;


%mend surv ;

* append all sub-datasets ;
data surv_outcome(keep = parameter label strata_1 strata_2 hr_ci p hrci_adjust p_adjust) ;
	set  ;
	if ProbChiSq NE . then do ;
		if ProbChiSq < 0.001 then p = "<0.001" ;
		if ProbChiSq >= 0.001 then p = put(ProbChiSq, 5.3) ;
	end ;
	if Probadjust NE . then do ;
		if Probadjust < 0.001 then p_adjust = "<0.001" ;
		if Probadjust >= 0.001 then p_adjust = put(Probadjust, 5.3) ;
	end ;
run ;

* output rtf table ;
ods listing close;
ods escapechar='^';
options papersize=letter orientation=portrait nocenter 
LEFTMARGIN=0.25in RIGHTMARGIN=0.25in TOPMARGIN=0.25in BOTTOMMARGIN=0.25in
nodate nonumber;
ods rtf file="" style=journal bodytitle;

title bold j=center font= 'Times New Roman' height=12pt "Table" ;


proc report data= surv_outcome nowindows split="*"
style(report)=[frame=Void bordertopwidth=1pt bordertopcolor=black]

style(header)=[font_face="times new roman" font_size=11pt font_style=roman fontweight=bold] 

style(column)=[font_face="times new roman" font_size=10pt] ;

	column label strata_2 strata_1 hr_ci p hrci_adjust p_adjust ;
	define label/ '' style=[cellwidth=1.3in just=left];
	define strata_1/ "trt1*(N = *[89.8%])" style=[cellwidth=1.2in just=center]; * order interested group first ;
	define strata_2/ "trt2*(N= *[10.2%])" style=[cellwidth=1.25in just=center]; 
	define hr_ci/ "Unadjusted HR*(95% CI)" style=[cellwidth=1.15in just=center];
	define p/ 'P Value' style=[cellwidth=0.75in just=center];
	define hrci_adjust/ "Adjusted HR*(95% CI)" style=[cellwidth=1.15in just=center];
	define p_adjust/ "P Value" style=[cellwidth=0.75in just=center];

	compute label;
	if label= " xxx" then call define(_row_,'style','style=[borderbottomwidth=1pt borderbottomcolor=black]'); 
	endcomp;

	compute after / style=[just=l font_size=8pt font_face="times new roman" fontweight=medium] ;
 		line @1 'Models are adjusted for .';
 		line @1 '';
 	endcomp;

run ;

ods rtf close ;
