suppressPackageStartupMessages({library(data.table);library(WeightIt);library(cobalt);library(ggplot2);library(sandwich)})
options(warn=1)
try(Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8'),silent=TRUE)
set.seed(20260905)
root <- Sys.getenv('ANALYSIS_ROOT','..'); raw <- Sys.getenv('ANALYSIS_RAW',file.path(root,'outputs','restricted_patient_level')); main <- Sys.getenv('ANALYSIS_MAIN',file.path(root,'outputs','aggregate_review','main')); supp <- Sys.getenv('ANALYSIS_SUPP',file.path(root,'outputs','aggregate_review','supplementary')); env_dir <- Sys.getenv('ANALYSIS_ENV',file.path(root,'outputs','aggregate_review','environment')); internal_fig <- Sys.getenv('ANALYSIS_INTERNAL_FIG',file.path(root,'outputs','aggregate_review','internal'))
dir.create(env_dir,recursive=TRUE,showWarnings=FALSE); dir.create(internal_fig,recursive=TRUE,showWarnings=FALSE)
d <- fread(file.path(raw,'baseline_covariates.csv'))
stopifnot(!anyDuplicated(d$subject_id),all(d$ivme_nci_48h>=0))
d[, ivme_group:=factor(ivme_group,levels=paste0('G',1:4))]
d[, procedure_class:=fcase(grepl('aortic',procedure_groups),'Aortic_any',
 grepl('CABG',procedure_groups)&grepl('valve',procedure_groups),'CABG_valve',
 procedure_groups=='CABG','CABG_only',procedure_groups=='valve','Valve_only',default='Other')]
factors <- c('gender','race_group','admission_group','procedure_class','anchor_year_group')
d[, (factors):=lapply(.SD,factor),.SDcols=factors]
hist <- grep('^prior_',names(d),value=TRUE)
rhs <- paste(c('age_at_admission','I(age_at_admission^2)',factors,hist),collapse=' + ')
psf <- as.formula(paste('ivme_group ~',rhs))
balance_formula <- psf
theme_set(theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face='bold'),legend.position='bottom'))
savefig <- function(p,name,folder=main,w=8,h=5) {
 ggsave(file.path(folder,paste0(name,'.png')),p,width=w,height=h,dpi=300,bg='white')
 ggsave(file.path(folder,paste0(name,'.pdf')),p,width=w,height=h,device=cairo_pdf)
}
labels <- c(G1='0',G2='>0-10',G3='>10-25',G4='>25')
calibrate_overlap <- function(dat,w0,form) {
  mm<-model.matrix(delete.response(terms(form)),dat)[,-1,drop=FALSE]
  keep<-apply(mm,2,sd)>1e-10;mm<-mm[,keep,drop=FALSE]
  mm<-scale(mm);qq<-qr(mm);mm<-mm[,qq$pivot[seq_len(qq$rank)],drop=FALSE]
  target<-colSums(w0*mm)/sum(w0);wc<-numeric(nrow(dat));audit<-list()
  for(g in levels(dat$ivme_group)){
    i<-which(dat$ivme_group==g);X<-mm[i,,drop=FALSE];b<-w0[i]
    fn<-function(l){u<-log(b)+as.vector(X%*%l);mx<-max(u);log(sum(exp(u-mx)))+mx-sum(target*l)}
    gr<-function(l){u<-log(b)+as.vector(X%*%l);p<-exp(u-max(u));p<-p/sum(p);as.vector(crossprod(X,p)-target)}
    oo<-optim(rep(0,ncol(X)),fn,gr,method='BFGS',control=list(maxit=2000,reltol=1e-11))
    if(oo$convergence!=0 || max(abs(gr(oo$par)))>1e-4) stop('Calibration failed for ',g)
    u<-log(b)+as.vector(X%*%oo$par);ww<-exp(u-max(u));wc[i]<-ww/sum(ww)
    audit[[g]]<-data.table(group=g,convergence=oo$convergence,max_moment_error=max(abs(gr(oo$par))))
  }
  list(weights=wc,audit=rbindlist(audit))
}
fitow <- function(dat,form=psf) {
 w <- weightit(form,data=dat,method='glm',estimand='ATO',multi.method='weightit',include.obj=TRUE)
 stopifnot(all(is.finite(w$weights)),all(w$weights>0))
 cc<-calibrate_overlap(dat,w$weights,form)
 m <- suppressWarnings(glm(incident_cam_positive~ivme_group,data=dat,family=binomial(),weights=cc$weights))
 stopifnot(m$converged,all(is.finite(coef(m))))
 list(w=w,calibrated_weights=cc$weights,calibration_audit=cc$audit,m=m,vc=vcovHC(m,type='HC0'))
}
summarize_fit <- function(obj,dat,name) {
 m<-obj$m; V<-obj$vc; b<-coef(m)
 X<-model.matrix(~ivme_group,data.frame(ivme_group=factor(paste0('G',1:4),levels=paste0('G',1:4))))
 eta<-as.vector(X%*%b); p<-plogis(eta); se<-sqrt(diag(X%*%V%*%t(X)))
 risks<-data.table(analysis=name,group=paste0('G',1:4),risk=p,lo=plogis(eta-1.96*se),hi=plogis(eta+1.96*se))
 pairs<-list(c(2,1),c(3,1),c(4,1),c(4,2))
 effects<-rbindlist(lapply(pairs,function(ij) {
  i<-ij[1];j<-ij[2];cvec<-X[i,]-X[j,];beta<-sum(cvec*b);s<-sqrt(as.numeric(t(cvec)%*%V%*%cvec))
  g<-p[i]*(1-p[i])*X[i,]-p[j]*(1-p[j])*X[j,];sr<-sqrt(as.numeric(t(g)%*%V%*%g))
  data.table(analysis=name,comparison=paste0('G',i,' vs G',j),OR=exp(beta),OR_lo=exp(beta-1.96*s),OR_hi=exp(beta+1.96*s),
   RD=p[i]-p[j],RD_lo=p[i]-p[j]-1.96*sr,RD_hi=p[i]-p[j]+1.96*sr,p_value=2*pnorm(-abs(beta/s)),n=nrow(dat),events=sum(dat$incident_cam_positive))
 }))
 list(risks=risks,effects=effects)
}
o<-fitow(d); d[,overlap_weight_raw:=o$w$weights];d[,overlap_weight:=o$calibrated_weights]
s<-summarize_fit(o,d,'Primary four-group ATO')
# Nonparametric bootstrap repeats the propensity model, calibration, and outcome estimation.
boot_one<-function(seed){
 set.seed(seed);bd<-d[sample.int(nrow(d),replace=TRUE)];bo<-try(fitow(bd),silent=TRUE)
 if(inherits(bo,'try-error'))return(NULL)
 bd[,bootstrap_weight:=bo$calibrated_weights]
 rr<-bd[,.(risk=weighted.mean(incident_cam_positive,bootstrap_weight)),by=ivme_group][order(ivme_group)]$risk
 if(length(rr)!=4||any(rr<=0|rr>=1))return(NULL)
 pp<-list(c(2,1),c(3,1),c(4,1),c(4,2))
 rbindlist(lapply(pp,function(ij)data.table(comparison=paste0('G',ij[1],' vs G',ij[2]),OR=(rr[ij[1]]/(1-rr[ij[1]]))/(rr[ij[2]]/(1-rr[ij[2]])),RD=rr[ij[1]]-rr[ij[2]])))[,`:=`(replicate=seed,risk_G1=rr[1],risk_G2=rr[2],risk_G3=rr[3],risk_G4=rr[4])]
}
boot_file<-file.path(raw,'primary_bootstrap_replicates.csv')
reuse_bootstrap <- identical(Sys.getenv('REUSE_BOOTSTRAP', '0'), '1')
if(reuse_bootstrap && file.exists(boot_file)) boot<-fread(boot_file) else {
 boots<-list();for(bi in 1:500){boots[[bi]]<-boot_one(2026090500L+bi);if(bi%%50==0)cat('bootstrap=',bi,'\n',sep='')}
 boot<-rbindlist(boots,fill=TRUE);fwrite(boot,boot_file)
}
if(uniqueN(boot$replicate)<450)stop('Too few successful bootstrap samples')
bci<-boot[,.(OR_lo=quantile(OR,.025),OR_hi=quantile(OR,.975),RD_lo=quantile(RD,.025),RD_hi=quantile(RD,.975),p_value=2*pnorm(-abs(mean(log(OR))/sd(log(OR))))),by=comparison]
s$effects[,c('OR_lo','OR_hi','RD_lo','RD_hi','p_value'):=NULL];s$effects<-merge(s$effects,bci,by='comparison',sort=FALSE)
riskboot<-unique(boot[,.(replicate,risk_G1,risk_G2,risk_G3,risk_G4)]);rb<-melt(riskboot,id.vars='replicate',variable.name='group',value.name='risk')
rb[,group:=sub('risk_','',group)];rci<-rb[,.(lo=quantile(risk,.025),hi=quantile(risk,.975)),by=group]
s$risks[,c('lo','hi'):=NULL];s$risks<-merge(s$risks,rci,by='group',sort=FALSE)
fwrite(s$effects,file.path(main,'Table2_primary_effects.csv'))
fwrite(s$risks,file.path(main,'adjusted_risks.csv'))
# E-values are calculated on the weighted risk-ratio scale because the outcome
# is not rare. The confidence limit is obtained from the full-pipeline
# bootstrap distribution and the limit closest to the null is used.
rr_point<-s$risks[group=='G4',risk]/s$risks[group=='G1',risk]
rr_boot<-boot[comparison=='G4 vs G1',risk_G4/risk_G1]
rr_ci<-quantile(rr_boot,c(.025,.975),na.rm=TRUE)
evalue_rr<-function(rr) rr+sqrt(rr*(rr-1))
evalue_result<-data.table(
 comparison='G4 vs G1',scale='weighted risk ratio',RR=rr_point,
 RR_lo=rr_ci[[1]],RR_hi=rr_ci[[2]],
 E_value_point=evalue_rr(rr_point),
 E_value_limit=evalue_rr(rr_ci[[1]])
)
fwrite(evalue_result,file.path(supp,'evalue_diagnostics.csv'))
rawbal<-bal.tab(balance_formula,data=d,weights=o$w$weights,method='weighting',un=TRUE,binary='std',s.d.denom='pooled')
rawbt<-as.data.table(rawbal$Balance,keep.rownames='covariate')
fwrite(rawbt,file.path(supp,'balance_raw_overlap.csv'))
bal<-bal.tab(balance_formula,data=d,weights=d$overlap_weight,method='weighting',un=TRUE,binary='std',s.d.denom='pooled')
bt<-as.data.table(bal$Balance,keep.rownames='covariate');fwrite(bt,file.path(supp,'balance_all_standardized.csv'))
ess<-d[,.(n=.N,events=sum(incident_cam_positive),ESS=sum(overlap_weight)^2/sum(overlap_weight^2),
 weight_min=min(overlap_weight),weight_max=max(overlap_weight),dose_median=median(ivme_nci_48h)),by=ivme_group]
fwrite(ess,file.path(supp,'weight_diagnostics.csv'))
fwrite(d,file.path(raw,'analysis_data_with_weights.csv'))
diag<-data.table(n=nrow(d),events=sum(d$incident_cam_positive),max_SMD_un=max(bt$Max.Diff.Un,na.rm=TRUE),max_SMD_raw_overlap=max(rawbal$Balance$Max.Diff.Adj,na.rm=TRUE),max_SMD_adj=max(bt$Max.Diff.Adj,na.rm=TRUE),successful_bootstrap=uniqueN(boot$replicate),inference='percentile bootstrap of full estimation pipeline')
fwrite(diag,file.path(supp,'primary_diagnostics.csv'));print(diag);print(s$effects)
# Table 1: transparent preweighting descriptions, standardized balance reported separately.
fmt<-function(x) sprintf('%.1f [%.1f, %.1f]',median(x,na.rm=TRUE),quantile(x,.25,na.rm=TRUE),quantile(x,.75,na.rm=TRUE))
rows<-list()
for(v in c('age_at_admission','ivme_nci_48h',grep('^baseline_',names(d),value=TRUE))) {
 rows[[length(rows)+1]]<-d[,.(value=fmt(get(v))),by=ivme_group][,variable:=v]
}
for(v in c(factors,hist)) for(z in sort(unique(d[[v]]))) rows[[length(rows)+1]]<-d[,.(value=sprintf('%d (%.1f%%)',sum(get(v)==z),100*mean(get(v)==z))),by=ivme_group][,variable:=paste(v,z,sep=': ')]
tab<-dcast(rbindlist(rows),variable~ivme_group,value.var='value');fwrite(tab,file.path(main,'Table1_baseline.csv'))
# Main figures.
p<-ggplot(s$risks,aes(group,risk))+geom_errorbar(aes(ymin=lo,ymax=hi),width=.12,color='#24566E')+
 geom_point(size=3.5,color='#24566E')+scale_x_discrete(labels=labels)+scale_y_continuous(labels=scales::percent_format())+
 labs(title='Adjusted CAM-ICU positivity by early opioid dose',x='Recorded ICU IVME in the first 48 h (mg)',y='Overlap-weighted risk with 95% CI')
savefig(p,'Figure2_adjusted_risks')
bl<-merge(
 bt[,.(covariate,Unweighted=Max.Diff.Un,`After entropy calibration`=Max.Diff.Adj)],
 rawbt[,.(covariate,`Initial overlap weighting`=Max.Diff.Adj)],
 by='covariate',all=TRUE)
bl<-melt(bl,id.vars='covariate',variable.name='Analysis stage',value.name='SMD')
bl[,`Analysis stage`:=factor(`Analysis stage`,levels=c('Unweighted','Initial overlap weighting','After entropy calibration'))]
covariate_labels<-c(
 'age_at_admission'='Age',
 'I(age_at_admission^2)'='Age squared',
 'gender_M'='Male sex',
 'race_group_Asian'='Asian race',
 'race_group_Black'='Black race',
 'race_group_Hispanic'='Hispanic ethnicity',
 'race_group_Other_or_unknown'='Other or unknown race',
 'race_group_White'='White race',
 'admission_group_Elective'='Elective admission',
 'admission_group_Other'='Other admission type',
 'admission_group_Urgent_or_emergency'='Urgent or emergency admission',
 'procedure_class_Aortic_any'='Aortic procedure',
 'procedure_class_CABG_only'='Isolated CABG',
 'procedure_class_CABG_valve'='CABG and valve procedure',
 'procedure_class_Other'='Other procedure',
 'procedure_class_Valve_only'='Isolated valve procedure',
 'anchor_year_group_2008 - 2010'='Calendar period 2008-2010',
 'anchor_year_group_2011 - 2013'='Calendar period 2011-2013',
 'anchor_year_group_2014 - 2016'='Calendar period 2014-2016',
 'anchor_year_group_2017 - 2019'='Calendar period 2017-2019',
 'anchor_year_group_2020 - 2022'='Calendar period 2020-2022',
 'prior_admission_observed'='Prior hospitalization observed',
 'prior_chf'='Prior heart failure',
 'prior_mi'='Prior myocardial infarction',
 'prior_cerebrovascular_disease'='Prior cerebrovascular disease',
 'prior_chronic_pulmonary_disease'='Prior chronic pulmonary disease',
 'prior_diabetes'='Prior diabetes',
 'prior_chronic_kidney_disease'='Prior chronic kidney disease',
 'prior_liver_disease'='Prior liver disease')
bl[covariate %chin% names(covariate_labels),covariate:=unname(covariate_labels[covariate])]
balance_table<-dcast(copy(bl),covariate~`Analysis stage`,value.var='SMD')
balance_table_order<-c(
 'Age','Age squared','Male sex','Prior hospitalization observed','Prior heart failure',
 'Prior myocardial infarction','Prior cerebrovascular disease','Prior chronic pulmonary disease',
 'Prior diabetes','Prior chronic kidney disease','Prior liver disease','Asian race','Black race',
 'Hispanic ethnicity','Other or unknown race','White race','Elective admission','Other admission type',
 'Urgent or emergency admission','Aortic procedure','Isolated CABG','CABG and valve procedure',
 'Other procedure','Isolated valve procedure','Calendar period 2008-2010','Calendar period 2011-2013',
 'Calendar period 2014-2016','Calendar period 2017-2019','Calendar period 2020-2022')
balance_table[,covariate:=factor(covariate,levels=balance_table_order)]
setorder(balance_table,covariate)
balance_table[,covariate:=as.character(covariate)]
setnames(balance_table,
 c('Unweighted','Initial overlap weighting','After entropy calibration'),
 c('unweighted_max_pairwise_abs_smd','initial_overlap_max_pairwise_abs_smd','entropy_calibrated_max_pairwise_abs_smd'))
fwrite(balance_table,file.path(supp,'TableS6_covariate_balance.csv'))
p<-ggplot(bl,aes(abs(SMD),reorder(covariate,abs(SMD),max)))+geom_vline(xintercept=.1,linetype=2,color='grey50')+
 geom_line(aes(group=covariate),color='grey78',linewidth=.35)+
 geom_point(aes(color=`Analysis stage`,shape=`Analysis stage`),size=2.3)+
 scale_color_manual(values=c('Unweighted'='#B16B40','Initial overlap weighting'='#7B6FA8','After entropy calibration'='#24566E'))+
 scale_shape_manual(values=c('Unweighted'=16,'Initial overlap weighting'=17,'After entropy calibration'=15))+
 labs(title='Balance across all four dose groups',subtitle='Maximum pairwise absolute SMD; dashed line marks 0.10',x='Absolute SMD',y=NULL,color='Analysis stage',shape='Analysis stage')+
 theme(axis.text.y=element_text(size=8),legend.position='bottom')
savefig(p,'FigureS1_love_plot',supp,w=9.5,h=8.5)
# Weight magnitude is arbitrary after group-specific normalization. Rescale to a
# within-group mean of one so that the four distributions are comparable.
d[,plot_weight:=overlap_weight/mean(overlap_weight),by=ivme_group]
weight_plot_summary<-d[,.(n=.N,minimum=min(plot_weight),q1=quantile(plot_weight,.25),median=median(plot_weight),q3=quantile(plot_weight,.75),maximum=max(plot_weight)),by=ivme_group]
fwrite(weight_plot_summary,file.path(supp,'FigureS2_weight_distribution_summary.csv'))
p<-ggplot(d,aes(ivme_group,plot_weight,fill=ivme_group))+geom_violin(trim=TRUE,alpha=.7)+geom_boxplot(width=.15,outlier.shape=NA,fill='white')+
 scale_x_discrete(labels=labels)+guides(fill='none')+labs(title='Distribution of calibrated overlap weights',subtitle='Weights rescaled to a within-group mean of 1',x='0-48 h IVME (mg)',y='Relative calibrated overlap weight')
savefig(p,'FigureS2_weights',supp)
# Figure S3 is generated once in 07_figures_tables.R, which also exports its
# exact plotting data. Keeping a single implementation prevents silent overwrite.
# Assessment coverage and nonparametric worst/best bounds for the wider unassessed-inclusive cohort.
a_path<-file.path(raw,'assessment_eligible_covariates.csv')
if(!file.exists(a_path)) stop('Run 02_extract_covariates.R before 03_overlap_weighting_bootstrap.R')
a<-fread(a_path)
ass<-a[,.(eligible=.N,assessed=sum(valid_cam_assessment_n>0),positive=sum(incident_cam_positive,na.rm=TRUE),
 median_assessments=as.numeric(median(valid_cam_assessment_n[valid_cam_assessment_n>0])),
 q1_assessments=as.numeric(quantile(valid_cam_assessment_n[valid_cam_assessment_n>0],.25)),q3_assessments=as.numeric(quantile(valid_cam_assessment_n[valid_cam_assessment_n>0],.75))),by=ivme_group]
ass[,`:=`(assessed_pct=100*assessed/eligible,unassessed=eligible-assessed,lower_risk_bound=positive/eligible,upper_risk_bound=(positive+eligible-assessed)/eligible)]
fwrite(ass,file.path(supp,'TableS_assessment_and_bounds.csv'))
p<-ggplot(ass,aes(ivme_group,assessed_pct/100))+geom_col(fill='#24566E',width=.6)+geom_text(aes(label=sprintf('%.1f%%\n%d/%d',assessed_pct,assessed,eligible)),vjust=-.15,size=3.5)+
 scale_y_continuous(labels=scales::percent_format(),limits=c(0,1))+scale_x_discrete(labels=labels)+labs(title='Post-landmark CAM-ICU assessment availability',subtitle=sprintf('Among patients with classifiable opioid exposure (N = %s)',format(sum(ass$eligible),big.mark=',')),x='0-48 h IVME (mg)',y='Proportion assessed')
savefig(p,'Retired_assessment_coverage',internal_fig)
# Missing-outcome sensitivity for the prespecified G4-versus-G1 contrast.
# q1 and q4 denote the assumed delirium risks among patients without an
# interpretable post-landmark CAM-ICU assessment. These are crude sensitivity
# calculations, not adjusted causal estimates or corrections for selection.
g1<-ass[ivme_group=='G1'];g4<-ass[ivme_group=='G4']
risk_missing<-function(g,q)(g$positive+g$unassessed*q)/g$eligible
q1_observed<-g1$positive/g1$assessed;q4_observed<-g4$positive/g4$assessed
q4_example<-.30
q1_tipping<-(risk_missing(g4,q4_example)*g1$eligible-g1$positive)/g1$unassessed
tipping_scenarios<-data.table(
 scenario=c('All unassessed patients event-free','Group-specific observed-risk imputation',
  'All unassessed patients with events','Contrast-worst case for G4',
  'Contrast-best case for G4','Tipping example: G4 missing risk 30%'),
 q1=c(0,q1_observed,1,1,0,q1_tipping),q4=c(0,q4_observed,1,0,1,q4_example))
tipping_scenarios[,`:=`(risk_G1=risk_missing(g1,q1),risk_G4=risk_missing(g4,q4))]
tipping_scenarios[,risk_difference:=risk_G4-risk_G1]
fwrite(tipping_scenarios,file.path(supp,'TableS4_tipping_point_scenarios.csv'))
tipping_grid<-CJ(q1=seq(0,1,.01),q4=seq(0,1,.01))
tipping_grid[,`:=`(risk_G1=risk_missing(g1,q1),risk_G4=risk_missing(g4,q4))]
tipping_grid[,risk_difference_pp:=100*(risk_G4-risk_G1)]
fwrite(tipping_grid,file.path(supp,'FigureS6_tipping_point_data.csv'))
p<-ggplot(tipping_grid,aes(q4,q1,fill=risk_difference_pp))+
 geom_raster()+geom_contour(data=tipping_grid,aes(q4,q1,z=risk_difference_pp),inherit.aes=FALSE,
  breaks=0,color='black',linewidth=.9)+
 geom_abline(slope=1,intercept=0,linetype=3,color='grey25',linewidth=.7)+
 geom_point(data=data.table(q4=q4_observed,q1=q1_observed),aes(q4,q1),inherit.aes=FALSE,
  shape=21,size=3,stroke=.8,fill='white',color='black')+
 scale_fill_gradient2(low='#B16B40',mid='white',high='#24566E',midpoint=0,
  name='G4 - G1 risk\ndifference (pp)')+
 scale_x_continuous(labels=scales::percent_format(accuracy=1),breaks=seq(0,1,.25),expand=c(0,0))+
 scale_y_continuous(labels=scales::percent_format(accuracy=1),breaks=seq(0,1,.25),expand=c(0,0))+
 coord_equal()+labs(title='Tipping-point analysis for unassessed outcomes',
  subtitle='Black contour: no crude risk difference; dotted line: equal missing-outcome risks',
  x='Assumed delirium risk among unassessed G4 patients',
  y='Assumed delirium risk among unassessed G1 patients')+
 theme(legend.position='right')
savefig(p,'FigureS6_tipping_point',supp,w=8,h=6.5)

# Descriptive comparison of assessed and unassessed patients using only
# variables available for the full classifiable-exposure cohort.
a[,assessed_status:=valid_cam_assessment_n>0]
a[,race_group:=fifelse(grepl('^WHITE',race),'White',fifelse(grepl('^BLACK',race),'Black',
 fifelse(grepl('^HISPANIC',race),'Hispanic',fifelse(grepl('^ASIAN',race),'Asian','Other_or_unknown'))))]
a[,admission_group:=fifelse(admission_type %chin% c('ELECTIVE','SURGICAL SAME DAY ADMISSION'),'Elective',
 fifelse(grepl('URGENT|EW|EMER',admission_type),'Urgent_or_emergency','Other'))]
a[,procedure_class:=fcase(grepl('aortic',procedure_groups),'Aortic_any',
 grepl('CABG',procedure_groups)&grepl('valve',procedure_groups),'CABG_valve',
 procedure_groups=='CABG','CABG_only',procedure_groups=='valve','Valve_only',default='Other')]
z<-a$assessed_status;n_assessed<-sum(z);n_unassessed<-sum(!z)
comparison_rows<-list()
add_cont_row<-function(label,x){
 m1<-mean(x[z],na.rm=TRUE);m0<-mean(x[!z],na.rm=TRUE);s1<-sd(x[z],na.rm=TRUE);s0<-sd(x[!z],na.rm=TRUE)
 den<-sqrt((s1^2+s0^2)/2);smd<-ifelse(den>0,abs(m1-m0)/den,0)
 comparison_rows[[length(comparison_rows)+1]]<<-data.table(characteristic=label,
  assessed=sprintf('%.1f (%.1f)',m1,s1),unassessed=sprintf('%.1f (%.1f)',m0,s0),absolute_SMD=smd)
}
add_bin_row<-function(label,x){
 x<-as.logical(x);p1<-mean(x[z],na.rm=TRUE);p0<-mean(x[!z],na.rm=TRUE)
 den<-sqrt((p1*(1-p1)+p0*(1-p0))/2);smd<-ifelse(den>0,abs(p1-p0)/den,0)
 comparison_rows[[length(comparison_rows)+1]]<<-data.table(characteristic=label,
  assessed=sprintf('%d/%d (%.1f%%)',sum(x[z],na.rm=TRUE),n_assessed,100*p1),
  unassessed=sprintf('%d/%d (%.1f%%)',sum(x[!z],na.rm=TRUE),n_unassessed,100*p0),absolute_SMD=smd)
}
add_cont_row('Age, years',a$age_at_admission)
add_bin_row('Male sex',a$gender=='M')
for(v in c('White','Black','Hispanic','Asian','Other_or_unknown'))
 add_bin_row(c(White='White race',Black='Black race',Hispanic='Hispanic ethnicity',Asian='Asian race',Other_or_unknown='Other or unknown race')[v],a$race_group==v)
for(v in c('Elective','Urgent_or_emergency','Other'))
 add_bin_row(c(Elective='Elective admission',Urgent_or_emergency='Urgent or emergency admission',Other='Other admission type')[v],a$admission_group==v)
for(v in c('CABG_only','Valve_only','CABG_valve','Aortic_any','Other'))
 add_bin_row(c(CABG_only='Isolated CABG',Valve_only='Isolated valve procedure',CABG_valve='CABG and valve procedure',Aortic_any='Aortic procedure',Other='Other cardiac procedure')[v],a$procedure_class==v)
for(v in sort(unique(a$anchor_year_group)))
 add_bin_row(paste('Calendar period',v),a$anchor_year_group==v)
add_bin_row('Prior hospitalization observed',a$prior_admission_observed==1)
prior_labels<-c(
 prior_chf='Prior heart failure',prior_mi='Prior myocardial infarction',
 prior_cerebrovascular_disease='Prior cerebrovascular disease',
 prior_chronic_pulmonary_disease='Prior chronic pulmonary disease',
 prior_diabetes='Prior diabetes',prior_chronic_kidney_disease='Prior chronic kidney disease',
 prior_liver_disease='Prior liver disease')
for(v in names(prior_labels)) add_bin_row(prior_labels[[v]],a[[v]]==1)
lab_labels<-c(
 baseline_creatinine='Pre-ICU creatinine',baseline_hemoglobin='Pre-ICU hemoglobin',
 baseline_sodium='Pre-ICU sodium',baseline_wbc='Pre-ICU white blood cell count',
 baseline_glucose='Pre-ICU glucose',baseline_lactate='Pre-ICU lactate',
 baseline_albumin='Pre-ICU albumin')
for(v in names(lab_labels)) add_bin_row(paste(lab_labels[[v]],'available'),!is.na(a[[v]]))
for(v in c('baseline_glucose','baseline_lactate')) add_cont_row(lab_labels[[v]],a[[v]])
for(v in paste0('G',1:4))add_bin_row(paste('IVME group',v),a$ivme_group==v)
assessed_comparison<-rbindlist(comparison_rows)
fwrite(assessed_comparison,file.path(supp,'TableS5_assessed_vs_unassessed.csv'))
# Composition is descriptive, not an estimated drug-specific treatment effect.
d[,drug_count:=(fentanyl_mcg_48h>0)+(hydromorphone_mg_48h>0)+(morphine_mg_48h>0)]
d[,drug_pattern:=fcase(drug_count==0,'No specified opioid',drug_count>1,'Mixed opioids',fentanyl_mcg_48h>0,'Fentanyl only',hydromorphone_mg_48h>0,'Hydromorphone only',default='Morphine only')]
comp<-d[,.(n=.N,events=sum(incident_cam_positive),median_IVME=median(ivme_nci_48h)),by=.(ivme_group,drug_pattern)]
fwrite(comp,file.path(supp,'TableS_drug_composition.csv'))
p<-ggplot(comp,aes(ivme_group,n,fill=drug_pattern))+geom_col(position='fill',width=.65)+scale_y_continuous(labels=scales::percent_format())+scale_x_discrete(labels=labels)+
 scale_fill_brewer(palette='Set2')+labs(title='Opioid composition within dose categories',x='0-48 h IVME (mg)',y='Proportion',fill=NULL)
savefig(p,'FigureS4_drug_composition',supp)
# Sensitivity 1: documented early CAM-ICU negative, refit weights in this restricted cohort.
neg<-droplevels(d[early_negative_n>0]);on<-fitow(neg);sn<-summarize_fit(on,neg,'Early CAM-ICU negative ATO')
bn<-bal.tab(psf,data=neg,weights=on$calibrated_weights,method='weighting',un=TRUE,binary='std');fwrite(as.data.table(bn$Balance,keep.rownames='covariate'),file.path(supp,'early_negative_balance.csv'))
# Sensitivity 2: expanded baseline laboratory adjustment; complete case for lactate and glucose.
# This changes the included population and is explicitly exploratory, not a missing-data correction.
lab<-droplevels(d[!is.na(baseline_lactate)&baseline_lactate>=0&!is.na(baseline_glucose)&baseline_glucose>0])
fl<-update(psf,.~.+log1p(baseline_lactate)+log(baseline_glucose));ol<-fitow(lab,fl);sl<-summarize_fit(ol,lab,'Pre-ICU lactate/glucose complete-case ATO')
fwrite(as.data.table(bal.tab(fl,data=lab,weights=ol$calibrated_weights,method='weighting',un=TRUE,binary='std')$Balance,keep.rownames='covariate'),file.path(supp,'laboratory_balance.csv'))
# Full-pipeline bootstrap for the two weighted sensitivity analyses.
boot_sensitivity<-function(dat,form,path,seedbase,R=300){
 if(reuse_bootstrap && file.exists(path)) b<-fread(path) else {
  out<-list();for(k in 1:R){set.seed(seedbase+k);bd<-dat[sample.int(nrow(dat),replace=TRUE)];oo<-try(fitow(bd,form),silent=TRUE);if(inherits(oo,'try-error'))next
   bd[,bw:=oo$calibrated_weights];rr<-bd[,.(risk=weighted.mean(incident_cam_positive,bw)),by=ivme_group][order(ivme_group)]$risk;if(length(rr)!=4||any(rr<=0|rr>=1))next
   pp<-list(c(2,1),c(3,1),c(4,1),c(4,2));out[[length(out)+1]]<-rbindlist(lapply(pp,function(ij)data.table(replicate=k,comparison=paste0('G',ij[1],' vs G',ij[2]),OR=(rr[ij[1]]/(1-rr[ij[1]]))/(rr[ij[2]]/(1-rr[ij[2]])),RD=rr[ij[1]]-rr[ij[2]])))
  };b<-rbindlist(out);fwrite(b,path)
 }
 if(uniqueN(b$replicate)<.9*R)stop('Sensitivity bootstrap failure rate too high')
 b[,.(OR_lo=quantile(OR,.025),OR_hi=quantile(OR,.975),RD_lo=quantile(RD,.025),RD_hi=quantile(RD,.975),p_value=2*pnorm(-abs(mean(log(OR))/sd(log(OR))))),by=comparison]
}
replace_ci<-function(point,ci){point[,c('OR_lo','OR_hi','RD_lo','RD_hi','p_value'):=NULL];merge(point,ci,by='comparison',sort=FALSE)}
sn$effects<-replace_ci(sn$effects,boot_sensitivity(neg,psf,file.path(raw,'early_negative_bootstrap.csv'),2026091000L))
sl$effects<-replace_ci(sl$effects,boot_sensitivity(lab,fl,file.path(raw,'laboratory_bootstrap.csv'),2026092000L))
# Sensitivity 3: remove the upper 5% of the positive-dose distribution. This is
# an influence analysis for the sparse high-dose tail, not a threshold search.
tail_cut<-as.numeric(quantile(d[ivme_nci_48h>0,ivme_nci_48h],.95,type=7))
tail<-droplevels(d[ivme_nci_48h<=tail_cut]);ot<-fitow(tail);st<-summarize_fit(ot,tail,'Positive-dose upper-tail trimmed ATO')
tail_boot_path<-file.path(raw,'upper_tail_trim_bootstrap.csv')
st$effects<-replace_ci(st$effects,boot_sensitivity(tail,psf,tail_boot_path,2026093000L))
tail_bal<-bal.tab(psf,data=tail,weights=ot$calibrated_weights,method='weighting',un=TRUE,binary='std',s.d.denom='pooled')
tail_bal_dt<-as.data.table(tail_bal$Balance,keep.rownames='covariate')
tail[,tail_weight:=ot$calibrated_weights]
tail_groups<-tail[,.(n=.N,events=sum(incident_cam_positive),crude_risk=mean(incident_cam_positive),
 weighted_risk=weighted.mean(incident_cam_positive,tail_weight),ESS=sum(tail_weight)^2/sum(tail_weight^2)),by=ivme_group]
tail_diag<-data.table(
 analysis='Positive-dose upper-tail trimmed ATO',positive_dose_percentile=.95,cutoff_mg=tail_cut,
 excluded_n=nrow(d)-nrow(tail),included_n=nrow(tail),included_events=sum(tail$incident_cam_positive),
 max_pairwise_abs_SMD=max(tail_bal_dt$Max.Diff.Adj,na.rm=TRUE),
 successful_bootstrap=uniqueN(fread(tail_boot_path)$replicate))
fwrite(tail_groups,file.path(supp,'TableS_upper_tail_trim_groups.csv'))
fwrite(st$effects,file.path(supp,'TableS_upper_tail_trim_effects.csv'))
fwrite(tail_bal_dt,file.path(supp,'upper_tail_trim_balance.csv'))
fwrite(tail_diag,file.path(supp,'upper_tail_trim_diagnostics.csv'))
fwrite(sn$effects,file.path(supp,'TableS_early_negative_effects.csv'));fwrite(sl$effects,file.path(supp,'TableS_laboratory_effects.csv'))
fwrite(rbindlist(list(s$effects,sn$effects,sl$effects,st$effects),use.names=TRUE,fill=TRUE),file.path(supp,'all_OW_effects.csv'))
saveRDS(list(primary=o,early_negative=on,laboratory=ol,upper_tail_trimmed=ot,rhs=rhs,psf=psf),file.path(raw,'models.rds'))
capture.output(sessionInfo(),file=file.path(env_dir,'R_sessionInfo.txt'))
cat('PRIMARY AND SENSITIVITY MODELS COMPLETE\n')
