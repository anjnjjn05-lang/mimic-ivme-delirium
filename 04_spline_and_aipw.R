suppressPackageStartupMessages({library(data.table);library(ggplot2);library(Hmisc);library(glmnet);library(sandwich)})
options(warn=1)
try(Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8'),silent=TRUE)
root<-Sys.getenv('ANALYSIS_ROOT','..');raw<-Sys.getenv('ANALYSIS_RAW',file.path(root,'outputs','restricted_patient_level'));main<-Sys.getenv('ANALYSIS_MAIN',file.path(root,'outputs','aggregate_review','main'));supp<-Sys.getenv('ANALYSIS_SUPP',file.path(root,'outputs','aggregate_review','supplementary'))
d<-fread(file.path(raw,'analysis_data_with_weights.csv'));saved<-readRDS(file.path(raw,'models.rds'));rhs<-saved$rhs
d[,ivme_group:=factor(ivme_group,levels=paste0('G',1:4))]
for(v in c('gender','race_group','admission_group','procedure_class','anchor_year_group'))d[,(v):=factor(get(v))]
theme_set(theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face='bold'),legend.position='bottom'))
savefig<-function(p,name,folder=main,w=8,h=5){ggsave(file.path(folder,paste0(name,'.png')),p,width=w,height=h,dpi=300,bg='white');ggsave(file.path(folder,paste0(name,'.pdf')),p,width=w,height=h,device=cairo_pdf)}
# Raw dose restricted cubic spline: three knots; two dose parameters.
pos<-droplevels(d[ivme_nci_48h>0]);knots<-as.numeric(quantile(pos$ivme_nci_48h,c(.1,.5,.9)))
basis<-function(z) Hmisc::rcspline.eval(z,knots=knots,inclx=TRUE)
bs<-basis(pos$ivme_nci_48h);pos[,`:=`(dose_linear=bs[,1],dose_nonlinear=bs[,2])]
fm<-as.formula(paste('incident_cam_positive~dose_linear+dose_nonlinear+',rhs));fl<-update(fm,.~.-dose_nonlinear);fn<-update(fm,.~.-dose_linear-dose_nonlinear)
fit<-glm(fm,pos,family=binomial());lin<-glm(fl,pos,family=binomial());nul<-glm(fn,pos,family=binomial())
stopifnot(fit$converged,!anyNA(coef(fit)))
V<-vcovHC(fit,type='HC0');br<-quantile(pos$ivme_nci_48h,c(.01,.99));grid<-exp(seq(log(br[1]),log(br[2]),length.out=160))
curve<-rbindlist(lapply(grid,function(z){
 nd<-copy(pos);b<-basis(rep(z,nrow(nd)));nd[,`:=`(dose_linear=b[,1],dose_nonlinear=b[,2])]
 X<-model.matrix(delete.response(terms(fit)),nd);p<-plogis(X%*%coef(fit));risk<-mean(p)
 g<-colMeans(X*as.numeric(p*(1-p)));se<-sqrt(as.numeric(t(g)%*%V%*%g));sl<-se/(risk*(1-risk))
 data.table(dose=z,risk=risk,lo=plogis(qlogis(risk)-1.96*sl),hi=plogis(qlogis(risk)+1.96*sl))
}))
test<-data.table(n=nrow(pos),events=sum(pos$incident_cam_positive),knot10=knots[1],knot50=knots[2],knot90=knots[3],
 overall_p=anova(nul,fit,test='Chisq')$`Pr(>Chi)`[2],nonlinear_p=anova(lin,fit,test='Chisq')$`Pr(>Chi)`[2],display_low=br[1],display_high=br[2])
fwrite(test,file.path(supp,'spline_diagnostics.csv'));fwrite(curve,file.path(main,'dose_response_curve.csv'))
p<-ggplot(curve,aes(dose,risk))+geom_ribbon(aes(ymin=lo,ymax=hi),fill='#88BBCD',alpha=.5)+geom_vline(xintercept=knots,linetype=3,color='grey45',linewidth=.45)+geom_line(color='#24566E',linewidth=1)+
 geom_rug(data=pos[ivme_nci_48h>=br[1]&ivme_nci_48h<=br[2]],aes(x=ivme_nci_48h),inherit.aes=FALSE,alpha=.1,sides='b')+
 scale_x_log10()+scale_y_continuous(labels=scales::percent_format())+labs(title='Continuous dose association among positive-dose patients',
 subtitle=sprintf('Raw-dose spline; nonlinear P = %.4f; dotted lines mark knots at %s mg',test$nonlinear_p,paste(format(knots,trim=TRUE),collapse=', ')),x='0-48 h IVME (mg, log axis)',y='Covariate-standardized risk with 95% CI')
savefig(p,'Figure3_dose_response')
saveRDS(list(model=fit,vcov=V,knots=knots),file.path(raw,'spline_model.rds'))
# Cross-fitted AIPW, G4 vs G1. Binary ATE, NOT four-group ATO.
q<-droplevels(d[ivme_group %in% c('G1','G4')]);a<-as.integer(q$ivme_group=='G4');y<-q$incident_cam_positive
# Use fixed age knots to avoid full-cohort learned preprocessing in outer folds.
xform<-as.formula(paste('~',gsub('splines::ns(age_at_admission,df=3)','splines::ns(age_at_admission,knots=c(50,70),Boundary.knots=c(18,100))',rhs,fixed=TRUE)))
X<-model.matrix(xform,q)[,-1,drop=FALSE]
fitprob<-function(xx,yy,xt){
 inner<-integer(length(yy));for(v in unique(yy)){ii<-which(yy==v);inner[ii]<-sample(rep(1:5,length.out=length(ii)))}
 cv<-cv.glmnet(xx,yy,family='binomial',alpha=1,foldid=inner,type.measure='deviance')
 as.numeric(predict(cv,xt,s='lambda.1se',type='response'))
}
runs<-list();preds<-list()
for(seed in c(20260905L,20260906L,20260907L)){
 set.seed(seed);fold<-integer(nrow(q));st<-interaction(a,y)
 for(z in levels(st)){ii<-which(st==z);fold[ii]<-sample(rep(1:5,length.out=length(ii)))}
 e<-m0<-m1<-rep(NA_real_,nrow(q))
 for(k in 1:5){te<-which(fold==k);tr<-which(fold!=k);t1<-tr[a[tr]==1];t0<-tr[a[tr]==0]
  e[te]<-fitprob(X[tr,,drop=FALSE],a[tr],X[te,,drop=FALSE]);m1[te]<-fitprob(X[t1,,drop=FALSE],y[t1],X[te,,drop=FALSE]);m0[te]<-fitprob(X[t0,,drop=FALSE],y[t0],X[te,,drop=FALSE])}
 eu<-pmin(.975,pmax(.025,e));phi1<-m1+a/eu*(y-m1);phi0<-m0+(1-a)/(1-eu)*(y-m0);p1<-mean(phi1);p0<-mean(phi0)
 stopifnot(p1>0,p1<1,p0>0,p0<1)
 i1<-phi1-p1;i0<-phi0-p0;rd<-p1-p0;sr<-sd(i1-i0)/sqrt(nrow(q));bo<-qlogis(p1)-qlogis(p0);so<-sd(i1/(p1*(1-p1))-i0/(p0*(1-p0)))/sqrt(nrow(q))
 runs[[length(runs)+1]]<-data.table(seed=seed,n=nrow(q),events=sum(y),risk_G1=p0,risk_G4=p1,RD=rd,RD_lo=rd-1.96*sr,RD_hi=rd+1.96*sr,OR=exp(bo),OR_lo=exp(bo-1.96*so),OR_hi=exp(bo+1.96*so),clipped=sum(eu!=e),ps_min=min(e),ps_max=max(e))
 preds[[length(preds)+1]]<-data.table(seed=seed,stay_id=q$stay_id,fold=fold,ps=e,m0=m0,m1=m1)
}
ml<-rbindlist(runs);fwrite(ml,file.path(supp,'TableS_ML_AIPW_ATE.csv'));fwrite(rbindlist(preds),file.path(raw,'ML_crossfit_predictions.csv'))
# Separate facets prevent presenting unlike target populations as one pooled estimate.
ef<-fread(file.path(supp,'all_OW_effects.csv'))[comparison=='G4 vs G1']
ef<-ef[,.(analysis,OR,OR_lo,OR_hi,target='Four-group overlap populations')]
ef<-rbind(ef,data.table(analysis=paste0('Cross-fitted AIPW seed ',ml$seed),OR=ml$OR,OR_lo=ml$OR_lo,OR_hi=ml$OR_hi,target='G1/G4 binary ATE population'))
display_order<-c('Primary four-group ATO','Early CAM-ICU negative ATO','Pre-ICU lactate/glucose complete-case ATO','Positive-dose upper-tail trimmed ATO',paste0('Cross-fitted AIPW seed ',sort(ml$seed)))
ef[,analysis:=factor(analysis,levels=rev(display_order))]
p<-ggplot(ef,aes(OR,analysis))+geom_vline(xintercept=1,linetype=2,color='grey50')+geom_errorbar(aes(xmin=OR_lo,xmax=OR_hi),width=.15,orientation='y')+
 geom_point(color='#24566E',size=2.5)+scale_x_log10()+facet_grid(target~.,scales='free_y',space='free_y')+labs(title='G4 versus G1 estimates across analysis specifications',x='Odds ratio with 95% CI',y=NULL)+theme(axis.text.y=element_text(size=9),strip.text.y=element_text(size=8))
savefig(p,'FigureS5_sensitivity_estimates',supp,w=10,h=5.5)
print(test);print(ml);cat('DOSE AND ML COMPLETE\n')
