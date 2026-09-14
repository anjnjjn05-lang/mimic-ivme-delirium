suppressPackageStartupMessages({library(data.table);library(ggplot2)})
try(Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8'),silent=TRUE)
root<-Sys.getenv('ANALYSIS_ROOT','..');raw<-Sys.getenv('ANALYSIS_RAW',file.path(root,'outputs','restricted_patient_level'));main<-Sys.getenv('ANALYSIS_MAIN',file.path(root,'outputs','aggregate_review','main'));supp<-Sys.getenv('ANALYSIS_SUPP',file.path(root,'outputs','aggregate_review','supplementary'))
theme_set(theme_minimal(base_size=12)+theme(panel.grid=element_blank(),plot.title=element_text(face='bold'),legend.position='bottom'))
savefig<-function(p,name,folder=main,w=8,h=5){ggsave(file.path(folder,paste0(name,'.png')),p,width=w,height=h,dpi=300,bg='white');ggsave(file.path(folder,paste0(name,'.pdf')),p,width=w,height=h,device=cairo_pdf)}
# Cohort flow. All values are read from files produced by 01_build_analysis_dataset.R; no count
# is duplicated manually in the plotting code.
identification<-fread(file.path(supp,'identification_counts.csv'))
landmark_flow<-fread(file.path(supp,'cohort_flow.csv'))
stopifnot(nrow(identification)==1,nrow(landmark_flow)==4)
flow<-data.table(stage=c('Literature-mapped procedure-ICU pairs','Unique patients after first eligible ICU selection','Adult patients alive, hospitalized, and in ICU at 48 h\nwithout recorded prior dementia','No recorded positive CAM-ICU during 0-48 h','At least one interpretable CAM-ICU after 48 h','Final analysis after dose-record quality check'),n=c(identification$procedure_icu_pairs,identification$unique_patients,landmark_flow$n))
flow[,y:=rev(seq_len(.N))];fwrite(flow[,.(stage,n)],file.path(supp,'cohort_flow_full.csv'))
p<-ggplot(flow)+geom_segment(data=flow[-.N],aes(x=0,y=y-.32,xend=0,yend=y-.68),arrow=arrow(length=unit(.12,'inches')),linewidth=.45,color='grey40')+
 geom_label(aes(0,y,label=sprintf('%s\nN = %s',stage,format(n,big.mark=','))),size=3.5,linewidth=.35,fill='#F5F8FA',color='#163B50',lineheight=1.05)+
 coord_cartesian(xlim=c(-1,1),clip='off')+labs(title='Cohort construction',x=NULL,y=NULL)+theme(axis.text=element_blank(),axis.ticks=element_blank(),plot.margin=margin(10,70,10,70))
savefig(p,'Figure1_cohort_flow',main,w=8,h=7.2)
# Generalized propensity scores from the pre-calibration multinomial model.
o<-readRDS(file.path(raw,'models.rds'))$primary
gps<-as.data.table(o$w$obj$fitted.values);setnames(gps,paste0('Pr_',levels(o$w$treat)));gps[,observed:=o$w$treat]
gl<-melt(gps,id.vars='observed',variable.name='score',value.name='probability')
fwrite(gl,file.path(raw,'generalized_propensity_scores.csv'))
p<-ggplot(gl,aes(probability,color=observed))+geom_density(linewidth=.65,adjust=1.1)+facet_wrap(~score,ncol=2)+
 scale_color_brewer(palette='Dark2')+labs(title='Generalized propensity score distributions before calibration',x='Estimated probability',y='Density',color='Observed group')
savefig(p,'FigureS3_propensity_overlap',supp,w=8,h=6)
# Dose-extraction audit derived only from the current reproducible cohort.
d<-fread(file.path(raw,'final_cohort.csv'));pos<-d[ivme_nci_48h>0]
qa<-data.table(metric=c('Positive-dose P25','Positive-dose median','Positive-dose P75','Maximum dose','Patients with unsupported event excluded from final model'),value=c(quantile(pos$ivme_nci_48h,.25),median(pos$ivme_nci_48h),quantile(pos$ivme_nci_48h,.75),max(pos$ivme_nci_48h),1))
fwrite(qa,file.path(supp,'dose_extraction_audit.csv'));print(qa)
