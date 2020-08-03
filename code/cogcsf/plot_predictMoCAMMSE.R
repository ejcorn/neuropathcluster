# relationship between path scores and cognition
rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'cogcluster/predictcog/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')
source('code/misc/trainfxns.R')

###########################################
### load data from spearman correlation ###
###########################################

p.meth <- 'cor'
load(file = paste0(savedir,'PredictCogFromPath_Model_',p.meth,'.RData'))
x.plot <- 'AllPathology'
coef.matrices <- lapply(y.names, function(y.name) region.by.item.matrix(m.all[[x.plot]][[y.name]]$r * (m.all[[x.plot]][[y.name]]$p<0.05))) # get fdr corrected correlations
names(coef.matrices) <- y.names
clim <- c(min(unlist(coef.matrices),na.rm = T),max(unlist(coef.matrices),na.rm = T))
p.list <- lapply(y.names, function(y.name) imagesc(coef.matrices[[y.name]],cmap='redblue_asymmetric',clim =clim ) + theme(legend.key.height=unit(0.3,'cm'),legend.key.width = unit(0.1,'cm'),legend.position = 'right')+
                   ggtitle(paste0(y.name)) + theme(plot.title = element_text(hjust=0.5,size=8),axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),plot.margin = unit(c(0,0,0,0),'cm')))
p.cor <- plot_grid(plotlist = p.list,nrow=1,align='hv')
ggsave(plot = p.cor,filename = paste0(savedir,'PredictCogFeatureWeights_',p.meth,'.pdf'),
       width = 20, height = 3, units = "cm")
save(coef.matrices,file = paste0(params$sourcedata.dir,'Fig5a_SourceData.RData'))
#######################################
### load data from model prediction ###
#######################################

x.labels <-c(AllPathology='All Pathology',ADNC='ADNC',LBD='LBD Type',Cluster='Cluster',EFA='EFA',PathologyType='Path Type Average',RegionalPathology='Regional Pathology Average')
for(p.meth in c('rf','glmnet')){
  load(file = paste0(savedir,'PredictCogFromPath_Model_',p.meth,'.RData'))
  comps <- combn(x.names,m = 2,simplify = FALSE)
  
  rsq <- lapply(y.names, function(y.name) 
    sapply(x.names, function(x.name) extract.bestTune.metric(m.all[[x.name]][[y.name]],metric = 'Rsquared')))
  names(rsq) <- y.names
  rsq.mean <- lapply(rsq, colMeans)
  rsq.plt <- lapply(rsq, function(X) collapse.columns(X,cnames = x.names))
  rsq.ord <- lapply(rsq.mean, function(X) names(sort(X,decreasing = T))) # order for plotting. sort by R^2
  ymax <- max(unlist(rsq)) # scale to overall max r squared observed
  p.list <- lapply(y.names, function(y.name) ggplot(rsq.plt[[y.name]],aes(x=names,y=values,fill=names)) + geom_boxplot() + theme_bw()+
                     scale_fill_brewer(palette='Blues',limits=names(sort(rsq.mean[[y.name]]))) + 
                     scale_y_continuous(limits=c(0,ymax)) + scale_x_discrete(limits=rsq.ord[[y.name]],labels=x.labels[rsq.ord[[y.name]]])+
                     xlab('')+ylab(expression(R^2)) + ggtitle(y.name) +
                     theme(text=element_text(size=8),plot.title = element_text(size=8,hjust=0.5),legend.position = 'none')+
                     theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) #+ stat_compare_means(size=1.75,comparisons = comps,position=5,method = "wilcox.test",label = 'p.signif')
                   )
  p.rsq <- plot_grid(plotlist = p.list,nrow=1,align='hv')
  ggsave(plot = p.rsq,filename = paste0(savedir,'PredictCogPerformance_',p.meth,'.pdf'),
         width = 18, height = 7, units = "cm")
  if(p.meth == 'rf'){save(rsq.plt,file = paste0(params$sourcedata.dir,'Fig5c_SourceData.RData'))}
  # when using all pathological features, plot coefficients to see regional structure-function mappings
  x.plot <- 'AllPathology'
  coef.matrices <- lapply(y.names, function(y.name) extract.coefs(m.all[[x.plot]][[y.name]],p.meth)) # get coefficients
  names(coef.matrices) <- y.names
  cmap <- ifelse(p.meth %in% c('lm','glmnet'),yes='redblue_asymmetric',no='OrRd')
  p.list <- lapply(y.names, function(y.name) imagesc(region.by.item.matrix(coef.matrices[[y.name]]),cmap=cmap) + theme(legend.key.height=unit(0.3,'cm'),legend.key.width = unit(0.1,'cm'),legend.position = 'right')+
                     ggtitle(paste0(y.name,': ',expression(R^2),' = ',signif(mean(rsq[[y.name]][,x.plot]),2))) + 
                     theme(plot.title = element_text(hjust=0.5,size=8),axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),plot.margin = unit(c(0,0,0,0),'cm')))
  p.fw <- plot_grid(plotlist = p.list,nrow=1,align='hv')
  ggsave(plot = p.fw,filename = paste0(savedir,'PredictCogFeatureWeights_',p.meth,'.pdf'),
         width = 20, height = 3, units = "cm")
  if(p.meth == 'rf'){save(coef.matrices,file = paste0(params$sourcedata.dir,'Fig5b_SourceData.RData'))}
  
  # p.all <- plot_grid(plotlist=list(p.fw,p.rsq),align='hv',nrow=2,rel_heights = c(1,1.2))
  # ggsave(plot = p.all,filename = paste0(savedir,'PredictCogAll',p.meth,'.pdf'),
  #        width = 18, height = 12, units = "cm")
} 