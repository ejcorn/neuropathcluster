##############
### Colors ###
##############

getClusterColors <- function(k){
  pal <- colorRampPalette(brewer.pal(name = 'Dark2',n=8))
  return(pal(8)[(1:k)+2])
}

################
### p-values ###
################

p.signif <- function(p){
  # take p values and make asterisks
  #   ns: p > 0.05
  # *: p <= 0.05
  # **: p <= 0.01
  # ***: p <= 0.001
  # ****: p <= 0.000001
  p.new <- rep('',length(p))
  p.new[p > 0.05] <- 'ns'
  p.new[p < 0.05 & p > 0.01] <- '*'
  p.new[p <= 0.01 & p > 0.001] <- '**'
  p.new[p <= 0.001 & p > 0.000001] <- '***'
  p.new[p <= 0.000001] <- '****'
  names(p.new) <- names(p)
  return(p.new)
  
}

p.signif.matrix <- function(p){
  # take matrix of p values and make asterisks
  #   ns: p > 0.05
  # *: p <= 0.05
  # **: p <= 0.01
  # ***: p <= 0.001
  # ****: p <= 0.000001
  p.new <- matrix(data = '',nrow = nrow(p),ncol=ncol(p))
  p.new[p > 0.05] <- 'ns'
  p.new[p < 0.05 & p > 0.01] <- '*'
  p.new[p <= 0.01 & p > 0.001] <- '**'
  p.new[p <= 0.001 & p > 0.000001] <- '***'
  p.new[p <= 0.000001] <- '****'
  return(p.new)
  
}

matrix.fdr.correct <- function(pvals){
 pvals.mat <- matrix(p.adjust(as.vector(as.matrix(pvals)), method='fdr'),ncol=ncol(pvals))
 colnames(pvals.mat) <- colnames(pvals)
 return(pvals.mat) 
}

list.vec.fdr.correct <- function(X){
  # for a list where each element is a vector of p-values
  # FDR correct over all p-values and return to list
  # correct
  p.fdr <- p.adjust(unlist(X),method='fdr')
  # get lengths --> cumulative indices
  X.ind <- cumsum(sapply(X,length))
  # initialize
  X.fdr <- list()
  X.names <- names(X)
  # start off list with first portion
  X.fdr[[X.names[1]]] <- p.fdr[1:X.ind[1]]
  # fill in list elements 2 to n
  for(i in 2:length(X.ind)){
    X.fdr[[X.names[i]]] <- p.fdr[(1+X.ind[i-1]):(X.ind[i])]
  }
  return(X.fdr)
}

list.fdr.correct <- function(X){
  # unlist a list, fdr correct over all values
  # relist the list in the same structure and return
  return(relist(flesh=p.adjust(unlist(X),method='fdr'),skeleton=X))
}

#############
### plots ###
#############

mult.comp.ggpubr <- function(p,i=2,method='fdr'){
  
  # build plot
  p <- ggplot_build(p)
  p.orig <- p$data[[i]]$annotation
  # it lists any v. low p-val as '< 2.2e-16'
  # to be conservative just change to 2.2e-16
  levels(p.orig)[levels(p.orig) == '< 2.2e-16'] <- 2.2e-16
  # in ggplot object p-value is listed 3 times for each line segment of the bars
  # choose every 3rd up to end
  sel <- seq(3,length(p.orig),by=3)
  # extract and adjust p values
  p.adj <- p.adjust(as.numeric(as.character(p.orig)[sel]),method=method)
  # replace p-values with adjusted ones
  p$data[[i]]$annotation <- as.factor(p.signif(rep(p.adj,each=3)))
  return(p)
}

fix.yscale <- function(p,ymax,n.breaks,i=2){
  # this function will take a ggpubr plot (after ggplot_build)
  # and set the y axis position for the annotations to accommodate a 
  # set y range
  # using this because there are 2 massive outliers in CSF TTau
  # I don't want to remove them from analysis but don't want to distort axes

  # compute difference between actual ymax and desired ymax for annotations
  #offset <- max(p$data[[i]]$yend) - ymax
  # compute scaling between actual y range and desired y range (assumes ymin=0)
  scl <- diff(p$layout$panel_params[[1]]$y.range) / ymax  
  scl <- scl*1.05 # makes the labels a little bit smaller just because
  # reduce y values by that amount
  p$data[[i]]$y <- p$data[[i]]$y/scl
  p$data[[i]]$yend <- p$data[[i]]$yend/scl

  # adjust plot y-axis for all facets
  for(f in 1:length(p$layout$panel_params)){
   p$layout$panel_params[[f]]$y.range <- c(0-0.05*ymax,1.05*ymax) 
   p$layout$panel_params[[f]]$y.labels <- as.character(seq(0,ymax,by=ymax/n.breaks))
   p$layout$panel_params[[f]]$y.major_source <- seq(0,ymax,by=ymax/n.breaks)
   p$layout$panel_params[[f]]$y.minor_source <- seq(0,ymax,by=ymax/(2*n.breaks))
  }
  return(p)

}

plot.dz.by.cluster <- function(diagnoses,partition,dz.short,clusterColors,leg.lab){
  # diagnoses: length n vector of diagnoses 
  # partition: length n vector of cluster membership
  # dz.short: vector of short labels of diseases
  # clusterColors: vector of colors for each of 1:k clusters in partition
  # leg.lab: title for legend

  k <- length(unique(partition))
  dx <- data.frame(NPDx = diagnoses,stringsAsFactors = F)
  dx <- dummy.data.frame(dx,names = 'NPDx')
  colnames(dx) <- gsub('NPDx','',colnames(dx))  # format column names
  k.dz <- lapply(1:k, function(k.i) colMeans(dx[partition==k.i,]))
  k.dz <- do.call('cbind',k.dz)

  df.plt <- data.frame(pr = as.vector(k.dz),
                       Dz = rep(rownames(k.dz),k),
                       Dz.short = rep(dz.short,k),
                       cl = as.vector(sapply(1:k, function(k.i) rep(paste('Cluster',k.i),ncol(dx)))))
  th <- 0.05 # only label diseases with > 0.5 % of a cluster
  df.plt$Dz.short[df.plt$pr < th] <- ''
  pal <- colorRampPalette(brewer.pal(name = 'Set3',n=12))
  dz.colors <- pal(length(unique(df.plt$Dz)))

  n.by.cluster <- count.ejc(partition)
  clusterNames <- sapply(1:k, function(k.i) paste('Cluster',k.i))
  p.k.dz <- ggplot(data=df.plt) + geom_col(aes(x=cl,y=pr,fill=Dz),position=position_stack(1)) +
    geom_text(aes(x=cl,y=pr,label=Dz.short,group = Dz),position=position_stack(vjust=.5),size=1.5) +
    annotate(geom='text',x=clusterNames,y=rep(1.1,k),label=paste('',n.by.cluster,sep='\n'),size=1.5)+
    scale_y_continuous(expand=c(0,0)) + xlab('') + ylab('Proportion') +
    scale_fill_manual(values = dz.colors,name=leg.lab) +
    theme_classic() + theme(text = element_text(size=6),legend.key.size = unit(0.01,'in')) +
    theme(axis.text.x=element_text(angle=90,vjust=0.5,color=clusterColors))
  return(p.k.dz)
}

plot.allele.beta.matrix <- function(b.mat,p.mat,g.i,a.i,min.beta,max.beta,clusterColors){
    # plots matrix of pairwise comparison beta weights
    clusterNames <- rownames(b.mat)
    melted_betas <- melt(t(b.mat))
    melted_betas$Var2 <- fliplr(melted_betas$Var2)  
    p <- p.signif.matrix(p.mat)
    melted_p <- melt(t(p))
    melted_p$Var2 <- fliplr(melted_p$Var2)
    p1 <- ggplot() + 
      geom_tile(data = melted_betas, aes(x=Var1, y=Var2, fill=value)) + xlab("") + ylab("") +
      geom_text(data = melted_p, aes(x=Var1,y=Var2,label=value,color=value!='ns'),size=2.5) +
      scale_fill_gradientn(colours = c('dark red','#c23b22','white','#779ecb','dark blue'),
      #scale_fill_gradientn(colours = c('white','#779ecb','dark blue'),   # blue
                           values = rescale(c(min.beta,0,max.beta)),
                           guide = "colorbar", limits=c(min.beta,max.beta),
                           na.value = 'white',name=expression(beta)) +
      scale_color_manual(values=c('black','white'),limits=c(F,T),guide='none')+
      ggtitle(paste(g.i,a.i,sep = ':')) +
      scale_y_discrete(limits=clusterNames,labels = fliplr(clusterNames),expand=c(0,0)) + 
      scale_x_discrete(limits=clusterNames,labels = clusterNames,expand=c(0,0),position='top') + coord_fixed() +
      theme_classic() + theme(text=element_text(size = 6), legend.key.size = unit(0.1,'inches'),
            legend.position = 'bottom',
            axis.line = element_blank(),
            axis.text.x=element_text(color=clusterColors,angle=90),
            axis.text.y = element_text(color=fliplr(clusterColors)),
            plot.title = element_text(size=8,hjust=0.5,face = 'bold'),
            plot.margin = unit(c(0, 0, 0, 0), "cm"),
            axis.ticks = element_blank())
    return(p1)
}

colorswap <- function(colors,met){
  # if input colors are insufficient, generate a new palette. otherwise return colors input
  if(length(colors) < length(met)){
    colors <- colorRampPalette(brewer.pal(name = 'Set3',n=12))(length(met))
  }
  return(colors)
}

plot.model.perf.met <- function(met,perf.met,colors,ttl=''){
  # met: list by class of vectors of a performance metric over n reps of train test splits
  # perf.met: name of performance metric for plot label, i.e. sensitivity, specificity
  # colors: character vector of hex list. length should equal number of classes
  # compute mean and 95% CI (1.96*standard error)
  df.plt <- lapply(met, function(a)
    c(mean(a,na.rm=T),mean(a,na.rm=T) + 1.96*sd(a,na.rm=T)/sqrt(length(a)),mean(a,na.rm=T) - 1.96*sd(a,na.rm=T)/sqrt(length(a))))
  df.plt <- as.data.frame(do.call('rbind',df.plt))
  df.plt$dz <- rownames(df.plt)
  names(df.plt) <- c('met.mean','met.ul','met.ll','dz')
  df.plt$lab <- paste(signif(df.plt$met.ll,2),'-',signif(df.plt$met.ul,2),sep='')

  colors <- colorswap(colors,met)

  p.dx.dz <- ggplot(data=df.plt) + geom_col(aes(x=dz,y=met.mean,fill=dz)) + 
    geom_errorbar(aes(x=dz,ymin=met.ll,ymax=met.ul)) +
    geom_text(aes(x=dz,y=0.25,label=paste(signif(met.mean,2),'\n',lab)),size=2,hjust=0.5) +
    scale_y_continuous(limits=c(0,1),expand=c(0,0)) + coord_flip() +
    scale_fill_manual(values=colors) +
    xlab('') + ylab(perf.met) + ggtitle(ttl) +
    theme_classic() + theme(text = element_text(size=8), plot.title = element_text(hjust=0.5)) +
    theme(legend.position = 'none',plot.margin = unit(c(0, 0, 0, 0), "cm"),
      axis.text.y = element_text(color = colors))

  return(p.dx.dz)
}

plot.model.roc <- function(met,colors,ttl=''){

  colors <- colorswap(colors,met)
  auc <- lapply(met, function(a) sapply(a, function(c) c$AUC))
  which.auc <- sapply(auc, function(a) which.min(abs(a-mean(a)))) # find AUC value closest to mean AUC  
  roc.objs <- mapply(function(X,Y) {list(X[[Y]]$roc.obj)}, X=met,Y=which.auc)
  p.roc <-  ggroc(roc.objs,alpha=0.6,size=1,legacy.axes=TRUE) +
    geom_abline(slope=1,intercept=0,size=0.5,linetype='dashed',alpha= 0.5) +
    theme_classic() + scale_color_manual(limits =names(met),values= colors) +
    ylab('TPR') + xlab('FPR') + ggtitle(ttl) + 
    theme(text = element_text(size=8), plot.title = element_text(hjust=0.5),
      legend.position = 'none', plot.margin = unit(c(0, 0, 0, 0), "cm"),
      axis.title.y=element_text(vjust=-1))
  return(p.roc)
}

remove.y.ticklabels <- function(p){
  # remove y ticklabels and tickmarks from ggplot object p
  p <- p + theme(axis.text.y = element_blank(),axis.ticks.y=element_blank())
  return(p)
}

get.continuous.mask <- function(x){
  # x is a data frame which may contain both continuous and interval/categorical predictors
  # For my purposes, a continuous predictor is one with more than 3 levels 
  # (allele tables have 3 levels, 0,1,2)
  return(sapply(1:ncol(x), function(i) length(unique(x[,i]))) > 3)
}

plot.featureweights.lm.contintsplit <- function(res,df,colors,ttl = 'Feature Weights'){
  # plot feature weight heat map, separating continuous predictors (i.e. CSF) 
  # from interval predictors (i.e. allele counts)

  # identify continuous predictors

  cont.mask <- get.continuous.mask(df[,-1]) # remove INDDID column, get continuous predictors
  n.cont <- sum(cont.mask)
  n.pred <- length(cont.mask)
  # make cluster-by-feature heat map of weights
  weight <- t(do.call('cbind',lapply(res, function(C) # [,-1] removes intercept
        rowMeans(do.call('cbind',lapply(C, function(C.i) as.matrix(C.i$FeatureImportance)))))))[,-1]

  colors <- colorswap(colors,res)
  # generate feature weight heat maps separately for continuous and interval predictors
  # the two are not comparable at all and thus should not be plotted on the same color axis

  if(n.cont < n.pred & n.cont > 0){ # if there are both interval and continuous predictors
    p.cont <- plot.featureweights.lm(weight[,cont.mask],colors,ttl)
    p.int <- plot.featureweights.lm(weight[,!cont.mask],colors,ttl)

    p.all <- plot_grid(plotlist = list(p.cont,p.int), align = 'hv',nrow=1,axis = 'b')
    return(p.all)
  } else if(n.cont == n.pred){ # if there are only continuous predictors
    # i.e. csf only
    return(plot.featureweights.lm(weight[,cont.mask],colors,ttl))
  } else if(n.cont == 0){ # if there are only interval predictors
    # i.e. gene only
    return(plot.featureweights.lm(weight[,!cont.mask],colors,ttl))
  }


}

plot.featureweights.lm <- function(weight,colors,ttl){
# plot feature weight heat map
# weight is dz label-by-feature heat map of beta weights

melted_w <- melt(t(weight))

min.w <- min(melted_w$value)
max.w <- max(melted_w$value)

p1 <- ggplot() + 
    geom_tile(data = melted_w, aes(x=Var1, y=Var2, fill=value)) + xlab("") + ylab("") +
    scale_fill_gradientn(colours = c('dark red','#c23b22','white','#779ecb','dark blue'),
                           values = rescale(c(min.w,0,max.w)),
                           guide = "colorbar", limits=c(min.w,max.w),
                           na.value = 'white',name=expression(beta)) +
    scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
    ggtitle(ttl) + theme_classic() + 
    theme(text = element_text(size=8), plot.title = element_text(hjust=0.5,size=8),
    legend.key.size=unit(0.3,'cm'),axis.text.x = element_text(hjust=0.5,vjust=0.5,angle=90,size=6),
    axis.text.y = element_text(color=colors,size=6),axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1))

return(p1)

}


plot.featureweights.rf <- function(res,colors,ttl = 'Feature Weights'){
  # plot feature weight heat map
  
  # weight is dz label-by-feature heat map of mean decrease in accuracy
  weight <- t(do.call('cbind',lapply(res, function(C)
      rowMeans(do.call('cbind',lapply(C, function(C.i) as.matrix(C.i$FeatureImportance)))))))
  melted_w <- melt(t(weight))

  min.w <- min(melted_w$value)
  max.w <- max(melted_w$value)
  p1 <- ggplot() + 
      geom_tile(data = melted_w, aes(x=Var1, y=Var2, fill=value)) + xlab("") + ylab("") +
      scale_fill_gradientn(colours = c('white','#779ecb','dark blue'),
                             values = rescale(c(0,max.w)),
                             guide = "colorbar", limits=c(0,max.w),
                             na.value = 'white',name=expression(delta)) +
      ggtitle(ttl) + theme_classic() + 
      theme(text = element_text(size=8), plot.title = element_text(hjust=0.5),
      legend.key.size=unit(0.4,'cm'),axis.text.x = element_text(hjust=0.5,vjust=0.5,angle=90),
      axis.text.y = element_text(color=colors), plot.margin = unit(c(0, 0, 0, 0), "cm"))

  return(p1)
}

plot.beta.bar <- function(p.vals,ul.ll,m,ttl,fl){

  # plot beta and 95% CI for bootstrapped regression coefficients, excluding intercept
  df <- data.frame(n=names(coefficients(m)[-1]),b=coefficients(m)[-1],ul=ul.ll[,'97.5%'],ll=ul.ll[,'2.5%'])  
  sig <- p.signif(p.vals)
  p <- ggplot(df,aes(x=n,y=b)) + geom_col(fill=fl,color='black',alpha=0.3) +
    geom_errorbar(aes(x=n,ymin=ll,ymax=ul),width=0.5)+
    ylab(expression(beta)) + xlab('') + ggtitle(ttl) + theme_classic() +
    annotate(geom='text',x=df$n,y= 1.1*max(df$ul),label=sig,color='black') +
    theme(text=element_text(size=8),plot.title = element_text(size=8,hjust=0.5))
  return(p)
}