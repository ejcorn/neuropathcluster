library(colorspace)
library(R.matlab)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
homedir <- '~/Dropbox/Neurodegeneration/PathCogClinDx/neuropathcluster/code/plot_brains/schematic'
setwd(homedir)

roi_names <- list('Anterior Cingulate Cortex','Superior-middle temporal cortex',
                  'Middle frontal cortex','Angular gyrus','Entorhinal Cortex',
                  'Thalamus','Caudate-putamen','Globus pallidus','Hippocampus','Amygdala',
                  'Midbrain','Pons','Medulla','Cerebellum')

cort_pal <- colorRampPalette(brewer.pal(name = 'Accent',n=5))
subcort_pal <- colorRampPalette(brewer.pal(name = 'Set1',n=9))
roi_hexcolors <- c(cort_pal(5),subcort_pal(9))
roi_hexcolors[roi_names == 'Caudate-putamen'] <- '#afeeee' # manually set this one to distinguish it
df <- data.frame(lab = unlist(roi_names), y = length(roi_names):1,col = roi_hexcolors,stringsAsFactors = F)
p <- ggplot(data=df) + geom_text(aes(x=0,y=y,label=lab,color=lab),size=2.25) +
  scale_color_manual(limits=df$lab,values = roi_hexcolors) + theme_classic()+
  theme(legend.position = 'none',line=element_blank(),text=element_blank())
p
ggsave(plot=p,filename = 'region_names.pdf',units='in',height=2,width=1.25)

roi_colors <- hex2RGB(roi_hexcolors) @ coords
# duplicate color for caudate-putamen so they can be plotted as the same color
roi_colors <- rbind(roi_colors[1:which(roi_names == 'Caudate-putamen'),],roi_colors[which(roi_names == 'Caudate-putamen'),],
                    roi_colors[(1+which(roi_names == 'Caudate-putamen')):length(roi_names),])
writeMat(con='regioncolors.mat',roi_colors=roi_colors)


hm.pal <- colorRampPalette(brewer.pal(name = 'Blues',n=5))
r <- 20
X <- lapply(1:7, function(c) round(runif(r,1,5))) # basic values for each pathological protein
X <- t(do.call('cbind',lapply(X, function(x) sapply(1:14, function(j) x + round(rnorm(r,mean=1,sd=0.4))))))
X.m <- melt(t(X))
X.m$value_discrete <- rep(NA,r)
X.m$value_discrete[X.m$value == 1] <- '0'
X.m$value_discrete[X.m$value == 2] <- 'Rare'
X.m$value_discrete[X.m$value == 3] <- '1+'
X.m$value_discrete[X.m$value == 4] <- '2+'
X.m$value_discrete[X.m$value == 5] <- '3+'

p <- ggplot(data=X.m,aes(x=Var1,y=Var2,fill=value_discrete)) + geom_tile() +
  scale_fill_manual(values=hm.pal(5)) + ylab('') + xlab('') + theme_classic()+
  theme(axis.line = element_blank(),axis.ticks=element_blank(),
        axis.text.x = element_blank(), axis.text.y=element_blank(),
        legend.title = element_blank(),legend.key.size = unit(0.1,'cm'),
        legend.text = element_text(size=6))
p
ggsave(plot=p,filename = 'mock_datamatrix.pdf',units='cm',width=5.75,height=4.75)

## mock similarity matrix

X.r <- melt(t(cor(X)))
p <- ggplot(data=X.r,aes(x=Var1,y=Var2,fill=value)) + geom_tile() +
  scale_fill_gradientn(colours = c('#8B0000','#c23b22','#ffffff','#779ecb','#00008b'),
                       values = rescale(c(-1,0,1)), guide = "colorbar",limits=c(-1,1),breaks = c(-1,0,1),name='r') +
  ylab('') + xlab('') + coord_fixed() +
  theme(axis.line = element_blank(),axis.ticks=element_blank(),
        axis.text.x = element_blank(), axis.text.y=element_blank(),
        legend.title = element_blank(),legend.key.size = unit(0.1,'cm'),
        legend.text = element_text(size=6))
p
ggsave(plot=p,filename = 'mock_similaritymatrix.pdf',units='cm',width=4,height=4)


## plot number of patients in each cluster
load(file='~/Dropbox/Neurodegeneration/PathCogClinDx/neuropathcluster/neuropathcluster_R0.2C0_010320/results_G6/analyzecluster/subjLouvainPartitionReordered.RData')
df <- data.frame(y=count.ejc(partition),x=names(count.ejc(partition)),stringsAsFactors = F)
clusterColors <- getClusterColors(max(df$x))
p <- ggplot(df) + geom_col(aes(x=x,y=y,fill=x)) + 
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=clusterColors)+ 
  annotate(geom='text',x = df$x,y=20,label=as.character(df$x),size=2)+theme_classic()+
  theme(axis.text = element_blank(),axis.ticks = element_blank(),legend.position = 'none',
        text=element_text(size=6)) +
  xlab('Cluster') + ylab('# of Patients')
ggsave(plot=p,filename = 'mock_ClusterMembership.pdf',units='cm',width=2,height=2)


## make table with region names
roi_abbrevs <- c('Cing','SMT','MF','Ang','CS','EC','Amyg','TS','CP','GP','SN','Med','CB','Pons','MB')
roi_long <- c('Anterior cingulate cortex','Superior-middle temporal cortex',
                  'Middle frontal cortex','Angular gyrus','CA1/Subiculum','Entorhinal cortex',
                  'Amygdala','Thalamus','Caudate-putamen','Globus pallidus','Substantia nigra',
                  'Medulla','Cerebellum','Pons','Midbrain')
x <- data.frame(Abbreviation=roi_abbrevs,Region=roi_long)
#rownames(x) <-NULL
print(xtable(x,caption='Definition of brain region abbreviations.',label='table:brainregions'), include.rownames=FALSE)