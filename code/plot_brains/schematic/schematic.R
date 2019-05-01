library(colorspace)
library(R.matlab)
library(ggplot2)
library(reshape2)
homedir <- '~/Dropbox/Neurodegeneration/PathCogClinDx/neuropathcluster/code/plot_brains/schematic'
setwd(homedir)

roi_names <- list('Cingulate Cortex','Occipital Cortex','Superior-middle temporal cortex',
                  'Middle frontal cortex','Angular gyrus','Entorhinal Cortex',
                  'Thalamus','Caudate-putamen','Globus pallidus','Hippocampus','Amygdala',
                  'Midbrain','Pons','Medulla','Cerebellum')

cort_pal <- colorRampPalette(brewer.pal(name = 'Accent',n=6))
subcort_pal <- colorRampPalette(brewer.pal(name = 'Set1',n=9))
roi_hexcolors <- c(cort_pal(6),subcort_pal(9))
df <- data.frame(lab = unlist(roi_names), y = length(roi_names):1,col = roi_hexcolors,stringsAsFactors = F)
p <- ggplot(data=df) + geom_text(aes(x=0,y=y,label=lab,color=lab),size=1.5) +
  scale_color_manual(limits=df$lab,values = roi_hexcolors) +
  theme(legend.position = 'none',line=element_blank(),text=element_blank())
p
ggsave(plot=p,filename = 'region_names.pdf',units='in',height=1.6,width=1)

roi_colors <- hex2RGB(roi_hexcolors) @ coords
# duplicate color for caudate-putamen so they can be plotted as the same color
roi_colors <- rbind(roi_colors[1:which(roi_names == 'Caudate-putamen'),],roi_colors[which(roi_names == 'Caudate-putamen'),],
                    roi_colors[(1+which(roi_names == 'Caudate-putamen')):length(roi_names),])
writeMat(con='regioncolors.mat',roi_colors=roi_colors)


hm.pal <- colorRampPalette(brewer.pal(name = 'Blues',n=5))
r <- 20
X <- t(do.call('cbind',lapply(1:9, function(c) round(runif(r,1,5)))))
X.m <- melt(t(X))
X.m$value_discrete <- rep(NA,r)
X.m$value_discrete[X.m$value == 1] <- '0'
X.m$value_discrete[X.m$value == 2] <- 'Rare'
X.m$value_discrete[X.m$value == 3] <- '1+'
X.m$value_discrete[X.m$value == 4] <- '2+'
X.m$value_discrete[X.m$value == 5] <- '3+'

p <- ggplot(data=X.m,aes(x=Var1,y=Var2,fill=value_discrete)) + geom_tile() +
  scale_fill_manual(values=hm.pal(5)) + ylab('') + xlab('') +
  theme(axis.line = element_blank(),axis.ticks=element_blank(),
        axis.text.x = element_blank(), axis.text.y=element_blank(),
        legend.title = element_blank(),legend.key.size = unit(0.1,'cm'),
        legend.text = element_text(size=6))
p
ggsave(plot=p,filename = 'mock_datamatrix.pdf',units='cm',width=10,height=6)

## mock similarity matrix

X.r <- melt(t(cor(X)))
p <- ggplot(data=X.r,aes(x=Var1,y=Var2,fill=value)) + geom_tile() +
  scale_fill_gradientn(colours = c('dark red','#c23b22','white','#779ecb','dark blue'),
                       values = rescale(c(-1,0,1)), guide = "colorbar",limits=c(-1,1),breaks = c(-1,0,1),name='r') +
  ylab('') + xlab('') + coord_fixed() +
  theme(axis.line = element_blank(),axis.ticks=element_blank(),
        axis.text.x = element_blank(), axis.text.y=element_blank(),
        legend.title = element_blank(),legend.key.size = unit(0.1,'cm'),
        legend.text = element_text(size=6))
p
ggsave(plot=p,filename = 'mock_similaritymatrix.pdf',units='cm',width=4,height=4)
