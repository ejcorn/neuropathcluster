rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'genecluster/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')
source('code/misc/plottingfxns.R')

microSample <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs
patientSample <- read.csv(paste(params$opdir,'processed/patientSample.csv',sep=''))[,-(1:2)] # Get rid of index column and INDDIDs

INDDIDs <- read.csv(paste(params$opdir,'processed/microSample.csv',sep=''))[,2]

if(sum(duplicated(INDDIDs))){
  break
}

#####################
### Load clusters ###
#####################

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))
INDDIDs <- remove.Disconnected.Subjects(INDDIDs,DisconnectedSubjects)
if(!is_empty(DisconnectedSubjects)){
  microSample <- microSample[-DisconnectedSubjects,]
  patientSample <- patientSample[-DisconnectedSubjects,]
}

clusterColors <- getClusterColors(k)

#################
### Gene data ###
#################

GOIs <- c('APOE','MAPTHaplotype')#,'C9orf72','LRRK2')
Genes <- read.csv(paste(params$opdir,'processed/Genetics_processed.csv',sep=''),stringsAsFactors = F)
Genes <- Genes[Genes$INDDID %in% INDDIDs,]
missing.mask <- rowSums(Genes[,GOIs] != '') == length(GOIs)
Genes <- Genes[missing.mask,]  # remove missing
Genes.df <- Genes[,GOIs] # select genes of interest

############################
### Process into alleles ###
############################

Alleles <- list(APOE=c('E2','E3','E4'),MAPTHaplotype=c('H1','H2'))
Genotypes <- lapply(Genes.df, function(X) unique(X))
# for each gene, get an Allele-by-Genotype table counting # of each allele per genotype
n.alleles.per.genotype <- list()
Allele.Tables <- list()
for(g.i in names(Genotypes)){
  n.alleles.per.genotype[[g.i]] <- 
    sapply(Genotypes[[g.i]], function(gt.i)
      sapply(Alleles[[g.i]], function(A) str_count(pattern = A,string = gt.i)))
  Allele.Tables[[g.i]] <- t(sapply(Genes.df[,g.i], function(gt.i) n.alleles.per.genotype[[g.i]][,gt.i]))
}

partitionSample <- partition[INDDIDs %in% unique(Genes$INDDID)]
clusterNames <- sapply(sort(unique(partitionSample)), function(i) paste('Cluster',i))

# Create pairwise comparisons using nx2 contingency tables, 
# where n is number of alleles, 2 is a pairwise intercluster comparison

results <- list()
for(g.i in names(Genotypes)){
  A <- Allele.Tables[[g.i]]
  results[[g.i]] <- list()
  for(k.1 in 1:k){
    results[[g.i]][[clusterNames[k.1]]] <- list()
    for(k.2 in 1:k){
      cont.tab <- rbind(colSums(A[partitionSample==k.1,]),colSums(A[partitionSample==k.2,]))
      rownames(cont.tab) <- clusterNames[c(k.1,k.2)]
      m <- fisher.test(cont.tab)
      results[[g.i]][[clusterNames[k.1]]][[clusterNames[k.2]]] <- m
    }
  }
}

# extract p-vals from unique tests and compute FDR threshold
cluster.pairs <- t(combn(1:k,2))
p.vals <- unlist(lapply(results, function(g.r) # loop through loci/genes
  as.vector(sapply(1:nrow(cluster.pairs), function(i) # loop through cluster comparisons
    g.r[[cluster.pairs[i,1]]][[cluster.pairs[i,2]]]$p.value))))
p.thrsh <-  compute.FDR(pvalue.vec = p.vals,q=0.05)

# extract difference in proportions for each allele from all tests 
p1p2 <- list()
pvals <- list()
for(g.i in names(results)){
  betas[[g.i]] <- list()
  pvals[[g.i]] <- list()
  g.r <- results[[g.i]]
  for(a.i in rownames(g.r[[1]][[1]])){
    # loop through cluster comparisons
    betas[[g.i]][[a.i]] <- t(sapply(1:k, function(k1) # transpose so that you have a matrix of log odds ratios, where ij is log odds cluster=i, not j
      sapply(1:k, function(k2) # set to 0 if p < FDR threshold
        g.r[[clusterNames[k1]]][[clusterNames[k2]]][a.i,'Estimate']*(g.r[[clusterNames[k1]]][[clusterNames[k2]]][a.i,'Pr(>|z|)'] < p.thrsh))))
    pvals[[g.i]][[a.i]] <- t(sapply(1:k, function(k1) # store p-vals to add significance asterisks
      sapply(1:k, function(k2) # set to 0 if p < FDR threshold
        g.r[[clusterNames[k1]]][[clusterNames[k2]]]$p.value)))
    colnames(betas[[g.i]][[a.i]]) <- 1:k
  }
}

# plot betas
max.beta <- max(as.vector(unlist(betas)),na.rm = T)
min.beta <- min(as.vector(unlist(betas)),na.rm=T)
plots <- list()
for(g.i in names(betas)){
  plots[[g.i]] <- list()
  for(a.i in names(betas[[g.i]])){
    t.mat <- betas[[g.i]][[a.i]]
    melted_betas <- melt(t(t.mat))
    melted_betas$Var2 <- fliplr(melted_betas$Var2)
    p <- p.signif.matrix(p = pvals[[g.i]][[a.i]]) # compute significance level
    melted_p <- melt(t(p))
    melted_p$Var2 <- fliplr(melted_p$Var2)
    p1 <- ggplot() + 
      geom_tile(data = melted_betas, aes(x=Var1, y=Var2, fill=value)) + xlab("") + ylab("") +
      geom_text(data = melted_p, aes(x=Var1,y=Var2,label=value,color=value!='ns'),size=1.5) +
      scale_fill_gradientn(colours = c('dark red','#c23b22','white','#779ecb','dark blue'),
                           values = rescale(c(min.beta,0,max.beta)),
                           guide = "colorbar", limits=c(min.beta,max.beta),
                           na.value = 'white',name=expression(beta)) +
      scale_color_manual(values=c('black','white'),limits=c(F,T),guide='none')+
      ggtitle(paste(g.i,a.i,sep = ':')) +
      scale_y_discrete(limits=1:k,labels = fliplr(clusterNames),expand=c(0,0)) + 
      scale_x_discrete(limits=1:k,labels = clusterNames,expand=c(0,0),position='top') + coord_fixed() +
      theme_classic() + theme(text=element_text(size = 6), legend.key.size = unit(0.1,'inches')) +
      theme(axis.line = element_blank(),axis.text.x=element_text(color=clusterColors),
            axis.text.y = element_text(color=fliplr(clusterColors)),
            plot.title = element_text(size=8,hjust=0.5,face = 'bold'),
            axis.ticks = element_blank())
    plots[[g.i]][[a.i]] <- p1
  }
  p.all <- plot_grid(plotlist = plots[[g.i]], align = 'hv',nrow=1)
  ggsave(filename = paste(savedir,g.i,'Betas.pdf',sep=''),plot = p.all,
         height = 2,width=2.5*length(plots[[g.i]]),units='in')
}