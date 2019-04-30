#########################
### Copath by cluster ###
#########################

pathItems.type <- list("Tau","Thio","Antibody","TDP43","Syn")
nTypes <- length(pathItems.type)
pathItems.index <- lapply(1:nTypes, function(i) grep(pathItems.type[[i]], colnames(microSample)))
type.by.patient <- lapply(1:nrow(microSample), function(Pt)
  do.call('cbind',lapply(pathItems.index, function(idx) as.numeric(microSample[Pt,idx]))))

copath.by.patient <- # flatten copath. cor-mat by row -> df
  do.call('rbind',lapply(type.by.patient, function(X) as.vector(t(cor(X)))))
copath.names <- sapply(pathItems.type, function(A) # name each element for combo
  sapply(pathItems.type, function(B) paste(A,'-',B,sep='')))
colnames(copath.by.patient) <- as.vector(t(copath.names))
#get rid of diagonal (1) and keep only utri (ltri is duplicate because symmetric cor-mat)
utri <- which(as.vector(t(upper.tri(matrix(1,nTypes,nTypes)))))
copath.by.patient <- as.data.frame(copath.by.patient[,utri])

df.plt <- data.frame(cl = as.character(as.vector(sapply(copath.by.patient, function(x) partition))),
  y = as.vector(as.matrix(copath.by.patient)),
  x = as.character(as.vector(sapply(colnames(copath.by.patient), function(X) rep(X,nrow(copath.by.patient))))))

p <- ggplot(data=df.plt) + geom_boxplot(aes(x=x,fill=cl,y=y),size=0.1) +
  theme_classic() + xlab('') + ylab('Copathology (r)') +
  theme(text=element_text(size=8),axis.text.x = element_text(angle=90,vjust=0.5,hjust=1))
p

ggsave(filename = paste(savedir,'CopathologyByClusterLouvain.pdf',sep=''),plot = p,
       height = 2,width=3,units='in')