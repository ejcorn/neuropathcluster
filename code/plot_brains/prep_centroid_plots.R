rm(list = setdiff(ls(), "params"))
homedir <- params$homedir
setwd(homedir)
savedir <- paste(params$resultsdir,'plot_brains/vals/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fxns.R')

######################
### Load centroids ###
######################

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))

pathItems.type <- list("NeuronLoss","Gliosis","Angiopathy","Ubiquitin","Thio","TDP43","Tau","Syn","Antibody")
names(pathItems.type) <- unlist(pathItems.type)
pathRegions.name <- list("Ci","OC","SM","MF","An","EC","CS","DG","Am","TS","CP","GP","Me","CB","Po","MB")
pathRegions.name <- list(Ci=1,OC=2,SM=3,MF=4,An=5,EC=6,CS=114,DG=114,Am=115,TS=109,CP=110,PU=111,GP=112,Me=236,Po=235,MB=234,CB=237)
pathRegions.prettylab <- list("Cing","OC","SM","MF","Ang","CA1/Sub","EC","DG","Amyg","TS","CP","GP","Med","CB","Pons","Mesenc.")

# Make empty matrix to hold amount of path for a given type for a given centroid
# with ordered regions. 2nd column is index in parcellation
C.region.plot <- matrix(NA,nrow=length(pathRegions.name),ncol=2)
rownames(C.region.plot) <- names(pathRegions.name)
colnames(C.region.plot) <- c('Region','Index')
C.region.plot[,2] <- unlist(pathRegions.name)

# iterate through clusters, extract weights on each region for each path feature

for(C in rownames(centroids)){
	centroid <- centroids[C,]
	centroid.by.type <- lapply(pathItems.type, function(p.type) centroid[grep(p.type,names(centroid))])
	for(p.type in names(centroid.by.type)){
		C.t <- centroid.by.type[[p.type]]
		# retrieve empty matrix
		C.region.plot.t <- C.region.plot
		for(p.region in names(pathRegions.name)){
			# store path for each region for given type in preset matrix
			# only search in first couple letters of each name
			# otherwise you get everything with Angiopathy for An (angular gyrus)
			region.path <- C.t[grep(p.region,substr(names(C.t),start=1,stop=2))]
			if(!is_empty(region.path)){
				C.region.plot.t[p.region,1] <- region.path
			}
			
		}
		# average DG and CA1/Subiculum to get one value for whole hippocampus	
		C.region.plot.t['DG',1] <- mean(c(C.region.plot.t['DG',1],C.region.plot.t['CS',1]))		
 		# remove CS, store whole hippocampus in DG
		C.region.plot.t <- C.region.plot.t[-grep('CS',rownames(C.region.plot.t)),]
		# duplicate caudate and store as putamen to plot together
		C.region.plot.t['PU',] <- c(C.region.plot.t['CP',1],111)
		print(sum(is.na(C.region.plot.t)))
		writeMat(con=paste(savedir,C,'_',p.type,'.mat',sep=''),C.region.plot.t=C.region.plot.t)
	}
}