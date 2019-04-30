rm(list=setdiff(ls(),'params'))

load(file = paste(params$resultsdir,'analyzecluster/subjLouvainPartitionReordered.RData',sep=''))

# Run brain plotting scricpts
mat.cmd <- paste(params$matlab.path,' -nodisplay -r "cd(\'',params$homedir,'code/plot_brains/\'); homedir = \'',params$homedir,'\'; ',
                 'opdir = \'',params$opdir,'\'; k = ',k,'; resultsdir = \'',params$resultsdir,
                 '\'; BCT_path =\'',params$BCT.path,'\'; run(\'plot_brains.m\'); exit;"',sep='')
print('Sending to shell:')
print(mat.cmd)
print('....')
system(mat.cmd)