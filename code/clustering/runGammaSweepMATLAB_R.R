rm(list=setdiff(ls(),'params'))

# Compute correlation matrix
mat.cmd <- paste(params$matlab.path,' -nodisplay -r "cd(\'',params$homedir,'code/clustering/\'); homedir = \'',params$homedir,'\'; ',
                 'opdir = \'',params$opdir,'\'; BCT_path =\'',params$BCT.path,'\'; run(\'computeSubjCorrMat.m\'); exit;"',sep='')
print('Sending to shell:')
print(mat.cmd)
print('....')
system(mat.cmd)

# Run Louvain clustering and consensus selection of partition
mat.cmd <- paste(params$matlab.path,' -nodisplay -r "cd(\'',params$homedir,'code/clustering/\'); homedir = \'',params$homedir,'\'; ',
                 'opdir = \'',params$opdir,'\'; BCT_path =\'',params$BCT.path,'\'; nreps = ',params$nreps_gammasweep,'; run(\'npsubjcommPartitionQuality.m\'); exit;"',sep='')
print('Sending to shell:')
print(mat.cmd)
print('....')
system(mat.cmd)