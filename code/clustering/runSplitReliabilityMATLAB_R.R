rm(list=setdiff(ls(),c('params','sampfrac')))

# Compute correlation matrix
mat.cmd <- paste(params$matlab.path,' -nodisplay -r "cd(\'',params$homedir,'code/clustering/\'); homedir = \'',params$homedir,'\'; ',
                 'opdir = \'',params$opdir,'\'; sampfrac = ',sampfrac,'; gamma = ',params$gamma.opt,'; BCT_path =\'',params$BCT.path,'\'; run(\'splitreliability_v2.m\'); exit;"',sep='')
print('Sending to shell:')
print(mat.cmd)
print('....')
system(mat.cmd)