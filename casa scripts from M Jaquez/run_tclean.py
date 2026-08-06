import time

t0 = time.time()
line = 'H38'
def run_tclean(msfile,n,r):
    if n == 0 :
       imagename = 'obs2hr_psf_model_n0_r{}'.format(r)
    else:
        imagename = 'obs2hr_disk_{}_line_r{}_n{}'.format(line,r,1)
        interactive = True
    #os.system('rm -rf '+imagename+'.*')
    tclean(vis=msfile, 
            datacolumn='', imagename=imagename, 
            imsize = [1200],startmodel='', 
            cell='0.5mas', specmode='cube',  
            gridder='standard', 
            deconvolver='multiscale',scales=[0,7,21], 
            weighting='briggs',robust=r, uvtaper='2.0mas', 
            threshold = '100.0uJy',niter=n, interactive=True)
    return n


#imsfile = 'disk_93GHz_line_{}_noisy.ms'.format(line)
msfile = ['obs1_disk_93GHz_line_H38_noisy.ms','obs4_disk_93GHz_line_H38_noisy.ms']
n =  0 #100000
r = - 0.5

run_tclean(msfile,n,r)

print ('the tclean finish in {} seconds or {} minuts'.format(time.time()-t0, (time.time()-t0)/60))

