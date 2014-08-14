from diagnostic_reading import ShellSlice
import matplotlib.pyplot as plt
from matplotlib import ticker
a = ShellSlice(filename='00280010',path='Shell_Slices/')

var_inds = [0,2,3]
rad_inds = [1,2,4]
#Tex can be enclosed in dollar signs within a string.  The r in front of the string is necessary for strings enclosing Tex
units = [r'cm s$^{-1}$', r'cm s$^{-1}$', r'erg g$^{-1}$ K$^{-1}$']  
vnames = [r'v$_r$', r'v$_\phi$', "S'"]
ncol = len(var_inds)
nrow = len(rad_inds)
nplots = ncol*nrow
        


import pylab as p 
import numpy as np
ind = 1
f1 = p.figure(figsize=(5.5*3, 5*3), dpi=80)

for  j in range(nrow):
    for i in range(ncol):
        slice = a.vals[:,:,rad_inds[j],var_inds[i]].reshape(a.nphi,a.ntheta)
        slice = slice-np.mean(slice)
        slice = np.transpose(slice)
        
        ax1 = f1.add_subplot(nrow,ncol,ind, projection="mollweide")
        
        ind = ind+1


        twosigma = 2*np.std(slice)
        print "scaling from +- ", twosigma
        contour_levels = twosigma*np.linspace(-1,1,256)
        image1 = ax1.imshow(slice,vmin=-twosigma, vmax=twosigma, extent=(-np.pi,np.pi,-np.pi/2,np.pi/2), clip_on=False, aspect=0.5, interpolation='bicubic')
        image1.axes.get_xaxis().set_visible(False)  # Remove the longitude and latitude axis labels
        image1.axes.get_yaxis().set_visible(False)
        image1.set_cmap('RdYlBu_r')  # Red/Blue map
        #image1.set_cmap('gnuplot2')  # "ASH-Like" yellow/violet map

        cbar = f1.colorbar(image1,orientation='horizontal', shrink=0.5)
        cbar.set_label(units[i])

        tick_locator = ticker.MaxNLocator(nbins=5)
        cbar.locator = tick_locator
        cbar.update_ticks()

pbottom = 0.05
pright = 0.95
pleft = 0.15
ptop = 0.95
phspace = 0.1
pwspace = 0.03
# These parameters can also be set interactively when using plt.show()
p.subplots_adjust(left = pleft, bottom = pbottom, right = pright, top = ptop, wspace=pwspace,hspace=phspace) 

# Add the radial information to the left margin
rspace = (ptop-pbottom)/nrow
rnorm = 6.96e10
lilbit = 0.05
xpos = 0.1*pleft
for i in range(nrow):
    ypos = ptop-rspace*0.4-i*(rspace+phspace/(nrow-1))
    ratio = float(a.radius[i]/rnorm)
    r_str = r'r/R$_\odot$ = %1.3f' % ratio
    f1.text(xpos,ypos,r_str,fontsize=16)
    ypos = ypos-lilbit
    ind_str = 'radial index: '+str(a.inds[i])
    f1.text(xpos,ypos,ind_str,fontsize=16)

cspace = (pright-pleft)/ncol
for i in range(ncol):
    ypos = ptop
    xpos = pleft+cspace*0.47+i*cspace 
    f1.text(xpos,ypos,vnames[i],fontsize=20)
p.savefig('shell_test.png')  
#plt.show()

