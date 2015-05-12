from diagnostic_reading import ShellSlice, AzAverage
import matplotlib.pyplot as plt
from matplotlib import ticker, font_manager
import pylab as p 
import numpy as np
from streamfunction import *
def get_lims(arr,boundstype='minmax',boundsfactor=1,themin=True):
    if (themin):
        if (boundstype == 'minmax'):
            val=arr.min()*boundsfactor
        elif (boundstype == 'rms'):
            val = -np.std(arr)*boundsfactor
    else:
        if(boundstype == 'minmax'):
            val=max(abs(arr.min()),arr.max())*boundsfactor
        elif (boundstype == 'rms'):
            val = np.std(arr)*boundsfactor
    return val
def plot_azav(fig,axis,field,radius,costheta,sintheta,r_bcz=0.71,mini=-1,maxi=-1,mycmap='jet',cbar=True, 
    boundsfactor = 1, boundstype = 'minmax', units = '',fontsize = 12, underlay = [0], nlevs = 6):
    #Modified version of Antoine Strukarek's routine
    r = radius/6.9599e10
    n_r=len(r)
    n_t=len(costheta)
    rtmp = r.reshape(1,n_r)
    cthtmp = costheta.reshape(n_t,1)
    sthtmp = sintheta.reshape(n_t,1)
    xr = p.dot(cthtmp,rtmp)
    yr = p.dot(sthtmp,rtmp)

    if (mini == -1):
        mini = get_lims(field,boundsfactor=boundsfactor,boundstype=boundstype,themin = True)
    if (maxi == -1):
        maxi = get_lims(field,boundsfactor=boundsfactor,boundstype=boundstype,themin = False)
    if (len(underlay)!=1):
        umini = get_lims(underlay,boundsfactor=boundsfactor,boundstype=boundstype,themin = True)
        umaxi = get_lims(underlay,boundsfactor=boundsfactor,boundstype=boundstype,themin = False)


    plt.hold(True)
    if (len(underlay) == 1):
        img = plt.pcolormesh(yr,xr,field,cmap=mycmap)
    else:
        img = plt.pcolormesh(yr,xr,underlay,cmap=mycmap)
    plt.plot(r_bcz*sintheta,r_bcz*costheta,'k--',[0,1],[0,0],'k--')
    plt.axis('equal')
    plt.axis('off')
    if (cbar):
        cbar = fig.colorbar(img,orientation='horizontal', shrink=0.5, aspect = 15)
        cbar.set_label(units)
        
        tick_locator = ticker.MaxNLocator(nbins=5)
        cbar.locator = tick_locator
        cbar.update_ticks()
        cbar.ax.tick_params(labelsize=fontsize)   #font size for the ticks

        t = cbar.ax.xaxis.label
        t.set_fontsize(fontsize)  # font size for the axis title
    if (len(underlay) == 1):
        plt.clim((mini,maxi))
    else:
        plt.clim((umini,umaxi))
    plt.xlim((0,1))
    plt.ylim((-1,1))
    plt.plot(r[0]*sintheta,r[0]*costheta,'k')
    plt.plot(r[n_r-1]*sintheta,r[n_r-1]*costheta,'k')
    plt.plot([0,0],[-r[n_r-1],r[n_r-1]],'k--')
    plt.plot([0,0],[-r[0],-r[n_r-1]],'k',[0,0],[r[n_r-1],r[0]],'k')
    plt.hold(False)

    plt.hold(True)
    levs=mini+np.linspace(1,nlevs,nlevs)/float(nlevs)*(maxi-mini)
    plt.contour(yr,xr,field,colors='w',levels=levs)
    plt.hold(False)


########################################################################
# Read in the data
b = AzAverage(filename='00095000',path='AZ_Avgs/')
print b.qv

       
ind = 1 # record number to grab



n_r = b.nr
n_t = b.ntheta
vphi = b.vals[:,:,b.lut[3],ind].reshape(n_t,n_r)
entropy = b.vals[:,:,b.lut[4],ind].reshape(n_t,n_r)

rhovr = b.vals[:,:,b.lut[13],ind].reshape(n_t,n_r)
rhovt = b.vals[:,:,b.lut[14],ind].reshape(n_t,n_r)


sintheta = b.sintheta
costheta = np.cos(np.arcsin(sintheta))


costheta[0:n_t/2] = -costheta[0:n_t/2] #sintheta is symmetric, so arcsin gives redundant values across equator
#costheta grid in Rayleigh runs from -1 to 1 (south pole to north pole..)



radius = b.radius
Omega=np.zeros((n_t,n_r))
GlobRate = 2.6e-6
for i in range(n_r):
    entropy[:,i]=entropy[:,i] - np.mean(entropy[:,i])
for i in range(n_r):
    Omega[:,i]=vphi[:,i]/(radius[i]*sintheta[:]) #+GlobRate
Omega = Omega*1.e9/(2.*np.pi) # s^-1 -> nHz

psi = streamfunction(rhovr,rhovt,radius,costheta,order=0)
rhovm = np.sqrt(rhovr**2+rhovt**2)*np.sign(psi)
print np.max(rhovm)
print np.max(psi)
print np.min(psi)
##############################################
#   Here we set up the actual figure.
#   We do a single row of 3 images 
#   Spacing is default spacing set up by subplot
f1 = p.figure(figsize=(5.5*3, 3*3), dpi=80)

tsize = 24
lsize = 15
ax1 = f1.add_subplot(1,3,1)
units = r'erg g$^{-1}$ K$^{-1}$'
plot_azav(f1,ax1,entropy,radius,costheta,sintheta,mycmap='RdYlBu_r',boundsfactor = 1.5, boundstype='rms', units=units, fontsize = lsize)
plt.title('Entropy',fontsize=tsize)

ax1 = f1.add_subplot(1,3,2)
units = 'nHz'
plot_azav(f1,ax1,Omega,radius,costheta,sintheta,mycmap='RdYlBu_r',boundsfactor = 0.05, units=units, fontsize = lsize)
plt.title(r'$\Omega$',fontsize=tsize)

ax1 = f1.add_subplot(1,3,3)
units = r'g cm$^{-2}$ s$^{-1}$'
plot_azav(f1,ax1,psi,radius,costheta,sintheta,mycmap='RdYlBu_r',boundsfactor = 1.5, boundstype='rms', units=units, fontsize = lsize, underlay = rhovm)
plt.title('Mass Flux',fontsize = tsize)

p.savefig('az_test.png')  
#plt.show()

