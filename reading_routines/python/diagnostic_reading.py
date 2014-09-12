import numpy as np
import os
class ReferenceState:
    """Rayleigh Reference State Structure
    ----------------------------------
    self.n_r     : number of radial points
    self.radius  : radial coordinates
    """

    def __init__(self,filename='none',path='./',esize=8):
        """filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename
           esize     : The size of each data element in bytes (default 8-byte doubles)"""
        if (filename == 'none'):
            the_file = path+'reference'
        else:
            the_file = path+filename
        fd = open(the_file,'rb')
        # We read an integer to assess which endian the file was written in...
        bs = check_endian(fd,314,'int32')
        
        nr = swapread(fd,dtype='int32',count=1,swap=bs)
        tmp = np.reshape(swapread(fd,dtype='float64',count=10*nr,swap=bs),(nr,10), order = 'F')
        self.nr = nr
        self.radius      = tmp[:,0]
        self.density     = tmp[:,1]
        self.dlnrho      = tmp[:,2]
        self.d2lnrho     = tmp[:,3]
        self.pressure    = tmp[:,4]
        self.temperature = tmp[:,5]
        self.dlnt        = tmp[:,6]
        self.dsdr        = tmp[:,7]
        self.entropy     = tmp[:,8]
        self.gravity     = tmp[:,9]
        self.ref = tmp
        self.names = ['radius', 'density', 'dlnrho', 'd2lnrho', 'pressure', 'temperature',
        'dlnt', 'dsdr','entropy','gravity']
        fd.close()

class GlobalAverage:
    """Rayleigh GlobalAverage Structure
    ----------------------------------
    self.n_r     : number of radial points
    self.radius  : radial coordinates
    """

    def __init__(self,filename='none',path='Shell_Avgs/',esize=8):
        """filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename
           esize     : The size of each data element in bytes (default 8-byte doubles)"""
        if (filename == 'none'):
            the_file = path+'00000001'
        else:
            the_file = path+filename
        fd = open(the_file,'rb')
        # We read an integer to assess which endian the file was written in...
        bs = check_endian(fd,314,'int32')
        version = swapread(fd,dtype='int32',count=1,swap=bs)
        nrec = swapread(fd,dtype='int32',count=1,swap=bs)
        nq = swapread(fd,dtype='int32',count=1,swap=bs)
        self.niter = nrec
        self.nq = nq
        self.qv = np.reshape(swapread(fd,dtype='int32',count=nq,swap=bs),(nq), order = 'F')
        self.vals  = np.zeros((nrec,nq),dtype='float64')
        self.iters = np.zeros(nrec,dtype='int32')
        self.time  = np.zeros(nrec,dtype='float64')
        self.version = version
        for i in range(nrec):
            tmp = np.reshape(swapread(fd,dtype='float64',count=nq,swap=bs),(nq), order = 'F')
            self.vals[i,:] = tmp
            self.time[i] = swapread(fd,dtype='float64',count=1,swap=bs)
            self.iters[i] = swapread(fd,dtype='int32',count=1,swap=bs)
        maxq = 250
        lut = np.zeros(maxq)+int(1000)
        self.lut = lut.astype('int32')
        for i,q in enumerate(self.qv):
            self.lut[q] = i
        fd.close()

class ShellAverage:
    """Rayleigh GlobalAverage Structure
    ----------------------------------
    self.n_r     : number of radial points
    self.radius  : radial coordinates
    """

    def __init__(self,filename='none',path='Shell_Avgs/',esize=8):
        """filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename
           esize     : The size of each data element in bytes (default 8-byte doubles)"""
        if (filename == 'none'):
            the_file = path+'00000001'
        else:
            the_file = path+filename
        fd = open(the_file,'rb')
        # We read an integer to assess which endian the file was written in...
        bs = check_endian(fd,314,'int32')
        version = swapread(fd,dtype='int32',count=1,swap=bs)
        nrec = swapread(fd,dtype='int32',count=1,swap=bs)
        nr = swapread(fd,dtype='int32',count=1,swap=bs)
        nq = swapread(fd,dtype='int32',count=1,swap=bs)

        self.version = version
        self.niter = nrec
        self.nq = nq
        self.nr = nr
        self.qv = np.reshape(swapread(fd,dtype='int32',count=nq,swap=bs),(nq), order = 'F')
        self.radius = np.reshape(swapread(fd,dtype='float64',count=nr,swap=bs),(nr), order = 'F')

        self.vals  = np.zeros((nr,nq,nrec),dtype='float64')
        self.iters = np.zeros(nrec,dtype='int32')
        self.time  = np.zeros(nrec,dtype='float64')

        for i in range(nrec):
            tmp = np.reshape(swapread(fd,dtype='float64',count=nq*nr,swap=bs),(nr,nq), order = 'F')
            self.vals[:,:,i] = tmp
            self.time[i] = swapread(fd,dtype='float64',count=1,swap=bs)
            self.iters[i] = swapread(fd,dtype='int32',count=1,swap=bs)
        maxq = 250
        lut = np.zeros(maxq)+int(1000)
        self.lut = lut.astype('int32')
        for i,q in enumerate(self.qv):
            self.lut[q] = i
        fd.close()

class AzAverage:
    """Rayleigh GlobalAverage Structure
    ----------------------------------
    self.n_r     : number of radial points
    self.radius  : radial coordinates
    """

    def __init__(self,filename='none',path='Shell_Avgs/',esize=8):
        """filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename
           esize     : The size of each data element in bytes (default 8-byte doubles)"""
        if (filename == 'none'):
            the_file = path+'00000001'
        else:
            the_file = path+filename
        fd = open(the_file,'rb')
        # We read an integer to assess which endian the file was written in...
        bs = check_endian(fd,314,'int32')
        version = swapread(fd,dtype='int32',count=1,swap=bs)
        nrec = swapread(fd,dtype='int32',count=1,swap=bs)
        nr = swapread(fd,dtype='int32',count=1,swap=bs)
        ntheta = swapread(fd,dtype='int32',count=1,swap=bs)
        nq = swapread(fd,dtype='int32',count=1,swap=bs)

        self.version = version
        self.niter = nrec
        self.nq = nq
        self.nr = nr
        self.ntheta = ntheta
        self.qv = np.reshape(swapread(fd,dtype='int32',count=nq,swap=bs),(nq), order = 'F')
        self.radius = np.reshape(swapread(fd,dtype='float64',count=nr,swap=bs),(nr), order = 'F')
        self.costheta = np.reshape(swapread(fd,dtype='float64',count=ntheta,swap=bs),(ntheta), order = 'F')
        self.sintheta = (1.0-self.costheta**2)**0.5
        self.vals  = np.zeros((nr,ntheta,nq,nrec),dtype='float64')
        self.iters = np.zeros(nrec,dtype='int32')
        self.time  = np.zeros(nrec,dtype='float64')

        for i in range(nrec):
            tmp = np.reshape(swapread(fd,dtype='float64',count=nq*nr*ntheta,swap=bs),(nr,ntheta,nq), order = 'F')
            self.vals[:,:,:,i] = tmp
            self.time[i] = swapread(fd,dtype='float64',count=1,swap=bs)
            self.iters[i] = swapread(fd,dtype='int32',count=1,swap=bs)
        maxq = 250
        lut = np.zeros(maxq)+int(1000)
        self.lut = lut.astype('int32')
        for i,q in enumerate(self.qv):
            self.lut[q] = i
        fd.close()

class ShellSlice:
    """Rayleigh GlobalAverage Structure
    ----------------------------------
    self.n_r     : number of radial points
    self.radius  : radial coordinates
    """

    def __init__(self,filename='none',path='Shell_Avgs/',esize=8):
        """filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename
           esize     : The size of each data element in bytes (default 8-byte doubles)"""
        if (filename == 'none'):
            the_file = path+'00000001'
        else:
            the_file = path+filename
        fd = open(the_file,'rb')
        # We read an integer to assess which endian the file was written in...
        bs = check_endian(fd,314,'int32')
        version = swapread(fd,dtype='int32',count=1,swap=bs)
        nrec = swapread(fd,dtype='int32',count=1,swap=bs)
        print bs, version, nrec
        ntheta = swapread(fd,dtype='int32',count=1,swap=bs)
        nphi = 2*ntheta
        nr = swapread(fd,dtype='int32',count=1,swap=bs)
        nq = swapread(fd,dtype='int32',count=1,swap=bs)

        self.niter = nrec
        self.nq = nq
        self.nr = nr
        self.ntheta = ntheta
        self.nphi = nphi

        self.qv = np.reshape(swapread(fd,dtype='int32',count=nq,swap=bs),(nq), order = 'F')
        self.radius = np.reshape(swapread(fd,dtype='float64',count=nr,swap=bs),(nr), order = 'F')
        self.inds = np.reshape(swapread(fd,dtype='int32',count=nr,swap=bs),(nr), order = 'F')
        self.costheta = np.reshape(swapread(fd,dtype='float64',count=ntheta,swap=bs),(ntheta), order = 'F')
        self.sintheta = (1.0-self.costheta**2)**0.5

        self.vals  = np.zeros((nphi,ntheta,nr,nq,nrec),dtype='float64')
        self.iters = np.zeros(nrec,dtype='int32')
        self.time  = np.zeros(nrec,dtype='float64')
        self.version = version
        for i in range(nrec):
            tmp = np.reshape(swapread(fd,dtype='float64',count=nq*nr*ntheta*nphi,swap=bs),(nphi,ntheta,nr,nq), order = 'F')
            self.vals[:,:,:,:,i] = tmp
            self.time[i] = swapread(fd,dtype='float64',count=1,swap=bs)
            self.iters[i] = swapread(fd,dtype='int32',count=1,swap=bs)
        maxq = 250
        lut = np.zeros(maxq)+int(1000)
        self.lut = lut.astype('int32')
        for i,q in enumerate(self.qv):
            self.lut[q] = i
        fd.close()

class ShellSpectra:
    """Rayleigh Shell Spectra Structure
    ----------------------------------
    self.n_r     : number of radial points
    self.radius  : radial coordinates
    """

    def __init__(self,filename='none',path='Shell_Spectra/',esize=8):
        """filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename
           esize     : The size of each data element in bytes (default 8-byte doubles)"""
        if (filename == 'none'):
            the_file = path+'00000001'
        else:
            the_file = path+filename
        fd = open(the_file,'rb')
        # We read an integer to assess which endian the file was written in...
        bs = check_endian(fd,314,'int32')
        version = swapread(fd,dtype='int32',count=1,swap=bs)
        nrec = swapread(fd,dtype='int32',count=1,swap=bs)
        print bs, version, nrec
        ntheta = swapread(fd,dtype='int32',count=1,swap=bs)
        ntheta = ntheta+1   #quick coding - ntheta is really lmax
        nphi = ntheta
        nr = swapread(fd,dtype='int32',count=1,swap=bs)
        nq = swapread(fd,dtype='int32',count=1,swap=bs)

        self.niter = nrec
        self.nq = nq
        self.nr = nr
        self.ntheta = ntheta
        self.nphi = nphi

        self.qv = np.reshape(swapread(fd,dtype='int32',count=nq,swap=bs),(nq), order = 'F')
        self.radius = np.reshape(swapread(fd,dtype='float64',count=nr,swap=bs),(nr), order = 'F')
        self.inds = np.reshape(swapread(fd,dtype='int32',count=nr,swap=bs),(nr), order = 'F')

        self.vals  = np.zeros((nphi,ntheta,nr,nq,2,nrec),dtype='float64')
        self.iters = np.zeros(nrec,dtype='int32')
        self.time  = np.zeros(nrec,dtype='float64')
        self.version = version
        for i in range(nrec):
            tmp = np.reshape(swapread(fd,dtype='float64',count=nq*nr*ntheta*nphi,swap=bs),(nphi,ntheta,nr,nq), order = 'F')
            self.vals[:,:,:,:,0,i] = tmp
            tmp = np.reshape(swapread(fd,dtype='float64',count=nq*nr*ntheta*nphi,swap=bs),(nphi,ntheta,nr,nq), order = 'F')
            self.vals[:,:,:,:,1,i] = tmp
            self.time[i] = swapread(fd,dtype='float64',count=1,swap=bs)
            self.iters[i] = swapread(fd,dtype='int32',count=1,swap=bs)
        maxq = 250
        lut = np.zeros(maxq)+int(1000)
        self.lut = lut.astype('int32')
        for i,q in enumerate(self.qv):
            self.lut[q] = i
        fd.close()


def swapread(fd,dtype='float64',count=1,swap=False):
        #simple wrapper to numpy.fromfile that allows byteswapping based on Boolean swap
        if (swap):
                val = np.fromfile(fd,dtype=dtype,count=count).byteswap()
        else:
                val = np.fromfile(fd,dtype=dtype,count=count)
        if (len(val) == 1):
                val = val[0]
        return val

def check_endian(fd,sig,sigtype):
    # returns False if first element read from file matches sig
    # True otherwise
    chk = np.fromfile(fd,dtype=sigtype,count=1)
    if (chk == sig):
        return False
    else:
        return True

def build_file_list(istart,iend,path = '.',diter = -1,ndig = 8,special=False):
    files = []
    if (diter < 1):
        # Examine the directory and grab all files that fall between istart and iend
        allfiles = os.listdir(path)
        allfiles.sort()
        for f in allfiles:
            if ( ('special' in f) and special ):
                fint = int(f[0:7])
                if ( (fint >= istart ) and (fint <= iend)  ):
                    files.append(path+'/'+f)
            if ( (not 'special' in f) and not(special) ):
                fint = int(f)
                if ( (fint >= istart ) and (fint <= iend)  ):
                    files.append(path+'/'+f)
    else:
        # Generate filename manually (no ls)
        i = istart
        digmod = "%0"+str(ndig)+"d"
        while (i <= iend):
            fiter = digmod % i           
            if (special):
                fiter=fiter+'_special'
            files.append(path+'/'+fiter)
            i = i+diter
    return files

########################################################
#  These routines allow us to time averages or compile multiple diagnostic files
#  and write out a single file in the same format (for use with viz routines for example).
def Compile_GlobalAverages(file_list,ofile):
    nfiles = len(file_list)
    #   We read the first file, assume that nrec doesn't change
    #   and use the nrecs + nq in the file to create our combined array
    a = GlobalAverage(file_list[0], path = '')
    nfiles = len(file_list)
    niter_estimate = a.niter*nfiles
    nq = a.nq
    combined = np.zeros((niter_estimate,a.nq),dtype='float64')
    time = np.zeros(niter_estimate,dtype='float64')
    iters = np.zeros(niter_estimate,dtype='int32')
    ncount = 0 # total number of iterations read so far
    ind = 0

    # We open the file that we want to store the compiled time traces into and write a header
    fd = open(ofile,'wb') #w = write, b = binary
    dims = np.zeros(4,dtype='int32')
    dims[0] = 314
    dims[1] = a.version
    dims[2] = a.niter   # We will fix this at the end
    dims[3] = a.nq
    dims.tofile(fd)
    a.qv.tofile(fd)

    tmp = np.zeros(a.nq,dtype='float64')
    simtime   = np.zeros(1,dtype='float64')
    iteration = np.zeros(1,dtype='int32')
    icount = np.zeros(1,dtype='int32')
    icount[0] = 0
    for i in range(nfiles):
        the_file = file_list[i]
        print the_file
        a = GlobalAverage(the_file,path='')
        nrec = a.niter
        for j in range(nrec):
            tmp[:] = a.vals[j,:]
            tmp.tofile(fd)
            iteration = a.iters[j]
            simtime = a.time[j]
            simtime.tofile(fd)
            iteration.tofile(fd)
            icount[0] = icount[0]+1
    fd.seek(8)
    icount.tofile(fd)   # insert the proper number of iterations
    fd.close()

def TimeAvg_AZAverages(file_list,ofile):
    nfiles = len(file_list)
    #   We read the first file, assume that nrec doesn't change
    #   and use the nrecs + nq in the file to create our combined array
    a = AzAverage(file_list[0], path = '')
    nfiles = len(file_list)


    nr = a.nr
    ntheta = a.ntheta
    nq = a.nq
    tmp = np.zeros((nr,ntheta,nq),dtype='float64')
    simtime   = np.zeros(1,dtype='float64')
    iteration = np.zeros(1,dtype='int32')
    icount = np.zeros(1,dtype='int32')
    ifinal = np.zeros(1,dtype='int32')
    tfinal = np.zeros(1,dtype='float64')
    icount[0] = 0
    i0 = a.iters[0]
    t0 = a.time[0]
    for i in range(0,nfiles):
        the_file = file_list[i]
        print the_file
        b = AzAverage(the_file,path='')
        nrec = b.niter
        for j in range(nrec):
            tmp[0:nr,0:ntheta,0:nq] += b.vals[0:nr,0:ntheta,0:nq,j].astype('float64')

            tfinal[0] = b.time[j]
            ifinal[0] = b.iters[j]
            icount[0] = icount[0]+1
    div = np.float(icount[0])
    print div
    tmp = tmp/div

    # We open the file that we want to store the compiled time traces into and write a header
    fd = open(ofile,'wb') #w = write, b = binary
    dims = np.zeros(6,dtype='int32')
    dims[0] = 314
    dims[1] = a.version
    dims[2] = 1
    dims[3] = a.nr
    dims[4] = a.ntheta
    dims[5] = a.nq
    dims.tofile(fd)
    a.qv.tofile(fd)
    a.radius.tofile(fd)
    a.costheta.tofile(fd)


    test = np.transpose(tmp)
    test.tofile(fd)
    t0.tofile(fd)
    i0.tofile(fd)
    # The final structure is identical to a normal az_average file save for the fact that final iteration adn final time are saved
    tfinal.tofile(fd)
    ifinal.tofile(fd)
    fd.close()
