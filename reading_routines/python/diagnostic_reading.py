import numpy as np

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
        print nrec, nq
        self.niter = nrec
        self.nq = nq
        self.qv = np.reshape(swapread(fd,dtype='int32',count=nq,swap=bs),(nq), order = 'F')
        self.vals  = np.zeros((nq,nrec),dtype='float64')
        self.iters = np.zeros(nrec,dtype='int32')
        self.time  = np.zeros(nrec,dtype='float64')
        self.version = version
        for i in range(nrec):
            tmp = np.reshape(swapread(fd,dtype='float64',count=nq,swap=bs),(nq), order = 'F')
            self.vals[:,i] = tmp
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
    print chk
    if (chk == sig):
        return False
    else:
        return True

