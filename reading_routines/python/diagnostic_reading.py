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

class ShellAverage:
    """Rayleigh ShellAverage Structure
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
        
        nr = swapread(fd,dtype='int32',count=1,swap=bs)
        nq = swapread(fd,dtype='int32',count=1,swap=bs)
        print nr,nq,bs
        qv = np.reshape(swapread(fd,dtype='int32',count=nq,swap=bs),(nq), order = 'F')
        rad = np.reshape(swapread(fd,dtype='float64',count=nr,swap=bs),(nr), order = 'F')
        tmp = np.reshape(swapread(fd,dtype='float64',count=nq*nr,swap=bs),(nr,nq), order = 'F')
        self.nr = nr
        self.nq = nq
        self.radius      = rad
        self.vals = tmp
        self.qv = qv
        fd.close()


class ShellSlice:
    """Rayleigh ShellAverage Structure
    ----------------------------------
    self.n_r     : number of radial points
    self.radius  : radial coordinates
    """

    def __init__(self,filename='none',path='Shell_Slices/',esize=8):
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
        
        nphi = swapread(fd,dtype='int32',count=1,swap=bs)
        ntheta = swapread(fd,dtype='int32',count=1,swap=bs)
        nr = swapread(fd,dtype='int32',count=1,swap=bs)
        nq = swapread(fd,dtype='int32',count=1,swap=bs)
        qv = np.reshape(swapread(fd,dtype='int32',count=nq,swap=bs),(nq), order = 'F')
        rad = np.reshape(swapread(fd,dtype='float64',count=nr,swap=bs),(nr), order = 'F')
        inds = np.reshape(swapread(fd,dtype='int32',count=nr,swap=bs),(nr), order = 'F')
        stheta = np.reshape(swapread(fd,dtype='float64',count=ntheta,swap=bs),(ntheta), order = 'F')
        tmp = np.reshape(swapread(fd,dtype='float64',count=nq*nr*ntheta*nphi,swap=bs),(nphi,ntheta,nr,nq), order = 'F')
        ctheta = np.reshape(swapread(fd,dtype='float64',count=ntheta,swap=bs),(ntheta), order = 'F')
        simtime = swapread(fd,dtype='float64',count=1,swap=bs)
        self.nr = nr
        self.nq = nq
        self.qv = qv
        self.radius      = rad
        self.inds = inds
        self.vals = tmp
        self.ntheta = ntheta
        self.nphi = nphi
        self.stheta = stheta
        self.ctheta = ctheta
        self.simtime = simtime
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

