import numpy as np
from astropy.io import fits
from scipy import ndimage
import h5py

def gaussian_umask(data, fwhm=10, zeroed=True):
    """
    fwhm in arcmin aka pixels
    """
    sigma = fwhm/(2*np.sqrt(2*np.log(2)))
    smoothdata = ndimage.filters.gaussian_filter(data, sigma=sigma)
    
    umask = data - smoothdata
    if zeroed:
        umask[np.where(umask < 0)] = 0
    return umask

def stack_slicedata(stackthese_data, stackon_data, cubenx=11, cubeny=11):
    """
    stack data
    """

    # shape of data to be stacked on
    maxny, maxnx = stackon_data.shape
    
    # square to be stacked into
    stackslice = np.zeros((cubeny, cubenx), np.float_)
    weightslice = np.zeros((cubeny, cubenx), np.float_)
    cubehalfx = np.int(np.floor(cubenx/2.0))
    cubehalfy = np.int(np.floor(cubeny/2.0))
    
    # Define stackable region -- just skip the edges
    stackable = np.copy(stackon_data) #np.ones(stackon_data.shape)
    row, col = np.indices(stackable.shape)
    stackable[np.where(row < cubehalfy)] = 0
    stackable[np.where(row > maxny - cubehalfy - 1)] = 0
    stackable[np.where(col < cubehalfx)] = 0
    stackable[np.where(col > maxnx - cubehalfx - 1)] = 0
    
    nonzeroyx = np.nonzero(stackable)
    nonzeroy = nonzeroyx[0]
    nonzerox = nonzeroyx[1]
    
    smallstarty = 0
    smallstopy = cubehalfy*2 + 1
    
    # Step through nonzero pixels
    for _ixy, (_y, _x) in enumerate(zip(nonzeroy, nonzerox)):
        
        # only stack if you're not on the edge
        stackyes = stackable[_y, _x]
        
        if stackyes:

            centerval = stackon_data[_y, _x]

            startx = np.int(max(0, _x - cubehalfx))
            stopx = np.int(min(maxnx, _x + cubehalfx + 1))
            starty = np.int(max(0, _y - cubehalfy))
            stopy = np.int(min(maxny, _y + cubehalfy + 1))
            try:
                stackslice[smallstarty:smallstopy, :] += centerval * stackthese_data[starty:stopy, startx:stopx]    
                weightslice[smallstarty:smallstopy, :] += centerval
            except:
                print(startx, stopx, starty, stopy, "small", smallstarty, smallstopy, stackslice[smallstarty:smallstopy, :].shape, stackthese_data[starty:stopy, startx:stopx].shape)

        else:
            pass
            
    return stackslice, weightslice

if __name__ == "__main__":
    
    # Get data
    dataroot = "../data/"
    WISE12 = fits.getdata(dataroot+"WISE_8192_on_ACT.fits")
    
    ACTfn = dataroot+"act_planck_s08_s18_cmb_f220_daynight_map_feb.fits"
    ACT220 = fits.getdata(ACTfn)
    act220hdr = fits.getheader(ACTfn)
    actpixelsize_arcmin = act220hdr["CDELT2"]*60 # pixel size in arcminutes
    
    ACTI = ACT220[0, :, :]
    ACTQ = ACT220[1, :, :]
    ACTU = ACT220[2, :, :]
    
    xstart = 8000
    xstop = 20000
    ystart = 100
    ystop = 8000

    WISEsub = WISE12[ystart:ystop, xstart:xstop]
    ACTIsub = ACTI[ystart:ystop, xstart:xstop]
    ACTQsub = ACTQ[ystart:ystop, xstart:xstop]
    ACTUsub = ACTU[ystart:ystop, xstart:xstop]
    ACTPsub = np.sqrt(ACTQsub**2 + ACTUsub**2)
    
    # Stack on unsharp masked WISE data
    WISEsub_umask30 = gaussian_umask(WISEsub, fwhm=30/actpixelsize_arcmin, zeroed=True)
    
    # Zero out anomalously high values
    clipdata = np.percentile(WISEsub_umask30[WISEsub_umask30 > 0], 99)
    WISEsub_umask30[WISEsub_umask30 > clipdata] = 0.
    
    size = 101
    stackI, weightI = stack_slicedata(ACTIsub, WISEsub_umask30, cubenx=size, cubeny=size)
    stackQ, weightQ = stack_slicedata(ACTQsub, WISEsub_umask30, cubenx=size, cubeny=size)
    stackU, weightU = stack_slicedata(ACTUsub, WISEsub_umask30, cubenx=size, cubeny=size)
    stackP, weightP = stack_slicedata(ACTPsub, WISEsub_umask30, cubenx=size, cubeny=size)
    
    outfn = dataroot + "stacked_ACTonWISEumask_x{}_{}_y{}_{}_size{}_WISEclip.h5".format(xstart, xstop, ystart, ystop, size)
    with h5py.File(outfn, 'w') as f:
        Istack = f.create_dataset(name='stackI', data=stackI)
        Qstack = f.create_dataset(name='stackQ', data=stackQ)
        Ustack = f.create_dataset(name='stackU', data=stackU)
        Pstack = f.create_dataset(name='stackP', data=stackP)
        Iweight = f.create_dataset(name='weightI', data=weightI)
        Qweight = f.create_dataset(name='weightQ', data=weightQ)
        Uweight = f.create_dataset(name='weightU', data=weightU)
        Pweight = f.create_dataset(name='weightP', data=weightP)
        
    
    
