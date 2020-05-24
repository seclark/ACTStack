import numpy as np
import healpy as hp
from astropy.io import fits
from reproject import reproject_from_healpix

if __name__ == "__main__":
        
    wisenside = 8192
    wiseroot = "/data/seclark/WISE12micron/"
    wisehpfn = wiseroot+"wssa_sample_{}.fits".format(wisenside)
    wisehpdata = fits.getdata(wisehpfn, hdu=0)
    
    print(wisehpdata.shape)
    
    reproj_wise, footprint = reproject_from_healpix((wisehpdata, wisehphdr), act220hdr2D, nested=False)