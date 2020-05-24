import numpy as np
import healpy as hp
from reproject import reproject_from_healpix

if __name__ == "__main__":
        
    wisenside = 8192
    wiseroot = "/data/seclark/WISE12micron/"
    wisehpfn = wiseroot+"wssa_sample_{}.fits".format(wisenside)
    wisehpdata = fits.getdata(wisehpfn, hdu=0)
    
    print(wisehpdata.shape)