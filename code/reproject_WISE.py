import numpy as np
import healpy as hp
from astropy.io import fits
from reproject import reproject_from_healpix
import time

if __name__ == "__main__":
        
    wisenside = 8192
    wiseroot = "/data/seclark/WISE12micron/"
    wisehpfn = wiseroot+"wssa_sample_{}.fits".format(wisenside)
    wisehphdr = fits.getheader(wisehpfn)
    wisehpdata = fits.getdata(wisehpfn, hdu=0)
    
    print(wisehpdata.shape)
    
    act220fn = "../data/act_planck_s08_s18_cmb_f220_daynight_map_feb.fits"
    act220hdr = fits.getheader(act220fn)
    
    act220hdr2D = fits.Header()
    act220hdr2D["SIMPLE"] = act220hdr["SIMPLE"]
    act220hdr2D["BITPIX"] = act220hdr["BITPIX"]
    act220hdr2D["NAXIS"] = 2
    act220hdr2D["NAXIS1"] = act220hdr["NAXIS1"]
    act220hdr2D["NAXIS2"] = act220hdr["NAXIS2"]
    act220hdr2D["WCSAXES"] = act220hdr["WCSAXES"]
    act220hdr2D["CRPIX1"] = act220hdr["CRPIX1"]
    act220hdr2D["CRPIX2"] = act220hdr["CRPIX2"]
    act220hdr2D["CDELT1"] = act220hdr["CDELT1"]
    act220hdr2D["CDELT2"] = act220hdr["CDELT2"]
    act220hdr2D["CUNIT1"] = act220hdr["CUNIT1"]
    act220hdr2D["CUNIT2"] = act220hdr["CUNIT2"]
    act220hdr2D["CTYPE1"] = act220hdr["CTYPE1"]
    act220hdr2D["CTYPE2"] = act220hdr["CTYPE2"]
    act220hdr2D["CRVAL1"] = act220hdr["CRVAL1"]
    act220hdr2D["CRVAL2"] = act220hdr["CRVAL2"]
    act220hdr2D["LONPOLE"] = act220hdr["LONPOLE"]
    act220hdr2D["LATPOLE"] = act220hdr["LATPOLE"]
    act220hdr2D["RADESYS"] = act220hdr["RADESYS"]
    
    print(wisehphdr)
    print(act220hdr2Dre)
    
    time0 = time.time()
    reproj_wise, footprint = reproject_from_healpix((wisehpdata, 'galactic'), act220hdr2D, nested=False)
    time1 = time.time()
    
    print("time = {} minutes".format((time1 - time0)/60.))
    
    fits.writeto("../data/WISE_8192_on_ACT.fits", reproj_wise, header=act220hdr2D)
    
    