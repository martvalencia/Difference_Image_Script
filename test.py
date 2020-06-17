from astropy.io import fits
from reproject import reproject_interp
from fits_align.ident import make_transforms
from fits_align.align import affineremap
from glob import glob
from numpy import shape
import glob
import ois
import numpy as np
import time
import os
start = time.time()
#
# this scirpt gets all the files in a directoty and appends them to a list to be used stating from all_files[0] as 1st file
#
image_list_prelim = []
for frame in glob.glob("*.fit"):
	image_list_prelim.append(frame)

image_list = [ fits.getdata(image) for image in image_list_prelim ]

compressed_image = np.sum(image_list, axis=0)
hdu_diff = fits.PrimaryHDU(diff_data, header = fits.getdata(image_list[0]))
# Read FITS data
ref_data = fits.getdata(image_list[0])
targ_data = fits.getdata(compressed_image)
#
# Read FITS headers
ref_hdr = ref_fits[0].header
targ_hdr = targ_fits[0].header
#		aligned_array = []
#
#subtaction of compressed image to reference image
#
diff_data, optimal_image, kernel, background = ois.optimal_system(image_list, ref_data, kernelshape = (11,11), method = "Bramich") 
hdu_diff = fits.PrimaryHDU(diff_data, header=ref_hdr)
hdu_diff.writeto("final_sub.fit", overwrite=True)
#
time = end - start
print("Ended in", "%.1f" % time, "seconds!")
