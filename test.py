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

# this scirpt gets all the files in a directoty and appends them to a list to be used stating from all_files[0] as 1st file

image_list = []
for frame in glob.glob("*.fit"):
	image_list.append(frame)

 #next we want to go through the list and align them
#
#ref_image = image_list[0]
#images_to_align = image_list[1:len(image_list)]

#identifications = make_transforms(ref_image, images_to_align)

#aligned_images = [ref_image]
#for id in identifications:
#	if id.ok:
#		alignedimg = affineremap(id.ukn.filepath, id.trans, outdir=tmpdir)
#aligned_images.append(alignedimg)

# finally we want to subtract all the files in the list 

#option 2 original data reworked for alignment of a list of images
aligned_list = []
counter = 1
a =len(image_list)
for image in image_list:
	#if counter > a
	#	exit
	#else

		ref_fits = fits.open(image_list[0])
		targ_fits = fits.open(image_list[counter])

# Read FITS data
		ref_data = fits.getdata(image_list[0])
		targ_data = fits.getdata(image_list[counter])

# Read FITS headers
		ref_hdr = ref_fits[0].header
		targ_hdr = targ_fits[0].header
		aligned_array = []

		#aligned_array, footprint = reproject_interp(targ_fits, ref_hdr) #aligned array is the new fits file aligned
		#print(type(aligned_array))
		aligned_hdu = fits.PrimaryHDU(aligned_array, header=ref_hdr)
		aligned_im = aligned_hdu.writeto("aligned.fit", overwrite=True)
		aligned_list.append(aligned_im)
		counter =+ 1


sub_images = []
nu_counter = 0

for image in aligned_list:
	ref_fits = fits.open(aligned_list[0])
	targ_fits = fits.open(aligned_list[counter])

# Read FITS data
	ref_data = fits.getdata(aligned_list[0])
	targ_data = fits.getdata(aligned_list[counter])

# Read FITS headers
	ref_hdr = ref_fits[0].header
	targ_hdr = targ_fits[0].header
	
	diff_data, optimal_image, kernel, background = ois.optimal_system(np.array(image_list[nu_counter]), np.array(ref_data), kernelshape = (1,1), method = "Bramich") 
	hdu_diff = fits.PrimaryHDU(diff_data, header=ref_hdr)
	sub_im = hdu_diff.writeto("sub.fit", overwrite = True)
	sub_images.append(sub_im)
	nu_counter =+ 1 

#stacking the final images to get our toal fits differential subtraction
final_image = np.sum(sub_images, axis=0)

time = end - start
print("Ended in", "%.1f" % time, "seconds!")