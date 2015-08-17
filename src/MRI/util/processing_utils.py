import numpy as np
from skimage import filters

def process_sum_square( fids, ndim=2  ):
	"""
	Utility function for process images from MRI fids and combining the channels
	using the sum of squares processing; reducing the diemsionality of the data by 1.

	Parameters
	----------
	fids:	   Numpy array containing the complex fids
                   Array must be ordered as [ num_channels, phe2, phe1, npts]

	ndim:      The dimensionality of the dataset.
		   (i.e. ndim = 2, means FT last two dimensions, ndim=3 tansform the last 3)
		  

	Returns
	-------
	
	A real ndarray encoding the images

	"""
	
	img = process( fids, ndim)
	img *= np.conj(img) 
	
	return np.squeeze(np.sqrt( np.sum( np.abs(img), axis=0) ))	


def  process( fids, ndim=2  ):
	"""
	Utility function for process images from MRI fids.

	Parameters
	----------
	fids:	   Numpy array containing the complex fids
                   Array must be ordered as [ num_channels, phe2, phe1, npts]

	ndim:      The dimensionality of the dataset.
		   (i.e. ndim = 2, means FT last two dimensions, ndim=3 tansform the last 3)
		  

	Returns
	-------
	
	A complex ndarray encoding the images

	"""
	
	img = np.empty_like( fids )
	ax = -1*(np.array( range(ndim) )+1)
	img = np.fft.fftshift( np.fft.ifftn( fids, axes=ax ), axes=ax )

	return np.squeeze(img)

def mask_background( img, mask=None, return_mask = False ):
	"""
	Utility function for masking MRI images to remove the background.

	Parameters
	----------
	img:	     Numpy array containing the CEST data images
                     Array must be ordered as [ phe2, phe1, npts]

	mask:       opitional mask matrix precuomputed elsewhere. 
		    Must be same size as a single image.
		    If not provided the skimage filter 'threshold_li' will be used
		    To compute the mask. 

	return_mask:  Bool defualt False. Return the mask or the masked image
		      If True this function returns the binary mask indicating
		      which pixels are part of the image.
                      If False this function returns the masked image.   

	Returns
	-------
	
	A masked MRI image or the image mask.

	"""
	if mask == None:
		
		if   len( img.shape ) == 3: t = img[0]
		elif len( img.shape ) == 2: t = img

		mask = img > filters.threshold_li(t)

	if return_mask: return mask

	else:  return np.ma.MaskedArray( np.multiply( img, mask ) )		
