import numpy as np
from skimage import filters

def process_sum_square( fids, ndim=2  ):
	"""
	Utility function for process images from MRI fids and combining the channels
	using the sum of squares processing; reducing the diemsionality of the data by 1.

	Applies ndim FFT and zero freq shift prior to combining channel data by sum-of-sqaures method.

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
        Applies ndim FFT and zero freq shift.

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
	
	img = np.fft.fftshift( np.fft.fftn( fids, axes=ax, ).astype( np.complex64), axes=ax )
	
	return np.squeeze(img)

def gen_background_mask( img ):
	"""
	Utility function for making MRI image masks to remove the background via threshold filter.

	Parameters
	----------
	img:	     Numpy array containing the CEST data images
                     Array must be ordered as [ phe2, phe1, npts]

	Returns
	-------
	
	A 2D image mask.

	"""
		
	if   len( img.shape ) == 3: t = img[0]
	elif len( img.shape ) == 2: t = img

	mask = img > filters.threshold_li(t)

	return mask

def apply_background_mask( mask, img ):
	return np.ma.MaskedArray( np.multiply( img, mask ) )		
