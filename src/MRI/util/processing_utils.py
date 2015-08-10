import numpy as np

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
	
	return np.sqrt( np.sum( np.abs(img), axis=0) )	


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

	return img
