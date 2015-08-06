

import numpy as np

def Process_MRI( fids, ndim=2, sum_sqare=0 ):
	ax = -1*(np.array( range(ndim) )+1)
	
	img = np.fft.fftshift( np.fft.ifftn( fids, axes=ax ), axes=ax )
	img *= np.conj(img) 
	
	return np.sqrt( np.sum( np.abs(img), axis=sum_sqare) )	
