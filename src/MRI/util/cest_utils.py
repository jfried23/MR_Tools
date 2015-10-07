import numpy as np
from scipy import interpolate

def eval_cest_roi( imgs, roi_mask ):
	"""
	Generate the CEST spectra over a ROI defined in a masked array.

	Parameters
	----------
	imgs :      3D numpy array
			The CEST like dataset  
	roi_mask  : 2D numpy boolean masked array
			Include point in the calculation? T/F

	Returns
	-------
	CEST spectra for the ROI

	"""
	if imgs.shape[-2:] != roi_mask.shape: raise ValueError
	return 	(imgs*roi_mask).mean( axis=-2).mean( axis=-1 )


def umt_cest_process( cest, frqs ):
	"""
	Reorder a uMT-CEST spectra

	Parameters
	----------
	cest : 1D numpy array
			The CEST spectra  
	frqs  : 1D numpy array
			The frequencey offset of each point

	Returns
	-------
	corrected CEST data

	"""
	resample_freqs = np.linspace( min( frqs), max(frqs), 1001 )	

	f = interpolate.interp1d( frqs, cest, kind='cubic', axis=0)
	hr = f( resample_freqs  )
	
	lh, rh =  hr[0:hr.size/2], hr[hr.size/2:]
	lmax, rmax = np.argmin(lh), np.argmin(rh)+500
	mid = ((rmax-lmax)/2) + lmax
	
	offset = resample_freqs[ mid ]
	
	#Now find the points in between these two limits
	left  = frqs > resample_freqs[lmax]
	right = frqs < resample_freqs[rmax]

	valids = np.multiply( left, right )

	cest_umt = cest[ valids ]
	freq_umt = frqs[ valids ]

	

	mdpt = cest_umt.size/2

	cest_umt = np.hstack( ( cest_umt[ mdpt: ], cest_umt[ 0:mdpt ] ) )

	return frqs[ valids ]-offset, cest_umt
		
