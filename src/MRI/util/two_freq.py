import numpy as np
from scipy import interpolate

import matplotlib.pyplot as plt


def order_umt_cest( cest, frqs ):
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

	"""
	plt.subplot(211)
	plt.plot( frqs, cest, 'ko')
	plt.plot( resample_freqs, hr,'b')
	
	plt.axvline( resample_freqs[rmax], color='r' )
	plt.axvline( resample_freqs[lmax], color='r' )
	plt.axvline( resample_freqs[mid], color='r' )

	plt.subplot(212)
	frq2 = np.hstack( ( resample_freqs[lmax:mid], resample_freqs[mid:rmax] ) )
	spec = np.hstack( (  hr[mid:rmax], hr[lmax:mid] ) )
	
	plt.plot( (frqs[ valids ]-offset)/300, cest_umt, 'o' )
	plt.plot( (frq2-offset)/300, spec )

		
	plt.show()	
	"""
	return frqs[ valids ]-offset, cest_umt



if __name__ == '__main__':
	import matplotlib.pyplot as plt
	import numpy as np

	frqs= np.linspace(-2500,2500,51)
	cest = np.load('/Users/josh/Documents/Data/Tumor_dataset/uMT_cest.p')

	#plt.plot( frqs, cest[0], 'o' )

	hr = order_umt_cest(cest[0], frqs) 
	plt.plot( hr[0], hr[1]); plt.show()
