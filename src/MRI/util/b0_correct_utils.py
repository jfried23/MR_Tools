import numpy as np
from scipy import interpolate

def calc_B0_map( wassr_img, wsr_offs, resolution = 1.0 ):
	"""
	Return the estimated the B0 map (in units Hz) from a WASSR dataset.

	Parameters
	----------
	wassr_img : 3D numpy array
			The WASSR dataset in shape [ #offsets, phase_enc, freq_encode]  
	wsr_offs  : 1D numpy array
			The frequencey offset of each plane in the WASSR dataset in units Hz.
	resolution: float, optional
			The frequencey resolution the of the resultant B0 map.


	Returns
	-------
	b0_map: 2D numpy array	

	"""
	if wsr_offs.shape[0] != wassr_img.shape[0]:
		raise IndexError("Input frequencies offsets do not match size of wassr_img!") 

	if min(wsr_offs) != wsr_offs[0]: wsr_offs = wsr_offs[::-1]

	resample_freqs = np.arange( np.min(wsr_offs), np.max(wsr_offs), resolution, dtype=np.float )
	

	f = interpolate.interp1d( wsr_offs, wassr_img, kind='cubic', axis=0)

	re_index = np.argmin( f( resample_freqs  ), axis=0 )
	
	b0_map = np.array([resample_freqs[i] for i in re_index.ravel()])
	
	return b0_map.reshape( wassr_img.shape[1:] )


def calc_B0_map_umtCEST(  umt_imgs, frq_off, resolution = 1. ):
	"""
	Estimate the B0 map from a uMT CEST dataset 

	Parameters
	----------
	umt_imgs : 3D numpy array
			The CEST dataset  
	frq_off  : 1D numpy array
			The stauration freqencey of each point in dimension 0 (in Hz)
	resolution: Resampling resolution (Hz)

	Returns
	-------
	A B0 map of the 

	"""
	delta = frq_off[1] - frq_off[0]
	
	resample_freqs = np.arange( np.min(frq_off), np.max(frq_off), resolution, dtype=np.float )

	f = interpolate.interp1d( frq_off, umt_imgs, kind='cubic', axis=0, bounds_error=False)
	hr = f( resample_freqs  )
	
	lh, rh =  hr[0:hr.shape[0]/2], hr[hr.shape[0]/2:]


	lmax, rmax = np.argmin(lh, axis=0), np.argmin(rh, axis=0) + len(resample_freqs)/2
	mid = ((rmax-lmax)/2) + lmax

	b0map = np.zeros_like( umt_imgs[0] )

	for i in range( mid.shape[-2] ):
		for ii in range( mid.shape[-1]):

			if umt_imgs[0].mask[i,ii]: continue
			else: b0map[i,ii] = hr[ mid[i,ii], i, ii ]


	return np.ma.masked_array( b0map, mask = umt_imgs.mask[0] )	
	

def correct_umt_ret_b0map( imgs, frqs, resolution =1.0 ):
	"""
	B0 correct a uMT-CEST dataset and return the estimated b0 map

	Parameters
	----------
	imgs : 1D numpy array
			The CEST spectra  
	frqs  : 1D numpy array
			The frequencey offset of each point

	Returns
	-------
	A B0 map

	"""

	resample_freqs = np.arange( min(frqs), max(frqs), resolution, dtype=np.float )
	out = np.zeros_like( imgs )

	f = interpolate.interp1d( frqs, imgs, kind='cubic', axis=0, fill_value=0.0)
	hr = f( resample_freqs  )
	
	lh, rh =  hr[0:hr.shape[0]/2], hr[hr.shape[0]/2:]

	lmax       = np.argmin(lh, axis=0)
	rmax	   = ((np.argmin(rh, axis=0) + len(resample_freqs)/2))

	mid = np.ma.masked_array( ((rmax-lmax)/2) + lmax, imgs.mask[0] )
	
	b0_map = np.zeros_like( imgs[0] )

	for i in range( imgs.shape[-2] ):
    		for ii in range( imgs.shape[-1] ):
        		if mid.mask[ i, ii]: continue
        
        		offset = resample_freqs[mid[i,ii]]
			b0_map[i,ii] = offset       			
 
        		left  = frqs > resample_freqs[lmax[i,ii]]
        		right = frqs < resample_freqs[rmax[i,ii]]
        
        		valids = left*right    
        
			#Now resample the points
			floc = interpolate.interp1d( frqs, imgs[:,i,ii], 
                                    kind='cubic', fill_value='NaN', bounds_error=False)

			imgs[:,i,ii] = floc( frqs - offset )

	return b0_map





def apply_B0_corr( imgs, b0_map, cest_offs ):
	"""
	Apply B0 correction to a CEST dataset using a B0 map

	Parameters
	----------
	imgs         : 3D numpy array
			The CEST dataset in shape [ #offsets, phase_enc, freq_encode]  
	b0_map    : 2D numpy array
			The estimated B0 map in unit Hz. See 'calc_B0_map' for details.
	cest_offs : 1D numpy array
			The frequencey offset of each plane in the CEST dataset (units Hz).

	Returns
	-------
	"""
	if len(cest_offs) != imgs.shape[0]:
		raise IndexError("Input frequencies offsets do not match size of cest_img!") 

	for i in range( b0_map.shape[-2] ):
		for ii in range( b0_map.shape[-1] ):
			if b0_map.mask[ i, ii]: continue
			f = interpolate.interp1d( cest_offs, imgs[:,i,ii], 
					kind='cubic', bounds_error=False,  fill_value='NaN')
			imgs[:,i,ii] =  f( cest_offs - b0_map[i,ii]  )


def calib_off_freq( b0_map, cest_offs ):
	"""
	Utility function.
	Generate the calibrated saturation offset frequencey of each pixel in MRI image

	Parameters
	----------
	b0_map    : 2D numpy array
			The estimated B0 map in unit Hz. See 'calc_B0_map' for details.
	cest_offs : 1D numpy array
			The frequencey offset of each plane in the CEST dataset (units Hz).

	Returns
	-------
	calib_freq: 3D numpy array
			The corrected offset frequencies for each point in the MRI
			[ cest_feq_dim, phase_encode, freq_encode]
	""" 

	ypts, xpts = b0map.shape
	calib_freq = np.empty( [ len(cest_offs), ypts, xpts ] )	
	for x in range( xpts ):
		for y in range( ypts ):
			calib_freq[:, y, x] = cest_offs - b0_map[y,x]

	return calib_freq


if __name__ == '__main__':
	import pickle
	import matplotlib.pyplot as plt

	test_b0    = 1
	test_apply = 0	

	if test_b0:
		wsr_ofs = np.arange (-500,520,20)
		wsr_img = pickle.load( open('/Users/josh/Desktop/wassr.data','rb') )
	
		b0_map = calc_B0_map( wsr_img, wsr_ofs )
		pickle.dump( b0_map, open('/Users/josh/Desktop/b0map.data','wb') )	
		
		

	if test_apply:
		b0map = pickle.load( open('/Users/josh/Desktop/b0map.data','rb') )
		cest = pickle.load( open('/Users/josh/Desktop/cest.data','rb') )

		cst_ofs = np.arange(-2500,2600,100)

		
		cest_b0_corr = apply_B0_corr( cest, b0map, cst_ofs )

		
