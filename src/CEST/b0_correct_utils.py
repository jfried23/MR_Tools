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
	if len(wsr_offs) != wassr_img.shape[0]:
		raise IndexError("Input frequencies offsets do not match size of wassr_img!") 

	if min(wsr_offs) != wsr_offs[0]: wsr_offs = wsr_offs[::-1]

	resample_freqs = np.arange( min(wsr_offs), max(wsr_offs), 1 )
	

	f = interpolate.interp1d( wsr_offs, wassr_img, kind='cubic', axis=0, assume_sorted=True)

	re_index = np.argmin( f( resample_freqs  ), axis=0 )
	
	b0_map = np.array([resample_freqs[i] for i in re_index.ravel()])
	
	return b0_map.reshape( wassr_img.shape[1:] )


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
	

def apply_B0_corr( cest_img, b0_map, cest_offs ):
	"""
	Apply B0 correction to a CEST dataset using a B0 map generated from WASSR dataset

	Parameters
	----------
	cest_img  : 3D numpy array
			The CEST dataset in shape [ #offsets, phase_enc, freq_encode]  
	b0_map    : 2D numpy array
			The estimated B0 map in unit Hz. See 'calc_B0_map' for details.
	cest_offs : 1D numpy array
			The frequencey offset of each plane in the CEST dataset (units Hz).

	Returns
	-------
	cest_b0_corr: 3D numpy array
			The B0 corrected CEST dataset.
	"""
	if len(cest_offs) != cest_img.shape[0]:
		raise IndexError("Input frequencies offsets do not match size of cest_img!") 

	cest_b0_corr = np.empty_like( cest_img )

	for x in range(cest_img.shape[-1]):
    		for y in range( cest_img.shape[-2]):
        		f = interpolate.interp1d( cest_offs-b0_map[y,x], cest_img[:,y,x], kind='cubic', bounds_error = False )
       			cest_b0_corr[:,y,x] = f(cest_offs) 
		
	return cest_b0_corr	



if __name__ == '__main__':
	import pickle
	import matplotlib.pyplot as plt

	test_b0    = 0
	test_apply = 1	

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

		
