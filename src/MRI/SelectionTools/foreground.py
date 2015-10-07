import numpy as np
import matplotlib.pyplot as plt

from skimage import feature, segmentation
from scipy import ndimage as ndi


class foreground( object ):

	def __init__(self): pass

	def getMask(self, currentImage, fac=1.0):
		if   len( currentImage.shape ) == 3: t = currentImage[0]
		elif len( currentImage.shape ) == 2: t = currentImage

		mask = currentImage > filters.threshold_li(t)*fac

		return mask



class skull_strip( object ):
	def __init__(self): pass

	def getMask( self, img. sigma=0.8 ):
		f = foreground()
		edges = feature.canny(f.getMask(img), sigma=sigma)	
		return ndi.binary_fill_holes(edges)
		

if __name__ == '__main__':

	import h5py


	f = h5py.File('/Users/josh/Documents/Data/Tumor_Sep15/20151005_MRIs.hdf5')
	img = f['meas_MID59_Gre_2D_Gauss_10X1440d_200FOV_2500Hz_FID3330'][...]
	
	f.close()

	f = h5py.File('/Users/josh/Documents/Data/Tumor_Oct10/20151005_MRIs.hdf5')
	img = f['meas_MID99_Ref_Gre_2D_12ms_9d_1shot_FID6464'][...]
	
	f.close()
	

	#f = h5py.File('/Users/josh/Documents/Data/Tumor_dataset/20150930_MRIs.hdf5')
	#img = f['meas_MID163_Ref_Gre_2D_12ms_9d_1shot_9_3s_FID4493.dat'][...]	



	plt.imshow( img, cmap=plt.cm.binary_r )

	"""
	ROI1 = foreground() #let user draw first ROI

	# show the image with the first ROI
	#plt.imshow(img[0], interpolation='nearest', cmap="Greys")

	msk = ROI1.getMask( img[0] ) 


	#plt.imshow( img[0],cmap=plt.cm.binary_r  )	
	plt.imshow( np.ma.masked_equal( (msk*img[0]), 0.0), alpha=1.0 )

	plt.show()
	#np.ma.make_mask( msk )

	"""

	from skimage import feature, segmentation
	from scipy import ndimage as ndi

	ROI1 = foreground()
	
	edges2 = feature.canny(ROI1.getMask(img), sigma=0.5)	
	#edges2 = segmentation.find_boundaries( ROI1.getMask(img) )


	fill_brain = ndi.binary_fill_holes(edges2)

	plt.imshow( img, cmap=plt.cm.binary) 
	plt.imshow( np.ma.masked_equal( fill_brain*img, 0) )
	plt.show()
