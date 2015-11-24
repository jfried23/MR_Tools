from numpy import shape, fft, zeros, cov, dot, reshape, linalg, sum, transpose, squeeze
from matplotlib.pyplot import imshow, plot, show, subplot

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from DiscreteSlider import DiscreteSlider

import sys
sys.path.insert(0, '../FileIO/')


class MRI( object ):
	def __init__( self, images ):

		#last two dimensions must be of imaging plane
		self._images = images
		self._shp = shape( self._images )
		self._current_slice = (len(self._shp)-2)*[0] #this images slice

		self._fig = plt.figure()
		self._display_image()

		
		axcolor = 'grey' #'lightgoldenrodyellow'	
		self._sliders=[]
		for d in range(len( self._current_slice )):
			if self._shp[d] == 1: continue				
			tmp = plt.axes([0.10, 0.07+(d*-0.03), 0.75, 0.02], axisbg=axcolor)
			self._sliders.append( DiscreteSlider(tmp, 'N'+str(d+1),  0, self._shp[d]-1, valinit=0,valfmt=u'%2i'))
			self._sliders[d].on_changed(self._update_slider)


		show()

	def _display_image( self, ax = None ):

		if ax == None:		
			gs = gridspec.GridSpec( 1, 1)
			ax    = subplot( gs[0:1,0:1] )

		ax.imshow( self._images[ tuple(self._current_slice) ], cmap = cm.Greys_r  )
		ax.set_xticks( [] )
		ax.set_yticks( [] )
		
		
			

	def _sel_data_set(self, lable ):
		self._curent_indx = self._radio_dic[ lable ]

	def _update_slider(self, val):
		
		for i, sl in enumerate(self._sliders):
			self._current_slice[i] = int( sl.val )
		
			
	
   		self._display_image()



		
if __name__ == '__main__':
	from SiemensReader import SiemensReader
	from CEST_MRI import CEST_MRI	
	from matplotlib.pyplot import imshow, plot, show
	from numpy import fft, sum, shape
	import pickle

	path1 = '/Users/josh/Documents/Data/test_sets/phantum-June10_2014/CEST720/meas_MID98_Gre_2D_P_delayGauss_10x720d_185V_st2us_FID25621.dat'
	path2 = '/Users/josh/Documents/Data/test_sets/phantum-June10_2014/CEST720/meas_MID99_Gre_2D_TD_ref_12ms_9d_FID25622.dat'  

	#ref = SiemensReader(path2)
	dat = SiemensReader(path1)
	
	#imgs = abs(dat.images())
	#imgs = imgs.astype(dtype='float')
	#imgs = squeeze(sum(imgs, axis=0 ))
	#print shape(imgs)
	#pickle.dump( imgs, open( "save.p", "wb" ) )
	

	imgs = pickle.load( open( "save.p", "rb" ) )
	myMRI = CEST_MRI( imgs )
	

	#out = myMRI.whiten( ref = ref.images() )

	#for i in range(51):
	#	imshow( out[i] )
	
	#	show()




