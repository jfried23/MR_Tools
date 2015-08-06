import sys
sys.path.insert(0, '../')

from FileIO.SiemensReader import SiemensReader
from MRI import MRI

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import RadioButtons, RectangleSelector

from numpy import sum, shape, squeeze, zeros, zeros_like, amax
from matplotlib.pyplot import imshow, plot, show, subplot

import pickle

class CEST_MRI( MRI ):
	def __init__( self, img ):
		self._radio_dic = {'red':3, 'blue':0, 'green':1,'black':2}
		self._c_indx = 0
		self._cest_spec = zeros( (4, shape(img)[0]), dtype='float')

		self._ref_img = zeros_like( img )
		ref =  amax( img, axis=0 ) 
		for i, d in enumerate(img ):
			self._ref_img[i] = img[i]/ref	
		super(CEST_MRI, self).__init__( squeeze(img) )

	def _display_image( self ):
		gs = gridspec.GridSpec( 5, 10)
		self.ax_cest    = subplot( gs[0:4,4:10] )
		self.ax_main    = subplot( gs[0:4,0:4] )

		axcolor = 'lightgoldenrodyellow'
		rax = plt.axes([0.05, 0.1, 0.15, 0.15], axisbg=axcolor)
		self._radio = RadioButtons(rax, self._radio_dic.keys() , active=0)
		self._radio.on_clicked( self._sel_data_set )

		self._selector = RectangleSelector(self.ax_main, self._selector_callback,
                             drawtype='box', useblit=False,
                             button=[1], # don't use middle button
                             minspanx=5, minspany=5,
                             spancoords='pixels')

		super(CEST_MRI, self)._display_image( self.ax_main )
			
	def _sel_data_set(self, lable ):
		self._c_indx = self._radio_dic[ lable ]	

	def _selector_callback(self,eclick, erelease):
    		x1, y1 = eclick.xdata, eclick.ydata
    		x2, y2 = erelease.xdata, erelease.ydata
		
		dx = [ eclick.xdata, erelease.xdata ]
		dy = [ eclick.ydata, erelease.ydata ]

		self._cest_spec[ self._c_indx] = sum( self._ref_img[:,min(dy):max(dy), min(dx):max(dx)], axis=(-1,-2) )/(abs(x1-x2)*abs(y1-y2))
		 
		self.ax_cest.clear()
		for i, spec in enumerate(self._cest_spec):
			k = { 3:'r', 0:'b', 1:'g', 2:'k'}[i]
			self.ax_cest.plot( spec, k )
		self._fig.canvas.draw()
	
		self.dump_data()

	def dump_data( self, name='./cest_data.p'):
		f = open( name, "wb" )
		pickle.dump( self._cest_spec, f )
		f.close()

		
if __name__ == '__main__':
	path = '/Users/josh/Documents/Data/test_sets/phantum-June10_2014/CEST720/meas_MID98_Gre_2D_P_delayGauss_10x720d_185V_st2us_FID25621.dat'
	#path = '/Users/josh/Documents/Data/test_sets/phantum-June10_2014/CEST720/meas_MID99_Gre_2D_TD_ref_12ms_9d_FID25622.dat' 

	path = '/Users/josh//Documents/Data/test_sets/phantum-June10_2014/TRECEST720/meas_MID100_Gre_2D_P1331Gauss_10x720d_st2us_FID25623.dat' 
	import numpy as np
	import Process_MRI

	c = SiemensReader(path)
	fids = c.read()

	img = Process_MRI.Process_MRI( fids )
	CEST_MRI(img)
