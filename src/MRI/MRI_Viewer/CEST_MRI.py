import sys
import os
#sys.path.insert(0, '../..')

import MRI

import ROI

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import RadioButtons, RectangleSelector


import numpy as np
from numpy import shape, squeeze, zeros, zeros_like, amax
from matplotlib.pyplot import imshow, plot, show, subplot


import pickle

class CEST_MRI( MRI.MRI ):
	def __init__( self, img ):		

		self._radio_dic = {'red':3, 'blue':0, 'green':1,' black':2}
		self._c_indx = 0
		self.cest = zeros( (4, shape(img)[0]), dtype='float')

		self.box_draw = {}
		self.select_map = {}

		self._ref_img = zeros_like( img )
		ref =  amax( img, axis=0 ) 
		for i, d in enumerate(img ):
			self._ref_img[i] = img[i]/ref	
		super(CEST_MRI, self).__init__( squeeze(img) )

	def __del__(self):
		self.dump_data()

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
	
	def eval_mask( self, msk = None ):
		"""
		INPUT a masked array of the same size as the last two dimensions
		of the CEST dataset

		Returns the averaged CEST signal over all points in the masked array	
		"""
		if msk == None: msk = self.select_map[  self._c_indx ]
		
		thr_msk = self._ref_img * msk.mask
		v = ( thr_msk.mean( axis=-2).mean(axis = -1 ))
		self.cest[ self._c_indx] = v

		return v	
		
	def _sel_data_set(self, lable ):
		self._c_indx = self._radio_dic[ lable ]	

	def _selector_callback(self,eclick, erelease):
		self.ax_main.clear()
		super(CEST_MRI, self)._display_image( self.ax_main )
    		x1, y1 = eclick.xdata, eclick.ydata
    		x2, y2 = erelease.xdata, erelease.ydata
		
		dx = [ eclick.xdata, erelease.xdata ]
		dy = [ eclick.ydata, erelease.ydata ]
	
		roi = ROI.ROI( self._ref_img  )
		rec_patch = roi.rectangle( (x1,y1), x2-x1 , y2-y1 )	
	
		self.select_map[  self._c_indx ] = roi
		self.box_draw[self._c_indx] = rec_patch

		self.eval_mask()
			
		 
		self.ax_cest.clear()
		for i, spec in enumerate(self.cest):
			k = { 3:'r', 0:'b', 1:'g', 2:'k'}[i]
			self.ax_cest.plot( spec, k )


		self.ax_main.clear()
		super(CEST_MRI, self)._display_image( self.ax_main )
		kc = { 3:'r', 0:'b', 1:'g', 2:'k'}

		for k in self.select_map.keys():
			s= self.box_draw[k]
			s.set_ec( kc[k] )
			self.ax_main.add_patch( s )
		self._fig.canvas.draw()


		self._selector = RectangleSelector(self.ax_main, self._selector_callback,
                             drawtype='box', useblit=False,
                             button=[1], # don't use middle button
                             minspanx=5, minspany=5,
                             spancoords='pixels')

		self.dump_data()

	
	def dump_data( self, name='./cest_data.p'):

		f = open( name, "wb" )
		pickle.dump( {'cest':self.cest,'sel':self.select_map}, f )
		f.close()

		
if __name__ == '__main__':
	path = '/Users/josh/Documents/Data/test_sets/phantum-June10_2014/CEST720/meas_MID98_Gre_2D_P_delayGauss_10x720d_185V_st2us_FID25621.dat'
	#path = '/Users/josh/Documents/Data/test_sets/phantum-June10_2014/CEST720/meas_MID99_Gre_2D_TD_ref_12ms_9d_FID25622.dat' 

	path = '/Users/josh//Documents/Data/test_sets/phantum-June10_2014/TRECEST720/meas_MID100_Gre_2D_P1331Gauss_10x720d_st2us_FID25623.dat' 
	import numpy as np
	import Process_MRI
	from FileIO.SiemensReader import SiemensReader
	

	c = SiemensReader(path)
	fids = c.read()

	img = Process_MRI.Process_MRI( fids )
	CEST_MRI(img)
