from FileReader import FileReader
from numpy import shape,rollaxis,fft, absolute, zeros_like, ndarray

import struct
import numpy as np

def Siemiens_Global_xml( path ):

	binary = open(path, 'rb' )
	offset, nevps = struct.unpack( '2I', binary.read(struct.calcsize('2I')))
	s=''
	for i in range( nevps):
        	name = _read_cstr(binary)
		size = struct.unpack( 'I', binary.read(struct.calcsize('I')))
		s += ''.join(map(chr, struct.unpack( '<%dB'%(size[0]), binary.read(struct.calcsize('<%dB'%(size[0]))))))
	binary.close()
	return s

	
class SiemensReader( FileReader ):

	def __init__(self, path):
		FileReader.__init__(self, path)

		self._glb_header_sz    = struct.unpack( 'I', self.binary.read(struct.calcsize('I')))[0]
		self._local_header_sz  = 128

		self.binary.seek( self._glb_header_sz )		
			
		self.read_headers( )

		

	def __iter__(self):
		self.binary.seek( self._glb_header_sz  )
		fid = np.empty( (self.pts_in_fid ) , dtype=complex )

		for i in range( self.npts ):
			if self.binary.tell() > self.file_size: break
			self.binary.seek( self._local_header_sz, 1 )	
			fid = np.fromfile( self.binary, dtype=np.complex64, count = self.pts_in_fid)

			yield fid


	def read_headers( self  ):
		self.binary.seek( self._glb_header_sz  )
		
		self.pts_in_fid  = 0		
	
		self.npts =0	
		while self.binary.tell() < self.file_size:
			hdr = SiemensHeaderObj( self.binary.read( self._local_header_sz ) )
			self.binary.seek( hdr.samples_in_scan*2*4, 1 )
			
			if self.npts == 0: 
				self.pts_in_fid = hdr.samples_in_scan
				self.hdr = hdr

			if hdr.is_end: break
			assert(  hdr.samples_in_scan == self.pts_in_fid  )
			self.npts +=1

	def read( self, shp = None, isr = None ):
		FileReader.read( self, shp, isr )
		self.fids = np.rollaxis( self.fids, -2 )
		
		print "\nPlacing channel dimension first. New shape: {0}".format( self.fids.shape ) 
		
		return self.fids

	def read_raw( self ):
		self.fids = np.empty( (self.npts, self.pts_in_fid ) , dtype=complex )
		for i,f in enumerate(self): self.fids[i] = f
		self.binary.close()	
		return self.fids		

	def guess_shape( self ):
		"""
		Return the shape in F_order indexing (last dim changing the fastest )
		"""
		k0       =  self.hdr.samples_in_scan
		k1       =  self.hdr.kspace_center * 2
		chns     =  self.hdr.used_channels
		k2       =  self.npts/(chns*k1)
	
		
		assert k1 * chns * k2 == self.npts 
		
		shape  = np.array([ k2, k1, chns, k0 ],dtype=int)
		is_real = [ True, False, True, False ]	
		return shape, is_real


class SiemensHeaderObj( object ):
	def __init__(self, hdr_binary):
		header = struct.unpack( '1H2B1i3I8B20H1f11H7f2H', hdr_binary )

		self.dma_length      	= header[0]
		self.dma_flags       	= header[1:3]
		self.userid          	= header[3]
		self.scan_num        	= header[4]
		self.timestamp       	= header[5]
		self.pmutimestamp    	= header[6]
		self.eval_info_mask  	= header[7:15]
		self.samples_in_scan 	= header[15]       #Real pts in FID
		self.used_channels   	= header[16]       #Recieve channels
		self.line            	= header[17]
		self.acquisition     	= header[18]
		self.slice	     	= header[19]
		self.partition       	= header[20]
		self.echo            	= header[21]
		self.phase           	= header[22]
		self.repetition      	= header[23]
		self.set             	= header[24]
		self.segment         	= header[25]
		self.ida	    	= header[26]
		self.idb             	= header[27]
		self.idc             	= header[28]
		self.idd                = header[29]
		self.ide                = header[30]
		self.pre	        = header[31]
		self.post               = header[32]
		self.kspace_center      = header[33]
		self.coil_select        = header[34]
		self.readout_off_center = header[35]
		self.time_since_last_rf = header[36]
		self.kspace_center_line = header[37]
		self.kspace_center_part = header[38]
		self.ice_parameters     = header[39:43]
		self.free_parameters    = header[43:47]
		self.saggital           = header[47]
		self.coronal		= header[48]
		self.transverse         = header[49]
		self.quaternion         = header[50:54]
		self.channel_id         = header[54]
		self.ptab_pos_neg       = header[55]

		self.is_end             = bool( header[7] & 0x01 )

	

if __name__ == '__main__':
	from SiemensReader import SiemensReader
	path = '/Users/josh/Documents/Data/test_sets/phantum-June10_2014/CEST720/meas_MID98_Gre_2D_P_delayGauss_10x720d_185V_st2us_FID25621.dat'
	#path = '/Users/josh/Documents/Data/test_sets/phantum-June10_2014/CEST720/meas_MID99_Gre_2D_TD_ref_12ms_9d_FID25622.dat'  

	
	import matplotlib.pyplot as plt
	from numpy import fft
	import numpy	
	import seaborn
	from skimage import feature

	s = SiemensReader(path)
	f = s.read()
	
	img = abs(fft.fftshift(fft.fft2(f[18,0])))
	plt.imshow( img, cmap = plt.cm.gray )
	plt.show()
	

	#print np.squeeze(f).shape

	#print f.shape	
	##edges1 = feature.canny(img, sigma=.5)
	
	#plt.imshow( f[22, 0], cmap = plt.cm.gray )
	#plt.show()

	#img = s.images()
	#print shape(img)
	#imshow( img[0][0] ) 

	
	
	#shp = s.guess_shape()
	#print shp	
	#s.sort( )
	#f = s.fids


	#imshow( abs(fft.fftshift(fft.fft2(f[20,0]))) )
	#show()
	#for one in range(24):
	#	imshow( f[one] )
	#imshow( f[0] )
	#show()
	#result = zeros((192,192))
	#for i in range( shp[2] ):
	#	result += numpy.absolute( (fft.fftshift(fft.fft2(f[0,:,i,:]))) )
	#imshow( result )
	#show()


