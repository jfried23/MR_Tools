from FileReader import FileReader

from string import split 
from math import floor

from struct import calcsize, unpack

import numpy as np

import os
import glob
import struct


class BrukerReader( FileReader ):

	def __init__(self, path):
		FileReader.__init__(self, path)

		acqs_name = glob.glob( os.path.join(os.path.dirname(path), 'acqu*s') ) 
		self.npts = 1
		#First build a dictonary for each dimension in the experiment
		self.acq={}
 
		for dim in acqs_name:
			if 'acqus'  in dim: n = 1 
			if 'acqu2s' in dim: n = 2
			if 'acqu3s' in dim: n = 3
			if 'acqu4s' in dim: n = 4
			if 'acqu5s' in dim: n = 5

	
			dic = {}	
			for line in open( dim, 'r' ):
				if   'TD='       in line:
					dic['TD'] = int(split(line)[-1])
					if n==1:
						dic['TD'] /= 2
						self.pts_in_fid = dic['TD']
					else: self.npts *= dic['TD']

				elif 'SW_h='     in line: dic['SW_h'] = float(split(line)[-1])
				elif 'PULPROG='  in line: dic['PULPROG']= split(line)[-1]
				elif 'DTYPA='    in line: dic['DTYPA'] = { 0:'i', 1:'d'}[ int(split(line)[-1])]
				elif 'BYTORDA='  in line: dic['BYTORDA']={ 0:'<', 1:'>'}[ int(split(line)[-1])]	
				elif 'GRPDLY='   in line: dic['GRPDLY']= int(split(line)[-1])
				elif 'DECIM='    in line: dic['DECIM']=  int(split(line)[-1])
				elif 'DSPFVS='   in line: dic['DSPFVS']= int(split(line)[-1])

			self.acq[n]=dic
			
		#Now build the format string for reading in the data
		self.fmt = '%s%d%s' % (self.acq[1]['BYTORDA'],self.pts_in_fid*2, self.acq[1]['DTYPA'])

		#In bruker ser, each block must end at integer of 1024 
		#check to see what correction is needed
		blk_sz = calcsize(self.fmt)
		if blk_sz % 1024:  self.cor = 1024 - (blk_sz % 1024)   #correction needed
		else: self.cor = 0			                  #no correction

		#Now prepare for the dreaded digital filtering
		if ('GRPDLY' in self.acq[1] ) and self.acq[1]['GRPDLY'] != -1:
			digshift = acq[1]['GRPDLY'] 
	
		else:
			indx1,indx2 = self.acq[1]['DECIM'], self.acq[1]['DSPFVS'] - 10
			digshift = _brukdigital[indx1][indx2]
		
		self.roll_npts  = int(floor(digshift))+1 
		self.phase = digshift - self.roll_npts 

	def __iter__(self):
		self.binary.seek( 0 )
		fid = np.zeros( (self.pts_in_fid/2 ), dtype=complex )

		for i in range( self.npts ):
			t = unpack( self.fmt, self.binary.read(calcsize(self.fmt ))  )
			self.binary.seek(self.cor, 1) #Now skip the correction factor

			fid = self.cplx( t[0::2], t[1::2] )
			fid[0:self.roll_npts] = fid[-1]
			fid = np.roll( fid, -1*self.roll_npts )

			fid = np.multiply( fid, np.exp(1.0j*self.phase*np.arange(len(fid))/len(fid)) )
			
			yield fid	
			
			
	def read_raw( self ):
		self.fids = np.empty( (self.npts, self.pts_in_fid ) , dtype=complex )
		for i, f in enumerate(self): self.fids[i] = f
		self.binary.close()
		return self.fids
	

	def guess_shape( self ):
		v = [ self.acq[k]['TD'] for k in sorted(self.acq.keys())[::-1] ]
		r = [ True for r in v]
		r[-1] = False
		return (v, r )	


########### BRUKER ARBITRARY NUMBERS DO NOT CHANGE!!!!!!!!!! ############################
# https://ucdb.googlecode.com/hg/application/ProSpectND/html/dmx_digital_filters.html
#             DECIM              DSPFVS 10       DSPFVS 11      DSPFVS 12 
_brukdigital= { 2:          (    44.7500,         46.0000,        46.311   ),
   		3:          (    33.5000,         36.5000,        36.530   ),
   		4:          (    66.6250,         48.0000,        47.870   ),
		6:          (    59.0833,         50.1667,        50.229   ),
		8:          (    68.5625,         53.2500,        53.289   ),
  	       12:          (    60.3750,         69.5000,        69.551   ),
  	       16:          (    69.5313,         72.2500,        71.600   ),
               24:          (    61.0208,         70.1667,        70.184   ),	
               32:          (    70.0156,         72.7500,        72.138   ),	
               48:          (    61.3438,         70.5000,        70.528   ),
               64:          (    70.2578,         73.0000,        72.348   ),	
               96:          (    61.5052,         70.6667,        70.700   ),
              128:          (    70.3789,         72.5000,        72.524   ),	
              192:          (    61.5859,         71.3333		   ),
              256:          (    70.4395,         72.2500		   ),
              384:          (    61.6263,         71.6667		   ),
              512:          (    70.4697,         72.1250		   ),
              768:          (    61.6465,         71.8333		   ),
             1024:          (    70.4849,         72.0625		   ),
             1536:          (    61.6566,         71.9167		   ),		
             2048:          (    70.4924,         72.0313		   )    }
#########################################################################################
    
if __name__ == '__main__':
	import matplotlib.pyplot as plt
	from numpy import fft
	path = '/Users/josh/Documents/Data/test_sets/BSA_Power_Dependence/3/ser'
	#path = '/Users/josh/Documents/Data/test_sets/BSA_Power_Dependence/7/ser'
	#path = '/Users/josh/Documents/Data/test_sets/BSA_Power_Dependence/1/fid'

	b = BrukerReader(path)
	#fid = b.getNext()
	
	#plot(fid)
	#show()
	#i=0
	for i, fid in enumerate( b.read() ):
		plt.plot(abs(fft.fftshift(fft.fft(fid))))


	#plt.plot( b.iter2(), 'r' )
	#	i+=1
	#	plt.plot(abs(fft.fftshift(fft.fft(fid))))
	#print i
	plt.show()

	#fid,b = read_bruker(path)
	

	#plot( fft.fftshift(fft.fft(fid[5])) )
	#show()
	#print shape(fid)
	#print b
	
