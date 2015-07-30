import abc

import os
import numpy as np


"""
FileReader is an ABC for reading in MR files.
Derived methods must impliment
	1) __iter__(self):
		That iterates through all the fids in the file
	2) guess_shape(self):
		The guesses the format of the input form the binary/accesory files
		The return value must list the dimension sizes in terms of [ slowest -> fastest]
		e.g. 
			Fastest is size of fid (#complex pts)   nps
			2nd fastest is ni in 2d dim             ni2
			3rd fastest is n2 in 3d dim             ni3
			
			out should be     [ ni3, ni2, nps]

The method
	1) sort(self, shape, dtypes)
		Assembles the fids in the object into the matrix defined by "shape" 
		(supplied from guess shape only overide if guess is wrong)
		
		inputs: 
		1) shape : an array listing the object dimensions in order of slowest changeing to				  fastest. The fastest, last, dimension should be the number of complex pts
			   in the fid.
		2) dtype: an array of length len(shape)-1 containing the keywords 'complex' or 'None'
			  (the last dimension #pts in fid is assumed to be complex)
			  If set to complex it assumes points in that dimension are collected as
			      1st -- real  
			      2nd -- immg
			      3rd -- real
			      ect. 
			if not just overide this in the concreete class!
"""

class FileReader( object ):

	def __init__(self, path,):
		self.npts       = 0	 #number of fids in the file		
		self.pts_in_fid = 0
		self.cplx = np.vectorize( complex ) 

		self.file_size = os.path.getsize(path)
		self.binary = open( path, 'rb')
		
		

	@abc.abstractmethod
	def __iter__(self):	
		"""
		self.__iter__()
		=========================================================  
		A generator returning a complex np.array containg each fid in the
		the order it was collected. 
		"""
		return

	
	@abc.abstractmethod
	def guess_shape( self ):
		"""
		self.guess_shape()
		=========================================================  
		Returns the best guess of the number of points in each dimension
		of the experiment. Function returns a tuple with two entries 
			1) The np.array containg the shape of the data matrix. 
			   In order of slowest changing dimension to the fastest 
                           (C-like index ordering)
			2) Boolean list, indicating if a given dimension is real (True)
			   or real+imaginary ( False ).
		"""
		return

	@abc.abstractmethod
	def read_raw( self ):
		"""
		self.read_raw()
		=========================================================  
		Returns a numpy array containing all the fids in the file.
		The shape is [npts, pts_in_each_fid].
		"""
		return


	def read( self, shp=None,  isr = None ):

		s,r = self.guess_shape()		
		if shp == None: shp = s
		if isr == None: isr = r	

		fids = np.empty( shp, dtype='complex')

		shp = np.delete( shp, -1 )
		indx=[ 0 for i in range(len(shp) ) ]
		
		for num, f in enumerate(self):
			for s, dsize in enumerate(shp[::-1]):	
				num, ans = divmod( num, dsize )

				if not isr[s]:
					if dsize % 2 != 0: raise ValueError
					off = ans // 2
					cntr = dsize/2

					if ans % 2 == 0: ans = cntr - off -1
					else:            ans = cntr + off 
					
				indx[s] = int(ans)
			
			indx = indx[::-1]
			fids[ tuple(indx) ] = f	

		return fids

