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

"""

class FileReader( object ):

	def __init__(self, path,):
		self.npts       = 0	 #number of fids in the file		
		self.pts_in_fid = 0
		self.cplx = np.vectorize( complex ) 

		self.file_size = os.path.getsize(path)
		self.path = path	
		self.fids = None
		

	@abc.abstractmethod
	def __iter__(self):	
		"""
		self.__iter__()
		=========================================================  
		A generator returning a complex np.array for each fid in the
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
		The shape is [npts, pts_in_each_fid]. Use "read" function
		to order and sort the fids for processing.
		"""
		return


	def read( self, shp=None,  isr = None ):
		"""
		self.read( shape = None, is_real = None )
		=========================================================  
		
		Parameters: 
				shape -- array of type 'int'
					 describes of the number of dimensions 
					 in the experiment. Slowest changing
					 dimension to the fastest.
					 eg [ ni2, ni1, pts_in_fid]
					 if 'None' the shape provided by "guess_shape" is used.

		              is_real -- array of type 'bool'
					 describes if a given dimension is real or complex.
				         if 'False' elements in that dimension will be
				         sorted as real + imaginary pairs in the output fid.
	                                 MUST BE SAME LENGTH AS SHAPE.
					 if 'None', is_real from "guess_shape" is used.
		"""


		s,r = self.guess_shape()
		
		if isinstance( shp, type(None) ): 
			shp = s 
			print "\nAssuming {0}D experiment, of shape {1}".format( len(shp), shp)
		if isinstance( isr, type(None) ): 
			isr = r 
			key = {True:'Real',False:'Real+Complex'}
			print "\nAssuming dimensions are {1}\n".format( len(shp), [key[v] for v in isr]) 

		isr = np.delete(isr,-1) #remember last axis is always complex

		self.read_raw()
		self.fids.resize( shp )

		#Go through each axis determin if it is real (True) or real+complex (False)

		for ax, v in enumerate(isr):
			if v == True: continue
			#only sort the complex dimensions
			
			i_index = range(0,shp[ax],2)
			r_index = range(1,shp[ax],2)[::-1]
			
			rls = np.take( self.fids, r_index, axis=ax)
			ims = np.take( self.fids, i_index, axis=ax)

			self.fids = np.concatenate( (rls, ims), axis=ax)

		return self.fids
