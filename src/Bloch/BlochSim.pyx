import sys
sys.path.append('../../src/')

import numpy as np
import scipy as sp

import scipy.linalg

cimport numpy as np


cimport cython

import math

import numpy as np

import PulseSeq.PulseElem
import PulseSeq.PulseSeq

import src.Spin
import src.util.Gamma



cdef class BlochSim( object ):
	cdef np.float __B0
	cdef list __spins
	cdef np.ndarray  __Rmtx

	def __init__( self, *arg, float B0 = 9.4 ):
		if not all( isinstance(i, src.Spin.Spin) for i in arg): 
			raise ValueError("A Non-spin object was passed to the BlochSim Constructor!")
		
		cdef np.ndarray[np.float, ndim=2] __Rmtx = self.__Rmtx

		self.__B0      = B0
		self.__spins   = arg
		self.__Rmtx    = np.zeros( (len(arg),len(arg) ), dtype=np.float )
	
	def add_kex(self, frm, to, kex ):
			
		frm, to = frm-1, to-1		
		
		cfct =  (self.__spins[frm]).c/(self.__spins[to]).c 
		
		self.__Rmtx[frm,to]= kex/( 1. + cfct )
		self.__Rmtx[to,frm]= kex - self.__Rmtx[frm,to]

		self.__Rmtx[frm,frm] = 0.0
		self.__Rmtx[to, to ] = 0.0	
		self.__Rmtx[frm,frm] = sum( self.__Rmtx[frm] )
		self.__Rmtx[to, to ] = sum( self.__Rmtx[to ] )

	@cython.boundscheck(False) 
	cdef np.ndarray[np.float64_t, ndim=2] dM( self, object pulse ):
		cdef int sz, i, ii, x, y, z
		cdef float w, R1, R2, c		
		cdef np.ndarray[np.float64_t, ndim=2] dM
		
		sz = len( self.__spins)
		dM  = np.zeros( ( 3*sz+1, 3*sz+1) ) #the transformation matrix 

		for i in range(sz): 
			x = (3*(i))
			y = (3*(i))+1
			z = (3*(i))+2
				
			w  = (self.__spins[i]).x0 * (self.__spins[i]).gamma * self.__B0
			R1 = (self.__spins[i]).R1
			R2 = (self.__spins[i]).R2
			c  = (self.__spins[i]).c
				
			if pulse.gamma  == (self.__spins[i]).gamma:
				B1 = 2 * math.pi * pulse.B1  
				dw = pulse.offset
			else:
				B1 = [0.0, 0.0]
				dw = 0.0

			# Terms describing spin under influence of B0 & B1
			dM[x,x] = -( R2 + self.__Rmtx[i,i] )
			dM[x,y] = -2*math.pi*( w - dw )
			dM[x,z] =  B1[1]
                
			dM[y,y] = -( R2 + self.__Rmtx[i,i] )
			dM[y,x] =  2*math.pi*( w - dw )
			dM[y,z] = -B1[0]
                
			dM[z,z] = -(R1 + self.__Rmtx[i,i] )
			dM[z,x] = -B1[1]
			dM[z,y] =  B1[0]

			dM[z,-1]=  c*R1


			#Terms describing kex between spins
			for ii in range(sz):
				if ii == i: continue
				
				ox = (3*ii)
				oy = (3*ii)+1
				oz = (3*ii)+2
				
				dM[ox, x] = self.__Rmtx[i,ii]
				dM[oy, y] = self.__Rmtx[i,ii]
				dM[oz, z] = self.__Rmtx[i,ii]

		return dM

									
	@cython.boundscheck(False) 
	cpdef run( self, object pulseSeq ):
		
		cdef np.ndarray[np.float64_t, ndim=1] I0
		cdef np.ndarray[np.float64_t, ndim=2] M, I
		cdef int sz, i

		I0 = np.concatenate( [ s.v for s in self.__spins] )
		I0 = np.concatenate( (I0,[1.0]))
		
		sz = len( self.__spins)
		
		I = np.zeros( (pulseSeq.ncyc, 	3*sz+1) )	

		M = np.identity( 3*sz+1 )


		i=0
		for p in pulseSeq:

			I0 = (p.M(self)).dot( I0 )
 	
			if isinstance( p,  PulseSeq.PulseElem.Acq ):
				I[i] = I0

				#reset I0
				I0 = np.concatenate( [ s.v for s in self.__spins] )
				I0 = np.concatenate( (I0,[1.0]))

				i+=1
	
		return np.squeeze(I)
