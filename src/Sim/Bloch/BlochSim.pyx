import sys
sys.path.append('../../')

import numpy as np
import scipy as sp

import scipy.linalg

cimport numpy as np
cimport cython

import array

import math
import numpy as np

import PulseSeq.PulseElem
import PulseSeq.PulseSeq

import Sim.Spin
import util.Gamma


cdef class BlochSim( object ):
	"""
	BlochSim objects are for preforming a bloch simulation of a pulse sequence.
	
	Constructor:  BlochSim( *arg, **kwargs )
		arg -- 	A variable number of Spin objects to model in the Bloch Sim
	     kwargs --  'B0' == B, the static magnetic field strength in Tesla. If
			           not provided B0 = 9.4 T (~400MHz) by default.

	Public Paramters
	----------------
	B0	: The static magnetic field strengh of the simulation in Tesla

	spins   : A list of spin objects to be inclded in the simulation


	Private Paramters
	----------------
	Rmtx	: A square matrix of length == number of spins in the simulation.
		  Desbribes the mass balanced exchnage rates between spins.

	"""



	cdef public float B0
	cdef public spins
	cdef np.float_t [:,::1] Rmtx

	def __init__( self, *arg, **kwargs):

		self.B0      = 9.4	
		self.spins   = arg
		self.Rmtx    = np.zeros( (len(arg),len(arg) ) )
	
		for k,v in kwargs.items(): 
			if k == 'B0': self.B0 = v
			else: raise ValueError


	def add_kex(self, int frm, int to, float kex ):
		"""
			add_kex( int from, int to, float kex)

			Establish a magnitization exchnage pathway between spin 'from'
			and spin number 'to' of kex s^-1. 
			Spin number is 1 indexed, not 0 and determined by the order in which
			it was inclded in the input args 
		"""
		cdef float cfct		
	
		frm, to = frm-1, to-1		
		
		cfct =  (self.spins[frm]).c/(self.spins[to]).c 
		
		self.Rmtx[frm,to]= kex/( 1. + cfct )
		self.Rmtx[to,frm]= kex - self.Rmtx[frm,to]

		self.Rmtx[frm,frm] = 0.0
		self.Rmtx[to, to ] = 0.0	
		self.Rmtx[frm,frm] = np.sum( self.Rmtx[frm] )
		self.Rmtx[to, to ] = np.sum( self.Rmtx[to ] )
		
		

	@cython.boundscheck(False) 
	cpdef np.ndarray[ np.float_t, ndim=2] dM( self, object pulse ):
		"""
			dM( PulseElem )

			Takes a PulseElem object and returns a dM matrix descrbing the evolution
			of the spins in the BlochSim object under the influence of the PulseElem

			returns a square matrix of size n x (Ix Iy Iz) + 1
		"""
		cdef int sz, i, ii, x, y, z
		cdef float w, R1, R2, c		
		cdef np.ndarray[ np.float_t, ndim=2] dM	

		sz = len( self.spins )
		dM  = np.zeros( ( 3*sz+1, 3*sz+1) ) #the transformation matrix 

		for i in range(sz): 
			x = (3*(i))
			y = (3*(i))+1
			z = (3*(i))+2
				
			w  = (self.spins[i]).x0 * (self.spins[i]).gamma * self.B0
			R1 = (self.spins[i]).R1
			R2 = (self.spins[i]).R2
			c  = (self.spins[i]).c
				
			if pulse.gamma  == (self.spins[i]).gamma:
				B1 = 2 * math.pi * pulse.B1  
				dw = pulse.offset
			else:
				B1 = [0.0, 0.0]
				dw = 0.0

			# Terms describing spin under influence of B0 & B1
			dM[x,x] = -( R2 + self.Rmtx[i,i] )
			dM[x,y] = -2*math.pi*( w - dw )
			dM[x,z] =  B1[1]
                
			dM[y,y] = -( R2 + self.Rmtx[i,i] )
			dM[y,x] =  2*math.pi*( w - dw )
			dM[y,z] = -B1[0]
                
			dM[z,z] = -(R1 + self.Rmtx[i,i] )
			dM[z,x] = -B1[1]
			dM[z,y] =  B1[0]

			dM[z,-1]=  c*R1


			#Terms describing kex between spins
			for ii in range(sz):
				if ii == i: continue
				
				ox = (3*ii)
				oy = (3*ii)+1
				oz = (3*ii)+2
				
				dM[ox, x] = self.Rmtx[i,ii]
				dM[oy, y] = self.Rmtx[i,ii]
				dM[oz, z] = self.Rmtx[i,ii]

		return dM

									
	@cython.boundscheck(False) 
	cpdef run( self, object pulseSeq ):
		"""
			run( PulseSeq )
			
			Runs a Bloch simulation defined by the supplied Pulse sequence.
			
			The result of the simulation is receored in the history member of each spin
			and is returned by this function

		"""

		
		cdef np.ndarray[ np.float_t, ndim=1]  I0
		cdef np.ndarray[ np.float_t, ndim=2]  M, I
		cdef int sz, i

		I0 = np.concatenate( [ np.asarray(s.v) for s in self.spins] )
		I0 = np.concatenate( (I0,[1.0]))
		
		sz = len( self.spins)
		
		I = np.zeros( (pulseSeq.ncyc, 	3*sz+1) )	

		M = np.identity( 3*sz+1 )
		

		i=0
		for p in pulseSeq:

			I0 = (p.M(self)).dot( I0 )
 	
			if isinstance( p,  PulseSeq.PulseElem.Acq ):
				I[i,:] = I0

				#reset I0
				I0 = np.concatenate( [ s.v for s in self.spins] )
				I0 = np.concatenate( (I0,[1.0]))

				i+=1


		h = np.array_split( I[:,0:-1], len(self.spins), axis=1 )
	
	
		for sp_num, s in enumerate( self.spins ):
			s.history = h[sp_num]

		return np.squeeze(I)
