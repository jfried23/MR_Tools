import sys
sys.path.append('../')

import util.Gamma 

import numpy as np

class Spin:
	"""
	The representation of a Spin species for a MR Simulator.

	Constructor Spin(R1=0, R2=0, x0=0, c=1, atm='1H')	

	Paramters
	----------
	R1	  : Spin lattice relaxation rate. (Units 1/sec) 

	R2	  : Spin-spin transverse relxation rate. (Units 1/sec)

	x0	  : Chemical shift relative to center of spectra. (Units 1.e-6 Hz -- ppm)

	c	  : The concentration or relative amount of signal (Abitrary units)

	atm	  : Name of atom type (string)


	Members
	---------

	history   : An np array descbing the evolution of the magnitization vector durring BlochSim
		    hisotry[:,0] --> evolution of x magnitzation
		    history[:,1] --> evolution of y magnitzation
		    history[:,3] --> evolution of z magnitzation
 
	__bounds  : A dictionary mapping the fitting parameters to tuples of
		    (lower_bound, upper_bound) 
		
	"""

	def __init__(self, R1=0, R2=0, x0=0, c=1, atm='1H'):
		self.__R1   = float(R1)
		self.__R2   = float(R2)
		self.__x0   = float(x0) 
		self.__c    = float(c)
		self.__gama  = util.Gamma.gamma_MHz[atm]
	
		self.__v    = np.array( [0.0, 0.0, self.__c], dtype=float )

		self.__bounds = { 'R1':None, 'R2':None, 'x0':None, 'c':None }

		self.history= None


	def __str__( self ):
		s  = 'Spin object:\n' 	
		s += '\tGy %f Mhz/T.\n' % (self.__gama)   
		s += '\tx0 %f ppm.\n' %(self.__x0)
		s += '\tR1:  %f\n\tR2:  %f\n\tc:   %f\n' %(self.__R1,self.__R2, self.__c)
		s += '\tv=[ %f, %f, %f]\n' %(self.__v[0],self.__v[1],self.__v[2])
		return s

	
	def set_lb_ub( self, key, lb = None, ub = None ):
		"""
		Make a spin paramater optimizable by estabishing upper and lower bonds
		for the solver.

		Parameters
		----------
		key       : String describing the spin member parameter to optimize in fitting
			    e.g. 'R1', 'R2', 'x0','c' 
		lb	  : Float -- the lower bound. If -Inf use None			
		ub        : Float -- the upper bound. If +Inf use None

		Returns
		-------
		---None---
		"""

		if key not in self.__bounds: raise ValueError
		self.__bounds[key] = ( lb, ub )	

	def get_opt_vals( self ):
		"""
		Return a list containg the current values of the 
		optimizable paramaters for this spin object. 

		Parameters
		----------
		---None---

		Returns
		-------
		x	: A list of floats containing the current values of the optimizable
			  paramters.
		"""

		x = []
		for key in sorted(self.__bounds):
			if self.__bounds[key] == None: continue		
			
			elif key == 'R1': x.append(self.__R1) 
			elif key == 'R2': x.append(self.__R2)
			elif key == 'x0': x.append(self.__x0) 
			elif key == 'c' : x.append(self.__c)
			else: raise ValueError 
		return x

	def set_opt_vals( self, vals ):
		"""
		Set the 'n' optimizable parameters of this spin object to the first
		'n' elements in the vals list.

		Parameters
		----------
		vals	: A list of floats containg the current parameter values for this spin obj.

		Returns
		-------
		---None--		
		"""
		
		if type(vals) == np.ndarray: vals = np.ndarray.tolist( vals )

		i=0
		for key in sorted(self.__bounds):
			if self.__bounds[key] == None: continue
			
			elif   key == 'R1': self.__R1 = vals.pop(0)
			elif   key == 'R2': self.__R2 = vals.pop(0)
			elif   key == 'x0': self.__x0 = vals.pop(0)
			elif   key == 'c' : 
				self.__c    = vals.pop(0)
				self.__v    = np.array( [0.0, 0.0, self.__c], dtype=float )
			else: raise ValueError
			i+=1


	def get_opt_limits( self ):
		"""
		Return a list of tuples descbing the (upper, lower) bounds for each paramter
		to optimized for this spin object.

		Parameters
		----------
		---None---
				
		Returns
		-------
		list of tuples ( upper_bound, lower_bound) for each fitting paramter		
		"""
		bounds=[]
		for key in sorted(self.__bounds):
			if self.__bounds[key] != None: bounds.append( self.__bounds[ key ] )
		return bounds
	

	#Getters
	@property
	def v ( self ) : return self.__v
	@property	
	def R1( self ) : return self.__R1
	@property
	def R2( self ) : return self.__R2
	@property
	def x0( self ): return self.__x0
	@property
	def c( self )  : return self.__c
	@property
	def gamma(self): return self.__gama

if __name__ == '__main__':
	h=Spin(1,1,3.5, 30)
	
	h.set_lb_ub( 'R2', 0.0, 200. )
	h.set_lb_ub( 'x0', 2, 5 )


	print h.get_opt_vals()
	print h.get_opt_limits()

	new=range(10)
	new[0]=-8

	h.set_opt_vals( new )

	print h

	print new
