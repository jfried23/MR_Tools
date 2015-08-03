import sys
sys.path.append('../')

import util.Gamma 

from numpy import array

class Spin:
	"""
		Spin Class
		
		A contaniner for representing a single Spin 
		
		history gets set by a Bloch simulation object
		self.history[0] = complex array of transverse magnitization
				  [ [ x + iy]_1
				    [ x + iy]_2
				     ...
				    [ x + iy]_n
				  ]
		self.history[1] = float of z magnitization
				  [ array_1
				    array_2
				     ...
				    array_n
                                   ] 
	"""
	def __init__(self, R1, R2, x0, c, atm='1H'):
		self.R1   = float(R1)
		self.R2   = float(R2)
		self.x0   = float(x0) 
		self.c    = float(c)
		self.gamma = util.Gamma.gamma_MHz[atm]
	
		self.v    = array( [0.0, 0.0, self.c], dtype=float )

		self.__bounds = { 'R1':None, 'R2':None, 'x0':None, 'c':None }

		self.history= None


	def __str__( self ):
		s  = 'Spin object:\n' 	
		s += '\tGy %f Mhz/T.\n' % (self.gamma)   
		s += '\tx0 %f ppm.\n' %(self.x0)
		s += '\tR1:  %f\n\tR2:  %f\n\tc:   %f\n' %(self.R1,self.R2, self.c)
		s += '\tv=[ %f, %f, %f]\n' %(self.v[0],self.v[1],self.v[2])
		return s

	
	def set_lb_ub( self, key, lb = None, ub = None ):
		if key in self.__bounds: self.__bounds[key] = ( lb, ub )	
		else: raise ValueError


	def get_opt_vals( self ):

		x   = []
		lb  = []
		ub  = []

		for key in sorted(self.__bounds.keys() ):
			if isinstance( self.__bounds[key], type(None) ): continue

			b = self.__bounds[key]

			x.append( getattr( self, key ) )
			lb.append( b[0] )
			ub.append( b[1] )

		return x, lb, ub 


	def set_opt_vals( self, x ):
		"""
		Takes a list of floats and pops from the front
		"""
		for key in sorted(self.__bounds.keys() ):
			if isinstance( self.__bounds[key], type(None) ): continue
			else: 
				setattr( self, key, x.pop(0) )
				

if __name__ == '__main__':
	h=Spin(1,1,3.5, 30)
	h.set_lb_ub( 'R1', 0.1, 5 )
	h.set_lb_ub( 'R2', None, 500 )

	x,lb,ub=h.get_opt_vals()

	x[1] = 9999
	x.append('Josh')
	x.append( 5555 )

