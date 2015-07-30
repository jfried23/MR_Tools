from util import Gamma 
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
		self.__R1   = float(R1)
		self.__R2   = float(R2)
		self.__x0  = float(x0) 
		self.__c    = float(c)
		self.__gama  = Gamma.gamma_MHz[atm]
	
		self.__v    = array( [0.0, 0.0, self.__c], dtype=float )

		self.__bounds = { 'R1':None, 'R2':None, 'x0':None, c:None }

		self.history= None


	def __str__( self ):
		s  = 'Spin object:\n' 	
		s += '\tGy %f Mhz/T.\n' % (self.__gama)   
		s += '\tx0 %f ppm.\n' %(self.__x0)
		s += '\tR1:  %f\n\tR2:  %f\n\tc:   %f\n' %(self.__R1,self.__R2, self.__c)
		s += '\tv=[ %f, %f, %f]\n' %(self.__v[0],self.__v[1],self.__v[2])
		return s

	
	def set_lb_ub( self, key, lb = None, ub = None ):
		if key not in self.__bounds.key(): raise ValueError
		self.__bounds[key] = ( lb, ub )	

	def get_opt_vals( self ):
		x = []
		for key in sorted(self.__bounds):
			if self.__bounds[key] == None: continue		
			
			elif key == 'R1': x.append(self.__R1) 
			elif key == 'R2': x.append(self.__R2)
			elif key == 'x0': x.append(self.__x0) 
			elif key == 'c' : x.append(self.__c)
			else: raise ValueError 
		return x

	def set_obt_vals( self, vals ):
		i=0
		for key in sorted(self.__bounds):
			if self.__bounds[key] == None: continue
			
			elif   key == 'R1': self.__R1 = vals[i]
			elif   key == 'R2': self.__R2 = vals[i]
			elif   key == 'x0': self.__x0 = vals[i]
			elif   key == 'c' : 
				self.__c    = vals[i]
				self.__v    = array( [0.0, 0.0, self.__c], dtype=float )
			else: raise ValueError
			i+=1

	def get_opt_limits( self ):
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
	print h
	print h.v
