import sys
sys.path.append('../../')

import numpy as np
from scipy import optimize

import Shapes

class SimLineshape( object ):

	def __init__( self, *arg, **kwargs):

		self.B0      = 7.0	
		self.spins   = arg
	
		for k,v in kwargs.items(): 
			if k == 'B0': self.B0 = v
			else: raise ValueError
	
	def run( self, ppm=np.linspace(-8,8,1000) ):

		ls=np.zeros_like( ppm) 
		for s in self.spins:
			ls += s.lineshape( ppm, s, B0 = self.B0 )

		return ls	
			
	def fit( self, ppm, data, solver = 'L-BFGS-B' ):
		x0=[]
		bounds = []
		for s in self.spins:
			x0.extend( s.get_opt_vals() )
			bounds.extend( s.get_opt_limits() )
		
		opt = optimize.minimize( self._err_func, x0, 
                                          args=(ppm, data), 
					  bounds = bounds, 
                                          method= solver, jac=False, 
                                          options={'maxiter':5e3,'maxfev':1e7} 
				       )
		return opt
		#x0 = np.ndarray.tolist( opt.x )
		#for one in self.spins: s.set_opt_vals( x0 )

	
	def _err_func( self, x0, ppm, data ):
		x = np.ndarray.tolist( x0 )
		for s in self.spins: 
			s.set_opt_vals( x )
		return np.linalg.norm( data - self.run( ppm )  )


if __name__ == '__main__':
	import Sim.Spin
	import matplotlib.pyplot as plt


	d = np.load('/Users/josh/Desktop/t.p.npy' )
	frq = d[0]
	cst = d[1]
	

	s1 = Sim.Spin.Spin( R2=1.0e3,  x0=0., c=-500)
	s2 = Sim.Spin.Spin( R1=0.97, R2=0)
	
	
	s1.lineshape = Shapes.Lorentzian
	s2.lineshape = Shapes.linear_baseline
	

	s1.set_lb_ub( 'R2', 1, None )
	s1.set_lb_ub( 'x0', -1, 1 )
	s1.set_lb_ub( 'c', None, 0 )

	s2.set_lb_ub( 'R1', 0, 2 )

	s = SimLineshape( s1, s2  )

	#plt.plot(s.run(),'r')

	s.fit( frq/300, cst) #, solver='Nelder-Mead'return  print s1
	print s2

	plt.plot( frq/300,cst,'o')
	plt.plot( np.linspace(-8,8, 1000), s.run(),'b' ); plt.show()
