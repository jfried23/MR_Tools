import sys
sys.path.append('../../')

import numpy as np
from scipy import optimize


import Lorentzian
import Sim.Spin


class LorenDiff( object ):
	
	def __init__(self, x_ppm, y_vals ):
		self.y_cest = y_vals

		self.water = Sim.Spin.Spin( R1=0, R2=20, x0=0.0, c = 100 )
		self.water.set_lb_ub( 'R2', None, None )
		self.water.set_lb_ub( 'c',  None, None )
		self.water.set_lb_ub( 'x0', None, None )

		self.baseline = y_vals[-1]
		self.lorenz = Lorentzian.Lorentzian( self.water, B0 = 7.0, ppm=x_ppm )

		
		x0 = self.water.get_opt_vals()
		x0.append( self.baseline )
		
		#Setup complete fire up the optimizer!
		opt = optimize.minimize( self.calc_rmsd, x0, method='Nelder-Mead', options={'maxiter':1e4,'maxfev':1e4} )	

		self.result = self.calc_diff()		
	
	def calc_diff( self ): return self.y_cest - (self.lorenz.run() + self.baseline )

	def calc_rmsd( self, x0 ):
		x0=np.ndarray.tolist( x0 )
		self.water.set_opt_vals( x0 )
		self.baseline = x0[-1] 	
	
		return np.linalg.norm( self.y_cest - (self.lorenz.run() + self.baseline ) )	

	
if __name__ == '__main__':
	import matplotlib.pyplot as plt 
	import numpy as np

	y_cest = np.load('/Users/josh/Documents/Programs/MR_Tools/cest_data.p' )[2]
	x_ppm  = np.linspace( -8.,8.,51 )

	 
	a= LorenDiff( x_ppm, y_cest )
	print a.water
	plt.plot(x_ppm, y_cest, 'bo' )
	plt.plot( x_ppm, a.lorenz.run() + a.baseline )

	plt.plot( x_ppm, a.result,'r' )
	plt.show()
