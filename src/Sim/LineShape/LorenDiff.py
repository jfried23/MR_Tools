import sys
sys.path.append('../../')

import numpy as np
from scipy import optimize


import Shapes 
import Sim.Spin

def _calc_RMSD( x0, y, ppm, spins ):

	x0 = np.ndarray.tolist(x0)
	ls = np.zeros_like(ppm)

	for i,s in enumerate(spins): 
		s.set_opt_vals( x0 )
		if i == 0: ls += Shapes.logistic_basline(ppm, s ) 
		else:      ls += Shapes.Gaussian(ppm, s) 
	return np.linalg.norm( y - ls  )


def LorenDiff( ppm, y ):

	baseline  = Sim.Spin.Spin( R1=y[-1], R2= 0,  x0=0.0, c = 1.0 )
	water     = Sim.Spin.Spin( R1=0    , R2=100 , x0=0.0, c = 1e3 )
	s1        = Sim.Spin.Spin( R1=0    , R2=1e4 , x0=-3.0, c = 1e3 )


	baseline.set_lb_ub( 'c', None, None )
	baseline.set_lb_ub( 'x0', None, None )
	baseline.set_lb_ub( 'R1', None, None )
	baseline.set_lb_ub( 'R2', None, None )

	water.set_lb_ub( 'R2', None, None )
	water.set_lb_ub( 'c',  None, None )
	water.set_lb_ub( 'x0', None, None )

	s1.set_lb_ub( 'R2', None, None )
	s1.set_lb_ub( 'c',  None, None )
	s1.set_lb_ub( 'x0',  None, None )


	spins = [baseline, water, s1]
	
	x0=[]
	for s in spins: 
		for v in s.get_opt_vals(): x0.append(v)

	opt = optimize.minimize( _calc_RMSD, x0, args=(y, ppm, spins), method='Nelder-Mead', jac=False, options={'maxiter':1e4,'maxfev':1e4} )
	
	ls = np.zeros_like(ppm)
	xxx = np.ndarray.tolist( opt.x )	
	for i,s in enumerate(spins): 
		s.set_opt_vals( xxx )
		if i == 0: ls += Shapes.logistic_basline(ppm, s )
		else: ls +=  Shapes.Guassian(ppm, s )
	return ls 


	
if __name__ == '__main__':
	import matplotlib.pyplot as plt 
	import numpy as np
	
	import h5py
	
	f=h5py.File('/Users/josh/Documents/Data/Tumor_dataset/20150930_MRIs.hdf5')

	d=f['meas_MID167_Gre_2D_Gauss_10X1440d_200FOV_2500Hz_FID4497.dat'][...]

	x_ppm  = np.linspace( -8.,8.,51 )

	ans = np.zeros_like( d )
	for i in range( d.shape[-1] ):
		for ii in range( d.shape[-2] ):
			ans[:,i,ii] =  d[:,i,ii] - LorenDiff( x_ppm, d[:,i,ii ] )
	

	#y_cest = np.load('/Users/josh/Documents/Programs/MR_Tools/cest_data.p' )[2]
	
	 
	#LorenDiff( x_ppm, y_cest )
	#print a.water
	#plt.plot(x_ppm, y_cest, 'bo' )
	#plt.plot( x_ppm, LorenDiff( x_ppm, y_cest ) )

	#plt.plot( x_ppm, a.result,'r' )
	#plt.show()
