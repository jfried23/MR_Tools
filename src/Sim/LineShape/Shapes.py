import sys
sys.path.append('../../')

import numpy as np
import scipy as sp

def Lorentzian( ppm, spinObj, B0=7.0 ):
	Delta = (ppm-spinObj.x0) * spinObj.gamma * B0
	Rho_inv = np.pi/spinObj.R2
	return spinObj.c*(2/(np.pi*np.sqrt(3)))*Rho_inv / ( 1 + (4./3)*( Delta*Rho_inv)**2)

def Gaussian( ppm, spinObj, B0=7.0 ):
	Delta = (ppm-spinObj.x0) * spinObj.gamma * B0
	Rho_inv = np.pi/spinObj.R2
	return  spinObj.c*(Rho_inv * np.sqrt( 2/np.pi) ) * np.exp( -2* (Delta * Rho_inv )**2 ) 

def Mixture( ppm, spinObj, B0=7.0, frac_Gaus = 0.5 ):
	frac_Loren = 1.0 - frac_Gaus
	
	shp  = frac_Loren * Lorentzian( ppm, spinObj, B0 )
	shp += frac_Gaus  * Gaussian(   ppm, spinObj, B0 )

	return shp

def SuperLorentzian( ppm, spinObj, B0=7.0 ):
	Delta = (ppm-spinObj.x0) * spinObj.gamma * B0
	rslt = np.empty_like( ppm )

	for i,v in enumerate( Delta ):
		rslt[i] = sp.integrate.quad( _eval_Sup, 0., np.pi/2, args=(v,1./spinObj.R2) )[0]
	return spinObj.c*rslt

def _eval_Sup( theta, delta, T2 ):
	fac = T2/np.abs( 3 * (np.cos( theta )**2) - 1 ) 
	e   = np.exp( -2 * ( 2*np.pi*delta*fac)**2 )
	return np.sqrt( 2/ np.pi ) * np.sin( theta ) * (fac) * e 


def logistic_basline( ppm, spinObj, B0=7.0 ): 
	return (spinObj.c/( 1 + np.exp( -spinObj.R2 * ( ppm - spinObj.x0 )* spinObj.gamma * B0 )) )+spinObj.R1 

def linear_baseline( ppm, spinObj, B0=7.0 ): 
	return spinObj.R1 + ppm*spinObj.R2*spinObj.gamma * B0 
	





if __name__ == '__main__':
	import Sim.Spin
	import matplotlib.pyplot as plt

	ppm = np.linspace(-80,80,3000 )

	s1 = Sim.Spin.Spin( R1 = 1.0, R2 = 200, x0=4, c=50 )
	s2 = Sim.Spin.Spin( R1 = 1.0, R2 = 200, x0=-4, c=50)

	from scipy import integrate
	
	#print np.trapz(Lorentzian(ppm,s1),ppm )
	#print np.trapz(Gaussian(ppm,s1),ppm )
	#print np.trapz(Mixture(ppm,s1),ppm )

	plt.plot(ppm, Lorentzian(ppm,s1) )
	plt.plot(ppm, Gaussian(ppm, s1), 'r' )
	plt.plot(ppm, Mixture(ppm, s1), 'g' )

	plt.plot(ppm, SuperLorentzian(ppm, s1), 'k' )

	plt.show()

