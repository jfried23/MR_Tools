import sys
sys.path.append('../../')

import Sim.Spin

import numpy as np
from math import pi

class Lorentzian( object ):

	def __init__( self, *arg  ):
		self.spins = arg
		self.B0 = 9.4
		self.Rmtx    = np.zeros( (len(arg),len(arg) ) )

	def add_kex(self, frm, to, kex ):	

		frm, to = frm-1, to-1		
		
		cfct =  (self.spins[frm]).c/(self.spins[to]).c 
		
		self.Rmtx[frm,to]= kex/( 1. + cfct )
		self.Rmtx[to,frm]= kex - self.Rmtx[frm,to]

		self.Rmtx[frm,frm] = 0.0
		self.Rmtx[to, to ] = 0.0	
		self.Rmtx[frm,frm] = np.sum( self.Rmtx[frm] )
		self.Rmtx[to, to ] = np.sum( self.Rmtx[to ] )		

	def run( self, ppm_range ):
		spec = np.zeros_like( ppm_range, dtype='float')
		for i,s in enumerate(self.spins):
			R2 = s2.R2 + self.Rmtx[i,i]
			spec +=  (1./pi)*(s.v[-1]*R2/( (R2**2) +((x-s.x0)*s.gamma*self.B0)**2))
		return spec

if __name__ == '__main__':
	import Sim.Spin
	import matplotlib.pyplot as plt

	s1 = Sim.Spin.Spin( R1 = 1.0, R2 = 50, x0=4, c=-50 )
	s2 = Sim.Spin.Spin( R1 = 1.0, R2 = 200, x0=-4, c=-1 )

	l = Lorentzian( s1, s2 )
	l.add_kex(1,2,0)
	x = np.linspace( -8,8,10000 )


	plt.plot( x, l.run(x) )
	plt.show()

