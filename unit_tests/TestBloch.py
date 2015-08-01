import sys
sys.path.insert(0, '../src/')
import unittest
from math import e, pi, cos, sin

import numpy as np

import Sim.Spin  as Spin
import util.Gamma as Gamma
import Sim.Bloch.BlochSim as BlochSim
import Sim.Bloch.PulseSeq.PulseElem as PulseElem
import Sim.Bloch.PulseSeq.PulseSeq as PulseSeq


class Test_Bloch(unittest.TestCase):

	def setUp(self): 
		self.cplx = np.vectorize( complex )

		self.g = Gamma.gamma_MHz[ '1H' ]



	def test_90(self):
		s1 = Spin.Spin(1.0, 1.0, 2.0, 1.0)
		b  = BlochSim.BlochSim( s1 )	
		ps = PulseSeq.PulseSeq()
		ps.add( PulseElem.Pulse( 1.e-9, 90., 0.) )

		v = b.run( ps )	
		err = np.linalg.norm( v.T - np.array([ 0.0, -1.0, 0.0, 1.0 ] ) )
		self.assertAlmostEqual(err, 0, 5)



	def test_R1(self):
		
		test_times = np.linspace( 0, 4, 10 )

		B0  =  9.4
		r1  =  1.0
		r2  =  0.0
		c   =  2.0

		s1=Spin.Spin(r1, r2, 4.0, c)

		b=BlochSim.BlochSim( s1, B0=9.5 )	
					
		ps = PulseSeq.PulseSeq()
		ps.add( PulseElem.Pulse( 1.e-9, 180., 0.) )
		ps.add( PulseElem.Delay( test_times ) )

		M = (b.run( ps ))[:,2]
		thry = [ c-((c+2.0)*e**(-t*r1)) for t in test_times] 
		
		err = np.linalg.norm( M - thry )

		self.assertAlmostEqual(err, 0.0, 6)



	def test_R2( self ):
		test_times = np.linspace( 0, 1, 10 )
		
		B0  =  9.4
		r1  =  0.0
		r2  =  5.0
		c   =  2.0

		s1=Spin.Spin(r1, r2, 0.0, c)
		b=BlochSim.BlochSim( s1, B0=9.5 )

		ps = PulseSeq.PulseSeq()
		ps.add( PulseElem.Pulse( 1.e-9, 90., 90.) )
		ps.add( PulseElem.Delay( test_times ) )

		M = (b.run( ps ))[:,0]

		thry = [ c*e**(-t*r2) for t in test_times] 

		err = np.linalg.norm( M - thry )

		self.assertAlmostEqual(err, 0.0, 6)



	def test_evolution( self ):

		test_times = np.linspace( 0, 100.e-3, 20 )

		B0  =  9.4
		r1  =  0.0
		r2  =  5.0
		ppm =  pi
		c   =  2.0
		
		w = B0*self.g*ppm

		s1=Spin.Spin(r1, r2, ppm, c)
		b=BlochSim.BlochSim( s1, B0=B0 )

		ps = PulseSeq.PulseSeq()
		ps.add( PulseElem.Pulse( 1.e-9, 90., 90.) )
		ps.add( PulseElem.Delay( test_times ) )

		M =  b.run( ps )[:,0:2]	
	
		M = self.cplx( M[:,0], M[:,1] )

		x = [ c*cos( 2*pi*w*t )*(e**(-r2*t)) for t in test_times]
		y = [ c*sin( 2*pi*w*t )*(e**(-r2*t)) for t in test_times]		
		
		thry = self.cplx( x, y  )
		
		err = np.linalg.norm( M - thry )
		
		self.assertAlmostEqual(err, 0.0, 3)


	def test_P1331( self ):
		
		ppm = 2.0
		B0=9.4

		tau = 1./(2. * self.g*B0 * ppm)


		s1=Spin.Spin(1.7,1.0, ppm, 1.0)
		s2=Spin.Spin(1.7,2.0, 0.0, 1.0)
		
		b=BlochSim.BlochSim( s1,s2, B0 = B0 )		

		ps = PulseSeq.PulseSeq()

		ps.add( PulseElem.Pulse( 10.e-9, 90./8, 0) )
		ps.add( PulseElem.Delay( tau ) )
		ps.add( PulseElem.Pulse( 30.e-9, 3*90./8, 180 ) )
		ps.add( PulseElem.Delay( tau ) )
		ps.add( PulseElem.Pulse( 30.e-9, 3*90./8, 0 ) )
		ps.add( PulseElem.Delay( tau ) )
		ps.add( PulseElem.Pulse( 10.e-9, 90./8, 180) )
		
		M = b.run(ps)
		err = np.linalg.norm( M.T - np.array([ 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0 ] ) )
	
		self.assertAlmostEqual(err, 0.0, 2)


if __name__ == '__main__':
	unittest.main()
