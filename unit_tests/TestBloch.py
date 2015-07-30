import sys
sys.path.insert(0, '../')
import unittest
from math import e, pi, cos, sin

import numpy as np

import src.Spin  as Spin
import src.util.Gamma as Gamma
import src.Bloch.BlochSim as BlochSim
import src.Bloch.PulseSeq.PulseElem as PulseElem
import src.Bloch.PulseSeq.PulseSeq as PulseSeq


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
"""
	def test_180(self):
		s1=Spin(1.0, 1.0, 2.0, 1.0)
		b=Bloch( s1 )	
		
		ps = PulseSeq( 9.4 )
		ps.add( pulse( 1.e-9, 180., 90) )
		ps.add( acq(0) )

		M = b.run( ps )	
		
		err  =  abs(s1.history[0] - complex( 0, 0 )) + abs( s1.history[1] + 1.0 )

		self.assertAlmostEqual(err, 0, 5)

	def test_R1(self):
		
		test_times = linspace( 0, 2, 20 )

		B0  =  9.4
		r1  =  3.0
		r2  = 50.0
		c   =  2.0
		g  = Gamma.gamma[ '1H' ]

		s1=Spin(r1, r2, 4.0, c)
		b=Bloch( s1 )	
					
		ps = PulseSeq(B0)
		ps.add( pulse( 1.e-9, 180., 0.) )
		ps.add( acq( test_times ) )

		M = b.run( ps )	

		thry = [ c-((c+2.0)*e**(-t*r1)) for t in test_times] 
		
		err = sum(abs( s1.history[1] - thry))		
		self.assertAlmostEqual(err, 0.0, 5)


	def test_R2(self):
		
		test_times = linspace( 0, 80.e-3, 40 )

		B0  =   9.4
		r1  =   0.0
		r2  =  50.0
		c   =   1.0
		x=3.0
		g  = Gamma.gamma[ '1H' ]

		s1=Spin(r1, r2, x, c)
		b=Bloch( s1 )	
			
		ps = PulseSeq(B0)
		ps.add( pulse( 1.e-9, 90., 0.) )
		ps.add( acq( test_times ) )

		M = b.run( ps )	
	
		thry = [c*sin( 2*pi*B0*g*x*t )*e**(-r2*t) for t in test_times ]
		err = sum( real( s1.history[0]) - thry)
	
		self.assertAlmostEqual(err, 0.0, 4)

		

	def test_evolution(self):
		
		test_times = linspace( 0, 100.e-3, 20 )

		B0  =  9.4
		r1  =  0.0
		r2  = 50.0
		ppm =  pi
		c   =  2.0

		g  = Gamma.gamma[ '1H' ]

		s1=Spin(r1, r2, ppm, c)
		b=Bloch( s1 )	
			
		w  = s1.ppm * g * B0


		ps = PulseSeq(B0)
		ps.add( pulse( 1.e-9, 90., 90.) )
		ps.add( acq( test_times ) )

		M = b.run( ps )	

		x_model = [ c*cos( 2*pi*w*t )*(e**(-r2*t)) for t in test_times]
		y_model = [ c*sin( 2*pi*w*t )*(e**(-r2*t)) for t in test_times]		
		
		thr = [ complex( x_model[i], y_model[i] ) for i in range(len( test_times))]
		err = sum(abs(s1.history[0] - thr))/len( test_times )
		
		self.assertAlmostEqual(err, 0.0, 5)


	def test_p1331(self):
		s1=Spin(0.7,1.0, 2.0, 1.0)
		s2=Spin(0.7,2.0, 0.0, 1.0)
		
		b=Bloch(s1,s2)

		tau = 1/(2. * 400 * 2)
		ps = PulseSeq( 9.4 )
		ps.add( pulse( 1.e-9, 90./8, 0) )
		ps.add( delay( tau ) )
		ps.add( pulse( 3.e-9, 3*90./8, 180) )
		ps.add( delay( tau ) )
		ps.add( pulse( 3.e-9, 3*90./8, 0) )
		ps.add( delay( tau ) )
		ps.add( pulse( 1.e-9, 90./8, 180) )
		ps.add( acq(0) )

		M = b.run( ps )
		err = abs( s1.history[0] -complex(0,1) ) + abs(s2.history[0] -complex(0,0))
		err += abs( s1.history[1]-complex(0,0)) + abs(s2.history[1]-complex(1,0))	
		self.assertAlmostEqual( err/4, 0.0, 1)

"""
