import sys
sys.path.append('../../../')

import util.Gamma
import numpy as np
import scipy as sp
import math

"""
Pulse Elements will not incrament unless/untill self.inc() is expliciatly called.
"""

class PulseElemBase( object ):
	def __init__(self, time_in_s, B1_in_Hz, offset, phase_in_deg, atm='1H' ):
		
		self._t     = np.atleast_1d( time_in_s)      #Durration        in units 's'
		self._B1    = np.atleast_1d( B1_in_Hz )      #Magnitude of B1  in units 'Hz'
		self._o     = np.atleast_1d( offset)         #Offset           in units 'Hz'
		self._phase = np.atleast_1d( phase_in_deg)   #Phase            in units 'rad'    
		self._gamma = util.Gamma.gamma_MHz[ atm ]	
		self._counter=0
		self._npts = max( map(len, [self._t, self._B1, self._o, self._phase]))

		self._M = None

	def M( self, BlochObj ):
		if isinstance(self._M, type(None) ): 
			self._M = sp.linalg.expm( self.time * BlochObj.dM( self ) )

		return self._M

	def inc(self):  
		self._counter += 1
		if self._npts != 1: self._M = None	

	def reset(self): 
		self._counter = 0
		self._M = None
	

	#####The accessors return the value of each paramter for *self._counter* counter value
	@property
	def npts(self): return self._npts

	@property
	def time(self): return float(self._t[ self._counter % len(self._t) ])

	@property
	def phase(self): return float(self._phase[ self._counter % len(self._phase) ])

	@property
	def offset(self): return float(self._o[ self._counter % len(self._o) ])

	@property 
	def gamma(self): return self._gamma

	@property
	def B1(self):
		pwr = self._B1[self._counter %  len(self._B1) ]
		return pwr * np.array([ math.cos( self.phase ), math.sin( self.phase ) ])

class Delay( PulseElemBase ):
	def __init__(self, time ):		
		PulseElemBase.__init__(self, time, 0., 0., 0., '1H')

class Acq( Delay ):
	def __init__(self, time = 0.0 ):
		Delay.__init__(self, time )	

class Pulse( PulseElemBase ):
	def __init__(self, time, angle_in_degrees, phase_in_degrees, offset=0.0, atm='1H' ):
		torad=np.vectorize( math.radians )
		phase      =  torad( np.atleast_1d(phase_in_degrees) )
		B1         =  (angle_in_degrees/360.)/time
		
		PulseElemBase.__init__(self, time, B1, offset, phase, atm)


class CW( PulseElemBase ):
	def __init__(self, time, B1_in_Hz, offset, phase=0, atm='1H'):
		PulseElemBase.__init__(self, time, B1_in_Hz, offset, phase, atm)


class CompositePulse( PulseElemBase ):
	"""
	CompositePulses are a container for multiple Pulse Elements.
	"""
	def __init__( self ):
		self._seq  = []	
		PulseElemBase.__init__(self, 0.0, 0.0, 0.0, 0.0, '1H')

	def inc( self ):
		for p in self._seq : p.inc()
		self._counter += 1
		if self._npts != 1 : self._M = None

	def reset( self ): 
		for p in self._seq : p.reset()	
		self._M = None	

	def M( self, BlochObj ):
		if isinstance(self._M, type(None) ):
			for p in self._seq:
				try: self._M = self._M.dot( p.M( BlochObj ) )
				except: self._M =  p.M( BlochObj )
		return self._M
		
	def add( self, pulse_elem):
		if not isinstance(pulse_elem, PulseElemBase ): 
			raise ValueError("A Non-pulse element was added to the pulseseqence!")

		self._seq.append( pulse_elem )

		self._t     += pulse_elem.time 
		self._npts   = max( self._npts, pulse_elem.npts )

class Loop( CompositePulse ):
	def __init__( self, ncyc ):
		self._seq  = []	
		self._ncyc = ncyc
		PulseElemBase.__init__(self, 0.0, 0.0, 0.0, 0.0, '1H')

	def inc( self ):
		for p in self._seq : p.inc()
		self._counter += 1
		if self._npts != 1 : self._M = None

	def reset( self ): 
		for p in self._seq : p.reset()
		self._M = None	
	
	def M( self, BlochObj ):
		if isinstance(self._M, type(None) ):
			for p in self._seq:
				try: self._M = self._M.dot( p.M( BlochObj ) )
				except: self._M =  p.M( BlochObj )
		self._M =  np.linalg.matrix_power( self._M, self._ncyc )
		return self._M

	@property
	def time(self): return self._t*self._ncyc


	def add( self, pulse_elem):
		if not isinstance(pulse_elem, PulseElemBase ): 
			raise ValueError("A Non-pulse element was added to the pulseseqence!")

		self._seq.append( pulse_elem )

		self._t     += pulse_elem.time 
		self._npts   = max( self._npts, pulse_elem.npts )


if __name__ =='__main__':
	x = np.linspace(1.e-6,1.e-3,101)
	
	p = Pulse( x, 90, 0.0, 0 )
	print p.B1
	p.inc()
	print p.B1

	acq = Acq( x )

