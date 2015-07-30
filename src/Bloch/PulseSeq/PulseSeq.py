import sys
sys.path.append('../../../')


import PulseElem

class PulseSeq( object ):
	def __init__( self ):
		self.__inc  = 0
		self.__ncyc = 0
		self.__seq = []		

	@property
	def ncyc( self ): return ( self.__ncyc)

	def __iter__(self):
		while self.__inc  < self.__ncyc:
			self.__inc += 1 
			for p in self.__seq: 
				yield p
				p.inc()
			yield PulseElem.Acq()
		
		self.reset()

	def reset( self ):
		for p in self.__seq: p.reset()
		self.__inc  = 0


	def add( self, pulse ):
		if isinstance(pulse, PulseElem.PulseElemBase ):  self.__seq.append( pulse )
		else: raise ValueError("%s object cannot be added to PulseSeq!" %(type(pulse)))

		self.__ncyc = max( self.__ncyc, pulse.npts ) 	
