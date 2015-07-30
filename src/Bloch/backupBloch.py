
from sim.Spin import *
from PulseSeq import *

from math import pi
from numpy import zeros, ones, dot,linspace, shape, transpose, squeeze, zeros_like
from scipy.linalg import expm


class Bloch:
	"""
	Members:
		spins: an array of Spin objects

		Rmtx: an 'nxn' matrix whoes elements (i,j) describe the the magnitizaion
		      exchnage rate between spin 'i' and spin 'j'. At object creation all		       		 elements are initalized to zero.
	Methods:
		add_kex( i,j, kex ): establish magnitization transfer is occuring between
			spins[i] -> spins[j] @ rate kex.

		dM( B0, pulse ):  returns the differential matrix descrbing evolution 				      of the system.
				     B0 static magnetic field strength in Tesla
				     pulse the PulseElem object
		run( pulse_seq, step=0.0):     runs the Bloch simulation
				             takes a pulse sequence object
						Returns the magnitization vector in the following
						format
							
						M[ inc, Si, ti ]
						Si = [S1_x        ti =[ 0.0
						      S1_y              1.e-9
						      S1_z 		2.e-9		
				M[0,:,:] =		.		 .
							.		 .
							.                .
						      SN_z           final_time
							1 ]	          ]
						
						
	"""

	def __init__(self, *arg):
		self.__spins = arg
		self.Rmtx = zeros( (len(arg),len(arg) ) )
		
		self.nspins = len( self.__spins )	

	def add_kex(self, frm, to, kex ):
		frm, to = frm-1, to-1		
		
		cfct =  (self.__spins[frm]).c/(self.__spins[to]).c 
		
		self.Rmtx[frm,to]= kex/( 1. + cfct )
		self.Rmtx[to,frm]= kex - self.Rmtx[frm,to]

		self.Rmtx[frm,frm] = 0.0
		self.Rmtx[to, to ] = 0.0	
		self.Rmtx[frm,frm] = 2*pi*sum( self.Rmtx[frm] )
		self.Rmtx[to, to ] = 2*pi*sum( self.Rmtx[to ] )
		
	def dM( self, B0, pulse ):
		sz = self.nspins
			
		dM  = zeros( ( 3*sz+1, 3*sz+1) ) #the transformation matrix
		
		for i in range(sz):  
			x = (3*(i));
                	y = (3*(i))+1;
                	z = (3*(i))+2;

			g  = (self.__spins[i]).gamma
			R1 = (self.__spins[i]).R1
			R2 = (self.__spins[i]).R2
			w  = (self.__spins[i]).ppm * g * B0 
			c  = (self.__spins[i]).c
			
			if pulse.atm()  == (self.__spins[i]).atm:
				B1 = pulse.B1()*2*pi
				dw = pulse.dw()
			else:
				B1 = [0.0, 0.0]
				dw = 0.0
				
			# Terms describing spin under influence of B0 & B1
			dM[x,x] = -( R2 + self.Rmtx[i,i] )
                	dM[x,y] = -2*pi*( w - dw )
                	dM[x,z] =  B1[1]
                
                	dM[y,y] = -( R2 + self.Rmtx[i,i] )
                	dM[y,x] =  2*pi*( w - dw )
                	dM[y,z] = -B1[0]
                
                	dM[z,z] = -( R1  + self.Rmtx[i,i] )
                	dM[z,x] = -B1[1]
                	dM[z,y] =  B1[0]

			dM[z,-1]=   c*R1 


			#Terms describing kex between spins
			for ii in range(sz):
				if ii == i: continue
				
				ox = (3*ii)
				oy = (3*ii)+1
				oz = (3*ii)+2
				
				dM[ox, x] = self.Rmtx[ii,i]
				dM[oy, y] = self.Rmtx[ii,i]
				dM[oz, z] = self.Rmtx[ii,i]
		return dM

	def run( self, pulse_seq ):
		#initalize M0
		
		M0 = ones( ( (3*self.nspins + 1), 1 ) ) 
		for i, s in enumerate(self.__spins): M0[ 3*i : (3*i)+ 3 ] = transpose([s.vec])
		#now find size of collection matrix
		if not isinstance( pulse_seq[-1], acq ):
			print 'Pulse sequence must terminate with an acq command!'
			raise SyntaxError

		ni     = pulse_seq.numcyc()
		vec_sz = len( M0 )
		time   = (pulse_seq[-1])._time
		np = len(time)
		record = zeros( (ni, vec_sz, np) )		
		
		for p in pulse_seq:
			dM = self.dM( pulse_seq.B0(), p )

			if isinstance( p, acq ):
				self._times = p._time
				for i, t in enumerate(p._time):
					if t == 0.0:
						record[ pulse_seq.thiscyc(),:, i] = transpose(M0)
					else:	
						r = transpose(dot(expm( t*dM ),  M0))
						record[ pulse_seq.thiscyc(),:, i] = r 
	
			else:	M0 = dot(expm( p.time()*dM ),  M0)
		

		
		for i, s in enumerate( self.__spins):
			rec_sig = [zeros( (ni,np), complex), zeros( (ni,np) )]
			for inc in range(ni):
				rec_sig[0][inc,:] = squeeze(array((record[ inc, (3*i)  , 0:np]) + 1j*(record[ inc, (3*i)+1, 0:np])))
				rec_sig[1][inc,:] = squeeze(array(record[ inc, (3*i)+2, 0:np]))
			self.__spins[i].history = squeeze(rec_sig)


		if ni == 1 and np > 1:
			#print 'Ni == 1!'
			pass

		elif ni > 1 and np == 1:
			print 'This is prob a z-spectra'

		elif ni == 1 and np ==1:
			return M0	
		
			
if __name__ == '__main__':
	t = linspace( 0, 60.e-3, 1000 )
	g  = Gamma.gamma[ '1H' ]

	s1=Spin(1.0, 5.0,   3.0, 1.0)
	s2=Spin(1.0, 15.0,  1.0, 1.0)
		
	b=Bloch(s1,s2)
	b.add_kex(1,2, 150 )	

	ps = PulseSeq( 9.4 )
	ps.add( pulse( 1.e-9, 90., 90) )
	ps.add( delay( 0.0 ) )
	ps.add( acq(t) )

	b.run( ps )

	from pylab import plot, show
	from numpy import fft, array, real, shape
	plot( s1.history[0] + s2.history[0] )	
	show()
		#print abs(s1.history[0]) + s1.history[1]
		#print abs(s2.history[0]) + s2.history[1]

		#print  s1.history[1]
		#print  s2.history[1]
	
	
