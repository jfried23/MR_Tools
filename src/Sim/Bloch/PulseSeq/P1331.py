import sys
sys.path.append('../../../')

import util.Gamma as Gamma
import Sim.Bloch.PulseSeq.PulseElem as PulseElem

def P1331( ppm, pw90 = 10.e-9, B0=9.4, atm='1H', cnt_pw = False ):
 
	tau = 1./(2. * Gamma.gamma_MHz[ atm ] * B0 * ppm)

	pw1 = pw90/8
	pw3 = 3*pw1

	if cnt_pw: pw3 = pw1
	
	cp = PulseElem.CompositePulse()
	
	cp.add( PulseElem.Pulse( pw1, 90./8, 0) ) 	
	cp.add( PulseElem.Delay( tau ) )
	cp.add( PulseElem.Pulse( pw3, 3*90./8, 180 ) )
	cp.add( PulseElem.Delay( tau ) ) 
	cp.add( PulseElem.Pulse( pw3, 3*90./8, 0 ) ) 
	cp.add( PulseElem.Delay( tau ) ) 
	cp.add( PulseElem.Pulse( pw1, 90./8, 180 ) )

	return cp

if __name__ == '__main__':
	h = P1331(2.0)
	print h
