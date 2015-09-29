import os
import sys

_file_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append( os.path.split( _file_path )[0] )


import argparse
import FileIO.SiemensReader
import MRI.MRI_Viewer.CEST_MRI
import MRI.MRI_Viewer.MRI


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Display Siemens MRI datasets.')
	
	parser.add_argument('-f','--file', nargs=1, required=True, 
                            help='Path(s) to the raw Siemens dataset(s).')
	parser.add_argument('-c', '--cest', default= False, action='store_true',
			    help='This is a cest dataset True/False.')


	
	args = parser.parse_args()
	
	data = FileIO.SiemensReader.read_siemens( os.path.abspath( args.file[0] ) )

	
	if args.cest:
		MRI.MRI_Viewer.CEST_MRI.CEST_MRI( data )
	else:
		MRI.MRI_Viewer.MRI.MRI( data )
