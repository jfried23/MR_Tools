import os
import sys

_file_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append( os.path.split( _file_path )[0] )


import argparse
import h5py
import FileIO.SiemensReader
import MRI.MRI_Viewer.CEST_MRI
import MRI.MRI_Viewer.MRI

def print_attrs(name, obj):
	if isinstance( obj, h5py._hl.group.Group ): return None
	else: print obj.name

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Display Siemens MRI datasets.')
	
	parser.add_argument('-f','--file', nargs=1, required=True, 
                            help='Path(s) to the raw Siemens dataset(s).')
	parser.add_argument('-c', '--cest', default= False, action='store_true',
			    help='This is a cest dataset True/False.')

	parser.add_argument('-r', '--raw', default= False, action='store_true',
			    help='This is time domain data True/False.')



	
	args = parser.parse_args()
	
	if args.raw:
		data = FileIO.SiemensReader.read_siemens( os.path.abspath( args.file[0] ) )
	else:
		f= h5py.File( os.path.abspath( args.file[0] ) )
		print 'Avialable files in this repository:'

		f.visititems(print_attrs)		

		a = raw_input("\nEnter the file name: ")
		data = f[a][...]
	
	if args.cest:
		MRI.MRI_Viewer.CEST_MRI.CEST_MRI( data )
	else:
		MRI.MRI_Viewer.MRI.MRI( data )
