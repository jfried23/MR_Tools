import os
import glob
import sys
import time

_file_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append( os.path.split( _file_path )[0] )

import argparse
import FileIO.SiemensReader
import h5py


def _full_path( partial_path ): return os.path.abspath( partial_path ) 
def _parent_dir( full_path ): return os.path.split(full_path)[0]
	


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Processes Siemens MRI datasets.')

	parser.add_argument('-f','--files', nargs='*', required=True, 
                            help='Path(s) to the raw Siemens dataset(s).')


	parser.add_argument('-g','--group', nargs=1, default=None, 
                           help='Name of the group in hdf5 repo. If not supplied no group is used.')


	parser.add_argument('-d','--dataset', nargs=1, default=None,
			    help='Path to an exisiting hdf5 repository.')


	
	args = parser.parse_args()
	
	if args.dataset	== None:
		#First make a defulat name just in case
		args.dataset = time.strftime("%Y%m%d_MRIs.hdf5")

		#Then look for an exisiting hdf5 file
		for filename in os.listdir('./'): 
			if filename.endswith(".hdf5"): 	args.dataset = filename
	else: args.dataset = args.dataset[0]

	

	hdf5 = h5py.File( _full_path(args.dataset) )
			 
	for filename in args.files:
		try:
			filepath = _full_path( filename )
			filename = os.path.split(filepath)[1]

			if args.group != None: filename = os.path.join( args.group[0], filename)

			data = FileIO.SiemensReader.read_siemens( filepath )

			hdf5.create_dataset( filename, data=data, shuffle=True, 
                                          	chunks=True, compression='gzip' )
		except: 
			print "File %s already exists! Skipping.." % ( filename )

	hdf5.flush()
	hdf5.close()



	
