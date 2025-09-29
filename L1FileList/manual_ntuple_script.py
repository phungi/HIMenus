import glob
import os
import time
import sys
import argparse


parser = argparse.ArgumentParser( description = 'Search in eos dir and subdirs for files on to use for the rate-estimation tool' )
parser.add_argument( 'input_directory', help = 'Input directory' )
args = parser.parse_args()
in_dir = args.input_directory

filelist = []
for path, subdirs, files in os.walk(in_dir):
	for name in files:
		if ".root" in name:
			filelist.append(os.path.join(path, name))
for in_file in filelist:
	print('root://eoscms.cern.ch/'+in_file)