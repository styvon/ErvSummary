
import os
os.environ['TF_CPP_MIN_LOG_LEVEL']='2'
import time
import sys
import getopt

import numpy as np
import pandas as pd

import ErvSummary

def main(argv):
	ervdir = None
	refdir = None
	outputdir = None
	nmer = 3
	center = 3

	try:
		opts, args = getopt.getopt(argv,"c:e:h:n:o:r:",["ervdir=","refdir=","outputdir=","nmer=","center="])
	except getopt.GetoptError:
		print('main.py -e <ervdir> -r <refdir> -o <outputdir> -n <nmer> -c <center>')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print('main.py -e <ervdir> -r <refdir> -o <outputdir> -n <nmer> -c <center>')
			sys.exit()
		elif opt in ("-e", "--ervdir"):
			ervdir = arg
		elif opt in ("-r", "--refdir"):
			refdir = arg
		elif opt in ("-o", "--outputdir"):
			outputdir = arg
		elif opt in ("-n", "--nmer"):
			nmer = arg
		elif opt in ("-c", "--center"):
			center = arg
	
	ervsummary = ErvSummary.ErvSummary(nmer=nmer, ervdir=ervdir, refdir=refdir, center=center)
	ervsummary.writeERV(outputdir)
	
if __name__ == "__main__":
	main(sys.argv[1:])


