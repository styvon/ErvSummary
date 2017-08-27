import sys
import os.path
import itertools
import pandas as pd
import numpy as np
import re
from time import sleep

pd.options.mode.chained_assignment = None

def is_int(s):
	try:
		val = int(s)
		if val < 2:
			return False
		else:
			return True
	except ValueError:
		return False

class ErvSummary:
	"""summary for nmer motifs
	requires pandas, numpy"""
	def _countERV(self, filename, nmer, center):
		"""count nERV from target file,
		center: center base position in the summary data (starting from 0)"""
		if not os.path.isfile(filename):
			raise ValueError('{} is not a valid directory'.format(filename))
		moves_start = list(range(0,-nmer,-1))
		moves_skip = list(range(0,nmer))
		steplen = (4**(nmer-1))
		
		with open(filename) as f:
			for line in itertools.islice(f, 1, None):

				s , m = str.split(line)[6:8]
				m = m[-5:]
				### if mutation type not in list, skip???###
				if m not in self.mtypes:
					print("Found an unknown mutation type", m)
					continue
				
				######################################

				for i in range(len(self.data)): 
					temp = s
					s = s[center+moves_start[i]:center+moves_start[i]+nmer] 
					s = s[:moves_skip[i]]+s[moves_skip[i]+1:] # skip center
					# print('i={},s={},m={}'.format(i,s,m))
					### if subtype not in list, skip???###
					if s not in self.subtypes:
						print("Found an unknown subtype", s, "in pattern", self.patterns[i])
						s = temp
						continue
					######################################
					self.data[i]['nERVs'][self.mtypes.index(m)*steplen+self.subtypes.index(s)]+=1
					s = temp

	def _countMOT(self, filename, nmer, center):
		"""count nMotifs from target file,
		center: center base position in the summary data (starting from 0)"""
		if not os.path.isfile(filename):
			raise ValueError('{} is not a valid directory'.format(filename))
		moves_start = list(range(0,-nmer,-1))
		moves_skip = list(range(0,nmer))
		
		#with open(filename,encoding='latin-1') as f:
		with open(filename) as f:
			for line in itertools.islice(f, 1, None):
				line = re.sub("\x08.", "", line)
				s, n  = line.split("\t")[1:3]
				n = int(float(n))

				for i in range(len(self.data)): 
					temp = s
					s = s[center+moves_start[i]:center+moves_start[i]+nmer] 
					s_type = s[:moves_skip[i]]+s[moves_skip[i]+1:] # strip subtype by skipping center posi
					s_cent = s[moves_skip[i]] # strip "center" site
					### if subtype not in list, skip current scan###
					if s_type not in self.subtypes:
						print("Found an unknown subtype", s, "in pattern", self.patterns[i])
						s = temp
						continue
					######################################
					# print('i={},s={},m={}'.format(i,s,m))
					self.data[i].loc[ (self.data[i].subtype == s_type) & (self.data[i].mtype.str[:1] == s_cent) , 'nMotifs'] += n
					s = temp
		
	def __init__(self, nmer, ervdir, refdir, center):
		if is_int(nmer) == False:
			raise ValueError('nmer should be integer greater than 1')
		else:
			nmer = int(nmer)
		center = int(center)
		ncol = 5

		self.patterns = list(set([''.join(i) for i in itertools.permutations('X'*(nmer-1)+'*')]))
		self.patterns.sort()
		self.mtypes = ['AT_CG', 'AT_GC', 'AT_TA', 'GC_AT', 'GC_CG', 'GC_TA']
		self.subtypes = [''.join(i) for i in itertools.product('ACGT', repeat = (nmer-1))]
		self.data = []
		for i in range(0, nmer):
			
			self.data.append(pd.DataFrame(np.zeros((4 ** (nmer-1) * 6, ncol),dtype=np.int32),
										columns=['pattern', 'mtype', 'subtype', 'nERVs', 'nMotifs']))										  
			self.data[i]['mtype'] = list(itertools.chain.from_iterable(itertools.repeat(x,4 ** (nmer-1)) for x in self.mtypes))
			self.data[i]['subtype'] = self.subtypes * 6
			self.data[i]['pattern'] = self.patterns[i]

		if ervdir is not None: 
			print('counting ERV...')
			if not ervdir.endswith('/'):
				ervdir = ervdir+'/'

			if not os.path.isdir(os.getcwd()+'//'+ervdir): # check valid directory for erv files
				raise ValueError('{} is not a valid directory'.format(ervdir))

			for ervfile in os.listdir(os.getcwd()+'//'+ervdir):
				self._countERV(os.getcwd()+'//'+ervdir + ervfile, nmer, center)
			print('counting ERV completed')
		else:
			print('erv not counted as ervdir is None')
		
		if refdir is not None: 
			print('counting motifs...')
			if not os.path.isdir(os.getcwd()+'/'+refdir): # check valid directory for ref motif files
				raise ValueError('{} is not a valid directory'.format(refdir))

			if not refdir.endswith('/'):
				refdir = refdir+'/'

			
			for reffile in os.listdir(os.getcwd()+'//'+refdir):	
				self._countMOT(os.getcwd()+'/'+refdir + reffile, nmer, center)

			print('counting motifs completed')
		else:
			print('reference motifs not counted as refdir is None')
			
		if ((ervdir is not None) & (refdir is not None)):
			print('counting relrate and wt ...')
			total_motifs = np.sum(self.data[0].nMotifs)
			for i in range(0, nmer):
				self.data[i]['ERV_rel_rate'] = self.data[i].nERVs / self.data[i].nMotifs
				# self.data[i]['wt'] = self.data[0].nERVs / total_motifs
			print('counting relrate and wt completed')

	def writeERV(self, dir):
		print('writing data to {}...'.format(dir))
		if os.path.isdir(dir)==False:
			raise ValueError('{} is not a directory'.format(dir))
		if dir.endswith('/')==False:
			dir = dir+'//'

		out = pd.concat(self.data)
		out.to_csv(dir+'{}mer.txt'.format(len(self.data)), sep=' ', index=False, header=True)
		print('writing data to {} complete'.format(dir))
		