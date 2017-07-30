import sys
import os.path
import itertools
import pandas as pd
import numpy as np

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
		print('counting ERV...')
		moves_start = list(range(0,-nmer,-1))
		moves_skip = list(range(0,nmer))
        
		with open(filename) as f:
			for line in itertools.islice(f, 1, None):
				s , m = str.split(line)[4:6]
				for i in range(len(self.data)): 
					temp = s
					s = s[center+moves_start[i]:center+moves_start[i]+nmer] # skip center@@
					s = s[:moves_skip[i]]+s[moves_skip[i]+1:]
					# print('i={},s={},m={}'.format(i,s,m))
					self.data[i]['nERVs'][self.mtypes.index(m)*6+self.subtypes.index(s)+1]+=1
					s = temp

        
	def __init__(self, nmer, ervfile, reffile, center):
		if is_int(nmer) == False:
			raise ValueError('nmer should be integer greater than 1')
		nmer = int(nmer)
		center = int(center)

		self.patterns = list([''.join(i) for i in itertools.permutations('X'*(nmer-1)+'*')])
		self.mtypes = ['AT_CG', 'AT_GC', 'AT_TA', 'GC_AT', 'GC_CG', 'GC_TA']
		self.subtypes = [''.join(i) for i in itertools.product('ACGT', repeat = (nmer-1))]
		self.data = []
		for i in range(0, nmer):
			#####################################
			## change column numbers with argc ###
			#####################################
			self.data.append(pd.DataFrame(np.zeros((4 ** (nmer-1) * 6, 5),dtype=np.int32),
	                                    columns=['pattern', 'mtype', 'subtype', 'nERVs', 'nMotifs']))                                          
			self.data[i]['mtype'] = list(itertools.chain.from_iterable(itertools.repeat(x,4 ** (nmer-1)) for x in self.mtypes))
			self.data[i]['subtype'] = self.subtypes * 6
			self.data[i]['pattern'] = self.patterns[i]
	    
		if center is None:
			center = 3

		if ervfile is not None:
			if os.path.isfile(ervfile)==False:
				raise ValueError('{} is not a file'.format(ervfile))
	            
			self._countERV(ervfile, nmer, center)
			print('counting ERV completed')
		else:
			print('erv not counted as ervfile is None')
	    
		if reffile is not None:
			#####################################
			#### count rel rate from relfile ####
			#####################################
			pass
		else:
			print('reference motifs not counted as reffile is None')
	        
		if ((ervfile is not None) & (reffile is not None)):
			print('counting relrate and wt ...')
			total_motifs = np.sum(self.data[0].nMotifs)
			for i in range(0, nmer):
				self.data[i]['ERV_rel_rate'] = self.data[i].nERVs / self.data[i].nMotifs
				self.data[i]['wt'] = self.data[0].nERVs / total_motifs
			print('counting relrate and wt completed')

	def writeERV(self, dir):
		print('writing data to {}...'.format(dir))
		if os.path.isdir(dir)==False:
			raise ValueError('{} is not a directory'.format(dir))
		if dir.endswith('/')==False:
			dir = dir+'/'

		out = pd.concat(self.data)
		out.to_csv(dir+'{}mer_.txt'.format(len(self.data)), sep=' ', index=False, header=True)
		print('writing data to {} complete'.format(dir))
        