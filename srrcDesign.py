# Matt @ WaveWalkerDSP.com
# https://wavewalkerdsp.com
# Twitter: @WaveWalkerDSP
#
#
# Copyright 2021, Matt @ WaveWalkerDSP.com
#
# Permission is hereby granted, free of charge, to any person obtaining a copy 
# of this software and associated documentation files (the "Software"), to deal 
# in the Software without restriction, including without limitation the rights to 
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
# of the Software, and to permit persons to whom the Software is furnished to do 
# so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all 
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
#
#
#
# A square-root raised cosine (SRRC) filter design script for Python. The
# math derivation is from Gaussian Waves:
#
# 	https://www.gaussianwaves.com/2018/10/square-root-raised-cosine-pulse-shaping/
#
#	The design equqations have been changed slightly such that the filter has
#	a gain of 0 dB at f=0.
#
# The SRRC design function requires NumPy to be installed:
#
#	$ sudo apt-get install python3-numpy
#
# To call the function in your Python code, first copy it to the location
# of your script. Then import the function and call it:
#
#	from srrcDesign import srrcDesign
#	span = 4
#	SPS = 2
#	beta = 0.2
#	srrcWeights = srrcDesign(SPS,span,beta)
#
# Parameters:
#
#	SPS: samples per symbol. this can also be considered the oversampling factor.
#		SPS must be an integer greater than 1.
#
#	span: the number of symbols in the filter to the left and right of the center tap.
#		the SRRC filter will have 2*span + 1 symbols in the impulse response.
#		span must be an integer greater than or equal to 1.
#
#	beta: the roll-off factor or excess bandwidth factor. a larger beta results in
#		more occupied bandwidth. beta must be greater than or equal to 0.
#
#

import numpy as np

def srrcDesign ( SPS, span, beta ):

	invalidParameter = False

	if (SPS <= 1):
		print('srrcDesign.py: SPS must be greater than 1')
		invalidParameter = True

	if ((int(SPS) - SPS) > 1e-10):
		print('srrcDesign.py: SPS must be an integer')
		invalidParameter = True

	if (span < 1):
		print('srrcDesign.py: span must be greater than or equal to 1')
		invalidParameter = True

	if ((int(span) - span) > 1e-10):
		print('srrcDesign.py: span must be an integer')
		invalidParameter = True

	if (beta < 0):
		print('srrcDesign.py: beta must be greater than or equal to 0')
		invalidParameter = True

	# check to see if a bad parameter was passed in
	if (invalidParameter):
		# return empty filter weights
		return []

	else: # parameters valid, proceed with design
		# build the time indexing
		nList = np.arange(-span*SPS,(span*SPS)+1)

		# pre-allocate memory for weights
		weights = np.zeros(len(nList))

		# compute the weights on a sample by sample basis
		for index in range(len(nList)):

			# select the time index
			n = nList[index]

			# design equations
			if (n == 0):
				weights[index] = (1/np.sqrt(SPS))*((1-beta) + (4*beta/np.pi))
			elif (np.abs(n*4*beta) == SPS):
				weights[index] = (beta/np.sqrt(2*SPS))*( (1+(2/np.pi))*np.sin(np.pi/(4*beta)) + (1-(2/np.pi))*np.cos(np.pi/(4*beta)) )
			else:
				weights[index] = (1/np.sqrt(SPS))*( (np.sin(np.pi*n*(1-beta)/SPS)) + (4*beta*n/SPS)*(np.cos(np.pi*n*(1+beta)/SPS)) ) / ( (np.pi*n/SPS) * (1 - (4*beta*n/SPS)**2) )

		# scale the weights to 0 dB gain at f=0
		weights = weights/np.sqrt(SPS)

		return weights




