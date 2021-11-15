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
# A half-band filter is efficient due the number of zero weights which do 
# not require computation when filtering. The function abstracts the Remez
# filter design function call for half band filter weights specifically.
# 
# The function requires both NumPy and SciPy. They can be installed by:
#
#	$ sudo apt-get install python3-numpy python3-scipy
#
#
# 
# The function can be called by:
#
# 	from halfBandDesign import halfBandDesign
# 	weights = halfBandDesign(filterLength,transitionBand)
#
# The filterLength parameter determines the length of the filter. filterLength
# must be an integer where filterLength+1 is divisible by 4 and greater than 6.
# For example, filterLength=7, 11, 15, 19, etc.
#
# The transitionBand parameter is the transition bandwidth of the filter
# and must be greater than 0 and less than 0.5.

import numpy as np
import scipy.signal

def halfBandDesign ( filterLength, transitionBand ):

	invalidInput = False

	# check if integer
	if (np.abs(filterLength - int(filterLength)) > 1e-10):
		print('halfBandDesign.py: filterLength must be an integer')
		invalidInput = True

	# check if too small
	if (filterLength < 7):
		print('halfBandDesign.py: filterLength must be larger than 6')
		invalidInput = True

	# check if proper length
	if (np.mod(filterLength+1,4) != 0):
		print('halfBandDesign.py: filterLength+1 must be divisble by 4')
		invalidInput = True

	# check range for transition band
	if (transitionBand <= 0 or transitionBand >= 0.5):
		print('halfBandDesign.py: transitionBand must be greater than 0 and less than 0.5')
		invalidInput = True

	if (invalidInput):
		return []

	else:

		# design a half band filter with remez
		cutoff = 0.25
		fPass = cutoff - (transitionBand/2)
		fStop = cutoff + (transitionBand/2)
		fVec = [0, fPass, fStop, 0.5]
		aVec = [1, 0]

		weights = scipy.signal.remez(filterLength,fVec,aVec)

		# force zero weights
		zeroWeightIndicesHalf = np.arange(2,(filterLength-1)/2,2,dtype=int)
		zeroWeightIndicesNegative = np.concatenate((-zeroWeightIndicesHalf[::-1],zeroWeightIndicesHalf))
		zeroWeightIndices = zeroWeightIndicesNegative - zeroWeightIndicesNegative[0] + 1

		weights[zeroWeightIndices] = 0

		return weights








