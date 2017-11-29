import numpy as np
from scipy import signal
import math
import cmath
from scipy import fftpack
from scipy.fftpack import fft, fftshift
import matplotlib.pyplot as plt
from scipy.signal import convolve

from utils import *

def box_window(length):
	return signal.boxcar(length)







def hamming_window(length):
	#x = np.ones(length)
	#i = np.linspace(0, length-1, length)
	#x = 0.54 - 0.46 * np.cos(2*math.pi * i / (length - 1))
	return signal.hamming(length)











def kaiserbessel_window(lobefrac, tolerance):
	w = int((1.0 / math.pi) * (1.0/lobefrac) * math.acosh(1.0/tolerance))
	beta = np.log(1.0/tolerance)
	return signal.kaiser(w, beta)


def gaussian_window(lobefrac, tolerance):
	w = int((2.0 / math.pi) * (1.0/lobefrac) * np.log(1.0/tolerance))
	std = (w/2.0) / math.sqrt(2.0*np.log(1.0/tolerance))
	x = signal.gaussian(w, std)
	return x


def chebyshev_window(lobefrac, tolerance):
	w = int((1.0 / math.pi) * (1.0/lobefrac) * np.log(1.0/tolerance))
	if ((w%2) != 0):
		w -= 1
	
	x0 = math.cosh(math.acosh(1.0/tolerance) / (w-1))
	w0 = 2.0*math.acos(1.0/x0)
	at = -20*math.log10(1.0/tolerance)
	x = signal.chebwin(w, at)
	return x





def make_flat_window(filtert, boxcar_length):
	
	# convolve filtert and a boxcar filter
	if(len(filtert) < 0):
		return None
	
	box = box_window(boxcar_length)
	
	
	
	return convolve(filtert, box)









class Filter():
	
	
	def __init__(self, x_t, x_f):
		
		
		self.sig_t = x_t
		self.sig_f = x_f
		
		
		
	


def make_multiple(x, n, b):
	
	#n=4096
	#tol = 1e-9
	#frac = 0.05
	#b=32
	#x = gaussian_window(frac, tol)
	
	w = x.shape[0]
	
	h = np.zeros(n)
	g = np.pad(x, (0, n- w), mode='constant', constant_values=0)
	g = left_shift(g, w/2)
	
	g = fftpack.fft(g, n)
	
	s = np.sum( g[:b])
	
	max_g = 0
	offset = b/2
	
	for i in range(n):
		h[(i + offset)%n] = s
		max_g = np.maximum(max_g, np.abs(s) )
		s = s + (g[(i + b)%n] - g[i])
	
	
	h = np.divide(h, max_g)
	
	
	
	offset_c = 1
	step = np.exp(-2.0*math.pi * complex(0,1) * (w/2) / n)
	
	for i in range(n):
		h[i] = np.multiply(h[i], offset_c)
		offset_c = np.multiply(offset_c, step)
	
	
	window = fftpack.fft(h, n)
	
	###window[w:n] = 0.0
	window = np.divide(window, n/2)
	
	x = window[:w]
	##x = np.pad(x, (0, n- w), mode='constant', constant_values=0)
	f = Filter(x, h)
	
	
	plt.plot(x)
	plt.show()
	
	
	plt.plot(h)
	plt.show()
	
	
	
	return f











