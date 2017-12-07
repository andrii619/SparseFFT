import numpy as np
from scipy import signal
import math
import cmath
from scipy import fftpack
from scipy.fftpack import fft, fftshift, ifft
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



#def I0(x):
#	ans = 1
#	curval = 1
#	for(int L = 1; curval > 0.001; L++)
#	{
#		curval = curval * (x / 2)*(x / 2) / L / L;
#		ans += curval;
#		//printf("State: %d %lf %lf\n", L, curval, ans);
#	}
#	return ans








def kaiserbessel_window(lobefrac, tolerance):
	w = int((1.0 / math.pi) * (1.0/lobefrac) * math.acosh(1.0/tolerance))
	beta = np.log(1.0/tolerance)
	return signal.kaiser(w, beta)



#def kaiserbessel_window(double lobefrac, double tolerance, int &w):
#	w = int((1.0 / math.pi) * (1.0/lobefrac) * math.acosh(1.0/tolerance))
#	B = np.log(1.0/tolerance)
#	##complex_t *x = (complex_t *)malloc(w*sizeof(*x))
#	x = np.zeros(w, dtype=complex)
#	for i in range(w):
#		tmp = (2.0 * (i - (w-1)/2.0)) / (w - 1)
#		x[i] = I0(B * math.sqrt(1 - tmp*tmp)) / I0(B)
#	
#	return x



def gaussian_window(lobefrac, tolerance):
	w = int((2.0 / math.pi) * (1.0/lobefrac) * np.log(1.0/tolerance))
	std = (w/2.0) / math.sqrt(2.0*np.log(1.0/tolerance))
	x = signal.gaussian(w, std)
	return x


def chebyshev_window2(lobefrac, tolerance):
	w = int((1 / math.pi) * (1/lobefrac) * math.acosh(1.0/tolerance))
	if ((w%2) != 0):
		w -= 1
	
	x0 = math.cosh(math.acosh(1.0/tolerance) / (w-1))
	w0 = 2.0*math.acos(1.0/x0)
	at = -20*math.log10(1.0/tolerance)
	x = signal.chebwin(w, at)
	return x


def Cheb(m, x):
	if(abs(x) <= 1):
		return math.cos(m * math.acos(x))
	else:
		return (cmath.cosh(m * cmath.acosh(x))).real




def chebyshev_window(lobefrac, tolerance):
	w = int((1.0 / math.pi) * (1.0/lobefrac) * math.acosh(1.0/tolerance))
	if (not (w%2)):
		w -= 1
	x = np.zeros(w)
	t0 = math.cosh(math.acosh(1.0/tolerance) / (w-1))
	for i in range(w):
		
		x[i] = Cheb(w-1, t0 * math.cos(math.pi * i / w)) * tolerance
		
	
	x = fftpack.fft(x, w)
	x = fftpack.fftshift(x)
	##x = shift(x, w/2)
	for i in range(w):
		x[i] = x[i].real
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
	
	x = x.astype(complex)
	w = x.shape[0]
	h = np.zeros(n, dtype=complex)
	g = np.pad(x, (0, n- w), mode='constant', constant_values= complex(0.0,0.0))
	g = left_shift(g, w/2)
	
	
	
	g = fftpack.fft(g, n)
	
	s = np.sum( g[:b])
	max_g = 0
	offset = int(b/2)
	
	for i in range(n):
		h[(i+n + offset)%n] = s
		max_g = np.maximum(max_g, abs(s) )
		s = s + (g[(i + b)%n] - g[i])
	
	h = np.divide(h, max_g)
	offset_c = complex(1,0)
	step = cmath.exp(-2.0*math.pi * complex(0,1) * (w/2) / n)
	
	for i in range(n):
		h[i] = (h[i] * offset_c)
		offset_c = (offset_c * step)
	
	
	window = fftpack.fft(h, n)/n
	window = window[::-1]
	window = np.roll(window, 1)
	
	
	#plt.plot(abs(window))
	#plt.show()
	
	#for l in range(n):
	#	
	#	print("g[%d] = %3.8f" % (l, window[l].real))
	#	
	
	
	####window = window[::-1]
	###window = left_shift(window, (n-w-1))
	###window = np.divide(window, n)
	x = window[:w]
	f = Filter(x, h)
	
	return f










#n = 8192
#k=50
#Bcst_loc = 4
#Bcst_est = 4
#tolerance_loc = 1e-8
#tolerance_est = 1e-8

#BB_loc = Bcst_loc*math.sqrt( (n*k) /math.log(n, 2) ) 
#BB_est = Bcst_est*math.sqrt( (n*k)/math.log(n, 2) ) 

#lobefrac_loc = 0.5/BB_loc
#lobefrac_est = 0.5/BB_est
#b_loc = int(1.2*1.1*(n/BB_loc))
#b_est = int(1.4*1.1*(n/BB_est))
#b=b_loc


#w = int((1.0 / math.pi) * (1.0/lobefrac_loc) * math.acosh(1.0/tolerance_loc))

#filter_loc = chebyshev_window(lobefrac_loc, tolerance_loc)
#x = filter_loc
#############################################
#x = x.astype(complex)
#w = x.shape[0]
#h = np.zeros(n, dtype=complex)
#g = np.pad(x, (0, n- w), mode='constant', constant_values= complex(0.0,0.0))
#g = left_shift(g, w/2)
#g = fftpack.fft(g, n)
#s = np.sum( g[:b])
#max_g = 0
#offset = int(b/2)

#for i in range(n):
#	h[(i+n + offset)%n] = s
#	max_g = np.maximum(max_g, abs(s) )
#	s = s + (g[(i + b)%n] - g[i])

#h = np.divide(h, max_g)
#offset_c = complex(1,0)
#step = cmath.exp(-2.0*math.pi * complex(0,1) * (w/2) / n)

#for i in range(n):
#	h[i] = (h[i] * offset_c)
#	offset_c = (offset_c * step)

#window = fftpack.fft(h, n)
#window = left_shift(window, (n-w))
#window = np.divide(window, n)
#x = window[:w]



#n = 8192
#k=10
#Bcst_loc =2
#tolerance_loc =1e-8
#BB_loc = math.floor( Bcst_loc*np.sqrt( (n*k) /np.log2(n) ) )
#lobefrac_loc = 0.5/BB_loc
#filter_t = chebyshev_window(lobefrac_loc, tolerance_loc)
#w_loc = filter_t.shape[0]

#b_loc = int(1.2*1.1*(n/BB_loc))


#filter_loc = make_multiple(filter_t, n, b_loc)

#plt.plot(abs(filter_loc.sig_t))
#plt.show()

#for l in range(w_loc):
#	
#	print("filt[%d] = %3.8f" %(l, filter_loc.sig_t[l].real*1000.0))
#	

























