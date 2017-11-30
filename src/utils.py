import fractions
import math
import cmath
import numpy as np
from scipy.stats import binom
from scipy import fftpack
from scipy.fftpack import fft, fftshift, ifft

import matplotlib.pyplot as plt


def gcd(a, b):
	return fractions.gcd(a,b)


def extended_euclid(a,b):
	if(b==0):
		return a,1,0
	else:
		d_p, x_p, y_p = ExtendedEuclid(b, a % b)
		d = d_p
		x = y_p
		y = x_p - (a/b)*y_p
		return d, x, y


def mod_inverse(a, m):
	d, x, y = ExtendedEuclid(a, m)
	if d != 1:
		raise Exception('modular inverse does not exist')
	else:
		return x % m


def phase(x):
	return cmath.phase(x)



def binomial_cdf(x, n, p):
	x = np.array(x)
	if(p<0 or p>1.0):
		raise Exception('Probability range is not [0,1]')
	if(n<0):
		raise Exception('Number of trials cannot be negative')
	if( (x<0).any() or (x>n).any() ):
		raise Exception('Cannot evaluate x outside [0,n]')
	return binom.cdf(x, n, p)




def shift(x, r):
	# x - complex array
	# r = shift amount
	# x[i] <- x[i-r]
	x = np.array(x)
	r = r % len(x)
	if(r>=len(x)):
		raise Exception('Shift amount if greater than lenght of x')
	
	x = np.roll(x, r, axis = 0)
	return x



def left_shift(x, r):
	
	x = np.array(x)
	r = r % len(x)
	if(r>=len(x)):
		raise Exception('Shift amount if greater than lenght of x')
	
	x = np.roll(x, len(x) - r, axis = 0)
	
	return x



def floor_to_pow2(x):
	ans = 1
	
	while(ans<=x):
		ans <<= 1
	
	return (ans / 2)



def cabs2(x):
	# square of abs value
	
	return abs(x)**2




def nth_element_immutable(x, n):
	# get the nth smallest element of the array
	x = np.array(x)
	if(n >= len(x) or n<0):
		raise Exception('index outside bounds')
	
	idx = np.argpartition(x, n)
	srtd = x[idx]
	
	return srtd[n]






def get_expermient_vs_N_parameters(N, WITH_COMB):
	
	params = {}
	
	if(N == 8192):
		if(WITH_COMB):
			##Bcst_loc =  2; Bcst_est =  2; Comb_cst = 32; comb_loops =8; est_loops =16; loc_loops =7; threshold_loops =6; tolerance_loc =1e-8; tolerance_est =1e-8;
			params['Bcst_loc'] = 2
			params['Bcst_est'] = 2
			params['Comb_cst'] = 32
			params['comb_loops'] = 8
			params['est_loops'] = 16
			params['loc_loops'] = 7
			params['threshold_loops'] = 6
			params['tolerance_loc'] = 1e-8
			params['tolerance_est'] = 1e-8
		else:
			##Bcst_loc =2; Bcst_est =  2; Comb_cst =1; comb_loops =1; est_loops =16; loc_loops =7; threshold_loops =6; tolerance_loc =1e-8; tolerance_est =1e-8;
			params['Bcst_loc'] = 2
			params['Bcst_est'] = 2
			params['Comb_cst'] = 1
			params['comb_loops'] = 1
			params['est_loops'] = 16
			params['loc_loops'] = 7
			params['threshold_loops'] = 6
			params['tolerance_loc'] = 1e-8
			params['tolerance_loc'] = 1e-8
		
		
	elif(N == 16384):
		if(WITH_COMB):
			##Bcst_loc =  4; Bcst_est =  4; Comb_cst = 32; comb_loops =8; est_loops =10; loc_loops =6; threshold_loops =5; tolerance_loc =1e-8; tolerance_est =1e-8;
			params['Bcst_loc'] = 4
			params['Bcst_est'] = 4
			params['Comb_cst'] = 32
			params['comb_loops'] = 8
			params['est_loops'] = 10
			params['loc_loops'] = 6
			params['threshold_loops'] = 5
			params['tolerance_loc'] = 1e-8
			params['tolerance_est'] = 1e-8
		else:
			params['Bcst_loc'] = 4
			params['Bcst_est'] = 4
			params['Comb_cst'] = 1
			params['comb_loops'] = 1
			params['est_loops'] = 10
			params['loc_loops'] = 6
			params['threshold_loops'] = 5
			params['tolerance_loc'] = 1e-8
			params['tolerance_est'] = 1e-8
			##Bcst_loc =4; Bcst_est =  4; Comb_cst =1; comb_loops =1; est_loops =10; loc_loops =6; threshold_loops =5; tolerance_loc =1e-8; tolerance_est =1e-8;
		
	elif(N == 32768):
		if(WITH_COMB):
			##Bcst_loc =  4; Bcst_est =  2; Comb_cst = 64; comb_loops =4; est_loops = 8; loc_loops =5; threshold_loops =4; tolerance_loc =1e-8; tolerance_est =1e-8;
			params['Bcst_loc'] = 4
			params['Bcst_est'] = 2
			params['Comb_cst'] = 64
			params['comb_loops'] = 4
			params['est_loops'] = 8
			params['loc_loops'] = 5
			params['threshold_loops'] = 4
			params['tolerance_loc'] = 1e-8
			params['tolerance_loc'] = 1e-8
			
		else:
			##Bcst_loc =4; Bcst_est =  2; Comb_cst =1; comb_loops =1; est_loops = 8; loc_loops =5; threshold_loops =4; tolerance_loc =1e-8; tolerance_est =1e-8;
			params['Bcst_loc'] = 4
			params['Bcst_est'] = 2
			params['Comb_cst'] = 1
			params['comb_loops'] = 1
			params['est_loops'] = 8
			params['loc_loops'] = 5
			params['threshold_loops'] = 4
			params['tolerance_loc'] = 1e-8
			params['tolerance_est'] = 1e-8
		
		
	elif(N == 65536):
		if(WITH_COMB):
			##Bcst_loc =  4; Bcst_est =  2; Comb_cst =128; comb_loops =6; est_loops =10; loc_loops =4; threshold_loops =2; tolerance_loc =1e-8; tolerance_est =1e-8;
			params['Bcst_loc'] = 4
			params['Bcst_est'] = 2
			params['Comb_cst'] = 128
			params['comb_loops'] = 6
			params['est_loops'] = 10
			params['loc_loops'] = 4
			params['threshold_loops'] = 2
			params['tolerance_loc'] = 1e-8
			params['tolerance_loc'] = 1e-8
		else:
			##Bcst_loc =4; Bcst_est =  2; Comb_cst =1; comb_loops =1; est_loops = 8; loc_loops =5; threshold_loops =4; tolerance_loc =1e-8; tolerance_est =1e-8;
			params['Bcst_loc'] = 4
			params['Bcst_est'] = 2
			params['Comb_cst'] = 1
			params['comb_loops'] = 1
			params['est_loops'] = 8
			params['loc_loops'] = 5
			params['threshold_loops'] = 4
			params['tolerance_loc'] = 1e-8
			params['tolerance_est'] = 1e-8
		
	elif(N == 131072):
		if(WITH_COMB):
			params['Bcst_loc'] = 1
			params['Bcst_est'] = 1
			params['Comb_cst'] = 8
			params['comb_loops'] = 2
			params['est_loops'] = 12
			params['loc_loops'] = 4
			params['threshold_loops'] = 3
			params['tolerance_loc'] = 1e-6
			params['tolerance_loc'] = 1e-8
			##Bcst_loc =  1; Bcst_est =  1; Comb_cst =  8; comb_loops =2; est_loops =12; loc_loops =4; threshold_loops =3; tolerance_loc =1e-6; tolerance_est =1e-8;
			
			
		else:
			##Bcst_loc =2; Bcst_est =  1; Comb_cst =1; comb_loops =1; est_loops =10; loc_loops =5; threshold_loops =4; tolerance_loc =1e-6; tolerance_est =1e-8;
			params['Bcst_loc'] = 2
			params['Bcst_est'] = 1
			params['Comb_cst'] = 1
			params['comb_loops'] = 1
			params['est_loops'] = 10
			params['loc_loops'] = 5
			params['threshold_loops'] = 4
			params['tolerance_loc'] = 1e-6
			params['tolerance_est'] = 1e-8
			
		
		
	elif(N == 262144):
		if(WITH_COMB):
			##Bcst_loc =  1; Bcst_est =  1; Comb_cst =  8; comb_loops =2; est_loops =14; loc_loops =5; threshold_loops =4; tolerance_loc =1e-6; tolerance_est =1e-8;
			params['Bcst_loc'] = 1
			params['Bcst_est'] = 1
			params['Comb_cst'] = 8
			params['comb_loops'] = 2
			params['est_loops'] = 14
			params['loc_loops'] = 5
			params['threshold_loops'] = 4
			params['tolerance_loc'] = 1e-6
			params['tolerance_loc'] = 1e-8
			
		else:
			##Bcst_loc =2; Bcst_est =0.5; Comb_cst =1; comb_loops =1; est_loops =14; loc_loops =4; threshold_loops =3; tolerance_loc =1e-6; tolerance_est =1e-8;
			params['Bcst_loc'] = 2
			params['Bcst_est'] = 0.5
			params['Comb_cst'] = 1
			params['comb_loops'] = 1
			params['est_loops'] = 14
			params['loc_loops'] = 4
			params['threshold_loops'] = 3
			params['tolerance_loc'] = 1e-6
			params['tolerance_est'] = 1e-8
		
		
	elif(N == 524288):
		
		if(WITH_COMB):
			##Bcst_loc =0.5; Bcst_est =0.5; Comb_cst =  8; comb_loops =1; est_loops =10; loc_loops =4; threshold_loops =3; tolerance_loc =1e-6; tolerance_est =1e-8;
			params['Bcst_loc'] = 0.5
			params['Bcst_est'] = 0.5
			params['Comb_cst'] = 8
			params['comb_loops'] = 1
			params['est_loops'] = 10
			params['loc_loops'] = 4
			params['threshold_loops'] = 3
			params['tolerance_loc'] = 1e-6
			params['tolerance_loc'] = 1e-8
			
			
		else:
			##Bcst_loc =1; Bcst_est =0.5; Comb_cst =1; comb_loops =1; est_loops =12; loc_loops =5; threshold_loops =4; tolerance_loc =1e-6; tolerance_est =1e-8;
			params['Bcst_loc'] = 1
			params['Bcst_est'] = 0.5
			params['Comb_cst'] = 1
			params['comb_loops'] = 1
			params['est_loops'] = 12
			params['loc_loops'] = 5
			params['threshold_loops'] = 4
			params['tolerance_loc'] = 1e-6
			params['tolerance_est'] = 1e-8
		
		
	elif(N == 1048576):
		
		if(WITH_COMB):
			##Bcst_loc =0.5; Bcst_est =0.5; Comb_cst =  8; comb_loops =2; est_loops =12; loc_loops =4; threshold_loops =2; tolerance_loc =1e-6; tolerance_est =1e-8;
			params['Bcst_loc'] = 0.5
			params['Bcst_est'] = 0.5
			params['Comb_cst'] = 8
			params['comb_loops'] = 2
			params['est_loops'] = 12
			params['loc_loops'] = 4
			params['threshold_loops'] = 2
			params['tolerance_loc'] = 1e-6
			params['tolerance_loc'] = 1e-8
		else:
			##Bcst_loc =2; Bcst_est =0.5; Comb_cst =1; comb_loops =1; est_loops =12; loc_loops =4; threshold_loops =3; tolerance_loc =1e-6; tolerance_est =1e-8;
			params['Bcst_loc'] = 2
			params['Bcst_est'] = 0.5
			params['Comb_cst'] = 1
			params['comb_loops'] = 1
			params['est_loops'] = 12
			params['loc_loops'] = 4
			params['threshold_loops'] = 3
			params['tolerance_loc'] = 1e-6
			params['tolerance_est'] = 1e-8
		
		
	elif(N == 2097152):
		
		if(WITH_COMB):
			##Bcst_loc =0.5; Bcst_est =0.2; Comb_cst =  8; comb_loops =1; est_loops =10; loc_loops =3; threshold_loops =2; tolerance_loc =1e-6; tolerance_est =1e-8;
			params['Bcst_loc'] = 0.5
			params['Bcst_est'] = 0.2
			params['Comb_cst'] = 8
			params['comb_loops'] = 1
			params['est_loops'] = 10
			params['loc_loops'] = 3
			params['threshold_loops'] = 2
			params['tolerance_loc'] = 1e-6
			params['tolerance_loc'] = 1e-8
			
		else:
			##Bcst_loc =2; Bcst_est =0.2; Comb_cst =1; comb_loops =1; est_loops =15; loc_loops =3; threshold_loops =2; tolerance_loc =1e-6; tolerance_est =1e-8;
			params['Bcst_loc'] = 2
			params['Bcst_est'] = 0.2
			params['Comb_cst'] = 1
			params['comb_loops'] = 1
			params['est_loops'] = 15
			params['loc_loops'] = 3
			params['threshold_loops'] = 2
			params['tolerance_loc'] = 1e-6
			params['tolerance_est'] = 1e-8
		
	elif(N == 4194304):
		
		if(WITH_COMB):
			##Bcst_loc =0.5; Bcst_est =0.2; Comb_cst =  8; comb_loops =1; est_loops = 8; loc_loops =3; threshold_loops =2; tolerance_loc =1e-6; tolerance_est =1e-8;
			params['Bcst_loc'] = 0.5
			params['Bcst_est'] = 0.2
			params['Comb_cst'] = 8
			params['comb_loops'] = 1
			params['est_loops'] = 8
			params['loc_loops'] = 3
			params['threshold_loops'] = 2
			params['tolerance_loc'] = 1e-6
			params['tolerance_loc'] = 1e-8
			
			
		else:
			##Bcst_loc =4; Bcst_est =0.2; Comb_cst =1; comb_loops =1; est_loops =10; loc_loops =3; threshold_loops =2; tolerance_loc =1e-6; tolerance_est =1e-8;
			params['Bcst_loc'] = 4
			params['Bcst_est'] = 0.2
			params['Comb_cst'] = 1
			params['comb_loops'] = 1
			params['est_loops'] = 10
			params['loc_loops'] = 3
			params['threshold_loops'] = 2
			params['tolerance_loc'] = 1e-6
			params['tolerance_est'] = 1e-8
		
	elif(N == 8388608):
		
		if(WITH_COMB):
			##Bcst_loc =0.5; Bcst_est =0.2; Comb_cst =  8; comb_loops =1; est_loops = 8; loc_loops =3; threshold_loops =2; tolerance_loc =1e-6; tolerance_est =1e-8;
			params['Bcst_loc'] = 0.5
			params['Bcst_est'] = 0.2
			params['Comb_cst'] = 8
			params['comb_loops'] = 1
			params['est_loops'] = 8
			params['loc_loops'] = 3
			params['threshold_loops'] = 2
			params['tolerance_loc'] = 1e-6
			params['tolerance_loc'] = 1e-8
			
			
		else:
			##Bcst_loc =2; Bcst_est =0.2; Comb_cst =1; comb_loops =1; est_loops = 8; loc_loops =3; threshold_loops =2; tolerance_loc =1e-6; tolerance_est =1e-8;
			params['Bcst_loc'] = 2
			params['Bcst_est'] = 0.2
			params['Comb_cst'] = 1
			params['comb_loops'] = 1
			params['est_loops'] = 8
			params['loc_loops'] = 3
			params['threshold_loops'] = 2
			params['tolerance_loc'] = 1e-6
			params['tolerance_est'] = 1e-8
			
		
		
	elif(N == 16777216):
		
		if(WITH_COMB):
			##Bcst_loc =0.5; Bcst_est =0.2; Comb_cst = 16; comb_loops =1; est_loops = 8; loc_loops =3; threshold_loops =2; tolerance_loc =1e-6; tolerance_est =1e-8;
			params['Bcst_loc'] = 0.5
			params['Bcst_est'] = 0.2
			params['Comb_cst'] = 16
			params['comb_loops'] = 1
			params['est_loops'] = 8
			params['loc_loops'] = 3
			params['threshold_loops'] = 2
			params['tolerance_loc'] = 1e-6
			params['tolerance_loc'] = 1e-8
			
		else:
			##Bcst_loc =4; Bcst_est =0.2; Comb_cst =1; comb_loops =1; est_loops = 8; loc_loops =3; threshold_loops =2; tolerance_loc =1e-6; tolerance_est =1e-8;
			params['Bcst_loc'] = 4
			params['Bcst_est'] = 0.2
			params['Comb_cst'] = 1
			params['comb_loops'] = 1
			params['est_loops'] = 8
			params['loc_loops'] = 3
			params['threshold_loops'] = 2
			params['tolerance_loc'] = 1e-6
			params['tolerance_est'] = 1e-8
			
		
		
	else:
		
		
		params['Bcst_loc'] = 1
		params['Bcst_est'] = 1
		params['Comb_cst'] = 2
		params['comb_loops'] = 1
		params['est_loops'] = 16
		params['loc_loops'] = 4
		params['threshold_loops'] = 3
		params['tolerance_loc'] = 1e-8
		params['tolerance_est'] = 1e-8
		
		
	
	return params



def get_expermient_vs_K_parameters(K, WITH_COMB):
	
	params= {}
	
	if(K == 50):
		
		if(WITH_COMB):
			#Bcst_loc =0.5; Bcst_est =0.2; Comb_cst = 16; comb_loops =1; est_loops =10; loc_loops =3; threshold_loops =2; tolerance_loc =1e-6; tolerance_est =1.0e-8;
			params['Bcst_loc'] = 0.5
			params['Bcst_est'] = 0.2
			params['Comb_cst'] = 16
			params['comb_loops'] = 1
			params['est_loops'] = 10
			params['loc_loops'] = 3
			params['threshold_loops'] = 2
			params['tolerance_loc'] = 1e-6
			params['tolerance_est'] = 1e-8
			
		else:
			
			##Bcst_loc =4; Bcst_est =0.2; Comb_cst =1; comb_loops =1; est_loops =10; loc_loops =3; threshold_loops =2; tolerance_loc =1e-6; tolerance_est =1.0e-8;
			params['Bcst_loc'] = 4
			params['Bcst_est'] = 0.2
			params['Comb_cst'] = 1
			params['comb_loops'] = 1
			params['est_loops'] = 10
			params['loc_loops'] = 3
			params['threshold_loops'] = 2
			params['tolerance_loc'] = 1e-6
			params['tolerance_est'] = 1e-8
			
		
		
		
	elif(K == 100):
		
		
		if(WITH_COMB):
			##Bcst_loc =0.5; Bcst_est =0.2; Comb_cst = 16; comb_loops =1; est_loops =12; loc_loops =4; threshold_loops =2; tolerance_loc =1e-6; tolerance_est =1.0e-8;
			params['Bcst_loc'] = 0.5
			params['Bcst_est'] = 0.2
			params['Comb_cst'] = 16
			params['comb_loops'] = 1
			params['est_loops'] = 12
			params['loc_loops'] = 4
			params['threshold_loops'] = 2
			params['tolerance_loc'] = 1e-6
			params['tolerance_est'] = 1e-8
			
		else:
			##Bcst_loc =2; Bcst_est =0.2; Comb_cst =1; comb_loops =1; est_loops =12; loc_loops =3; threshold_loops =2; tolerance_loc =1e-6; tolerance_est =1.0e-8;
			params['Bcst_loc'] = 2
			params['Bcst_est'] = 0.2
			params['Comb_cst'] = 1
			params['comb_loops'] = 1
			params['est_loops'] = 12
			params['loc_loops'] = 3
			params['threshold_loops'] = 2
			params['tolerance_loc'] = 1e-6
			params['tolerance_est'] = 1e-8
			
		
		
		
	elif(K == 200):
		
		
		if(WITH_COMB):
			##Bcst_loc =  0.5; Bcst_est =  0.5; Comb_cst = 32; comb_loops =1; est_loops = 8; loc_loops =4; threshold_loops =3; tolerance_loc =1e-6; tolerance_est =0.5e-8;
			params['Bcst_loc'] = 0.5
			params['Bcst_est'] = 0.5
			params['Comb_cst'] = 32
			params['comb_loops'] = 1
			params['est_loops'] = 8
			params['loc_loops'] = 4
			params['threshold_loops'] = 3
			params['tolerance_loc'] = 1e-6
			params['tolerance_est'] = 0.5e-8
		else:
			##Bcst_loc =4; Bcst_est =0.5; Comb_cst =1; comb_loops =1; est_loops =10; loc_loops =3; threshold_loops =2; tolerance_loc =1e-6; tolerance_est =0.5e-8;
			params['Bcst_loc'] = 4
			params['Bcst_est'] = 0.5
			params['Comb_cst'] = 1
			params['comb_loops'] = 1
			params['est_loops'] = 10
			params['loc_loops'] = 3
			params['threshold_loops'] = 2
			params['tolerance_loc'] = 1e-6
			params['tolerance_est'] = 0.5e-8
			
		
		
		
	elif(K == 500):
		
		
		if(WITH_COMB):
			##Bcst_loc =0.5; Bcst_est =0.5; Comb_cst = 64; comb_loops =1; est_loops =10; loc_loops =4; threshold_loops =3; tolerance_loc =1e-6; tolerance_est =0.5e-8;
			params['Bcst_loc'] = 0.5
			params['Bcst_est'] = 0.5
			params['Comb_cst'] = 64
			params['comb_loops'] = 1
			params['est_loops'] = 10
			params['loc_loops'] = 4
			params['threshold_loops'] = 3
			params['tolerance_loc'] = 1e-6
			params['tolerance_est'] = 0.5e-8
		else:
			###Bcst_loc =2; Bcst_est =  1; Comb_cst =1; comb_loops =1; est_loops =12; loc_loops =4; threshold_loops =3; tolerance_loc =1e-6; tolerance_est =0.5e-8;
			params['Bcst_loc'] = 2
			params['Bcst_est'] = 1
			params['Comb_cst'] = 1
			params['comb_loops'] = 1
			params['est_loops'] = 12
			params['loc_loops'] = 4
			params['threshold_loops'] = 3
			params['tolerance_loc'] = 1e-6
			params['tolerance_est'] = 0.5e-8
			
		
		
		
	elif(K == 1000):
		
		if(WITH_COMB):
			##Bcst_loc =  1; Bcst_est =  1; Comb_cst =128; comb_loops =3; est_loops =12; loc_loops =4; threshold_loops =3; tolerance_loc =1e-6; tolerance_est =0.5e-8;
			params['Bcst_loc'] = 1
			params['Bcst_est'] = 1
			params['Comb_cst'] = 128
			params['comb_loops'] = 3
			params['est_loops'] = 12
			params['loc_loops'] = 4
			params['threshold_loops'] = 3
			params['tolerance_loc'] = 1e-6
			params['tolerance_est'] = 0.5e-8
		else:
			##Bcst_loc =2; Bcst_est =  1; Comb_cst =1; comb_loops =1; est_loops =12; loc_loops =5; threshold_loops =4; tolerance_loc =1e-6; tolerance_est =1.0e-8;
			params['Bcst_loc'] = 2
			params['Bcst_est'] = 1
			params['Comb_cst'] = 1
			params['comb_loops'] = 1
			params['est_loops'] = 12
			params['loc_loops'] = 5
			params['threshold_loops'] = 4
			params['tolerance_loc'] = 1e-6
			params['tolerance_est'] = 1e-8
			
		
		
		
		
	elif(K == 2000):
		
		
		if(WITH_COMB):
			
			##Bcst_loc =  1; Bcst_est =  1; Comb_cst =512; comb_loops =3; est_loops =16; loc_loops =4; threshold_loops =3; tolerance_loc =1e-7; tolerance_est =0.2e-8;
			params['Bcst_loc'] = 1
			params['Bcst_est'] = 1
			params['Comb_cst'] = 512
			params['comb_loops'] = 3
			params['est_loops'] = 16
			params['loc_loops'] = 4
			params['threshold_loops'] = 3
			params['tolerance_loc'] = 1e-7
			params['tolerance_est'] = 0.2e-8
		else:
			
			#Bcst_loc =2; Bcst_est =  1; Comb_cst =1; comb_loops =1; est_loops =16; loc_loops =5; threshold_loops =4; tolerance_loc =1e-7; tolerance_est =0.5e-8;
			params['Bcst_loc'] = 2
			params['Bcst_est'] = 1
			params['Comb_cst'] = 1
			params['comb_loops'] = 1
			params['est_loops'] = 16
			params['loc_loops'] = 5
			params['threshold_loops'] = 4
			params['tolerance_loc'] = 1e-7
			params['tolerance_est'] = 0.5e-8
			
		
		
	elif(K == 2500):
		
		if(WITH_COMB):
			
			##Bcst_loc =  1; Bcst_est =  1; Comb_cst =512; comb_loops =3; est_loops =16; loc_loops =4; threshold_loops =3; tolerance_loc =1e-7; tolerance_est =0.2e-8;
			params['Bcst_loc'] = 1
			params['Bcst_est'] = 1
			params['Comb_cst'] = 512
			params['comb_loops'] = 3
			params['est_loops'] = 16
			params['loc_loops'] = 4
			params['threshold_loops'] = 3
			params['tolerance_loc'] = 1e-7
			params['tolerance_est'] = 0.2e-8
		else:
			
			##Bcst_loc =2; Bcst_est =  1; Comb_cst =1; comb_loops =1; est_loops =16; loc_loops =5; threshold_loops =4; tolerance_loc =1e-7; tolerance_est =0.5e-8;
			params['Bcst_loc'] = 2
			params['Bcst_est'] = 1
			params['Comb_cst'] = 1
			params['comb_loops'] = 1
			params['est_loops'] = 16
			params['loc_loops'] = 5
			params['threshold_loops'] = 4
			params['tolerance_loc'] = 1e-7
			params['tolerance_est'] = 0.5e-8
			
		
		
		
	elif(K == 4000):
		
		if(WITH_COMB):
			##Bcst_loc =  1; Bcst_est =  2; Comb_cst =512; comb_loops =3; est_loops =14; loc_loops =8; threshold_loops =7; tolerance_loc =1e-8; tolerance_est =0.5e-8;
			params['Bcst_loc'] = 1
			params['Bcst_est'] = 2
			params['Comb_cst'] = 512
			params['comb_loops'] = 3
			params['est_loops'] = 14
			params['loc_loops'] = 8
			params['threshold_loops'] = 7
			params['tolerance_loc'] = 1e-8
			params['tolerance_est'] = 0.5e-8
		else:
			##Bcst_loc =2; Bcst_est =  2; Comb_cst =1; comb_loops =1; est_loops =14; loc_loops =6; threshold_loops =5; tolerance_loc =1e-8; tolerance_est =1.0e-8;
			params['Bcst_loc'] = 2
			params['Bcst_est'] = 2
			params['Comb_cst'] = 1
			params['comb_loops'] = 1
			params['est_loops'] = 14
			params['loc_loops'] = 6
			params['threshold_loops'] = 5
			params['tolerance_loc'] = 1e-8
			params['tolerance_est'] = 1e-8
			
		
		
		
	else:
		
		
		params['Bcst_loc'] = 1
		params['Bcst_est'] = 1
		params['Comb_cst'] = 2
		params['comb_loops'] = 1
		params['est_loops'] = 16
		params['loc_loops'] = 4
		params['threshold_loops'] = 3
		params['tolerance_loc'] = 1e-8
		params['tolerance_est'] = 1e-8
		
	
	
	return params














def generate_random_signal(n, k):
	
	x_f = np.zeros(n)
	large_freq = np.random.randint(0, n, k)
	# generate k random indecies for the sparce k frequencies
	x_f[large_freq] = 1.0
	
	#plt.stem(x_f)
	#plt.show()
	
	
	
	x = fftpack.fft(x_f, n)
	
	return x, x_f






































