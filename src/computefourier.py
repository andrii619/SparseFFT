import math
import cmath
import numpy as np
from scipy import fftpack
from scipy.fftpack import fft, fftshift, ifft
from utils import *
from filters import *

def Comb_Filt(origx, n, num, W_Comb):
	# randomly samples W_Comb elements from input array. Then performs fft on it and finds the largest frequencies. 
	# returns Comb_Approved which contains the num largest frequencies.
	if(n%W_Comb):
		
		
		raise Exception('W_Comb is not divisible by N, which algorithm expects')
	
	
	x_sampt = np.zeros(W_Comb, dtype = complex)
	
	sigma = n/W_Comb;
	offset = int(np.random.rand(1)[0] * sigma)
	
	
	for i in range(W_Comb):
		x_sampt[i] = origx[offset + i*sigma]
		
	
	##idx = np.arange(W_Comb)
	##x_sampt = origx[offset + idx*sigma]
	
	
	x_sampt = fftpack.fft(x_sampt, W_Comb)
	samples = np.zeros(W_Comb, dtype = float)
	
	for i in range(W_Comb):
		
		samples[i] = cabs2(x_sampt[i])
		
		
	
	
	Comb_Approved = find_largest_indices(num,samples)
	
	return Comb_Approved


def inner_loop_locate(origx, n, current_filter, num, B, a, ai, b):
	
	if (n % B):
		print("Warning: n is not divisible by B")
	
	x_sampt = np.zeros(B, dtype = complex)
	index=b
	
	for i in range(current_filter.sig_t.shape[0]):
		
		x_sampt[i%B] += origx[index] * current_filter.sig_t[i]
		index = (index+ai) %n
		
		
	
	x_samp_i = fftpack.fft(x_sampt, B)
	samples = np.power( np.absolute(x_samp_i), 2)
	
	###print("num %d, samples len %d" %(num, len(samples)))
	
	J = find_largest_indices(num, samples)
	return x_samp_i, J




def inner_loop_filter_regular(J, n, num, B, a, ai, b, loop_threshold, score, hits, hits_found):
	
	
	for i in range(num):
		
		low = (int(math.ceil((J[i] - 0.5) * n / B)) + n)%n
		high = (int(math.ceil((J[i] + 0.5) * n / B)) + n)%n
		loc = (low*a)%n
		
		j=low
		while(j!= high):
			
			score[loc] += 1
			if(score[loc] == loop_threshold):
				hits[hits_found] = loc
				hits_found +=1
			
			loc = (loc + a)%n
			
			j = (j+1) %n
		
		
	
	return hits_found




def inner_loop_filter_Comb(J, n, num, B, a, ai, b, loop_threshold, score, hits, hits_found, Comb_Approved, num_Comb, W_Comb):
	
	permuted_approved = np.zeros( (num_Comb, 2), dtype = int )
	
	for m in range(num_Comb):
		
		prev = (Comb_Approved[m] * ai) % W_Comb
		
		permuted_approved[m][0] = prev
		permuted_approved[m][1] = (prev * a) % n
		
	
	permuted_approved = permuted_approved[ np.argsort( permuted_approved[:,0] ) ]
	
	
	
	for i in range(num):
		
		
		low = (int(math.ceil((J[i] - 0.5) * n / B)) + n)%n
		high = (int(math.ceil((J[i] + 0.5) * n / B)) + n)%n
		
		##comp_value = np.array([low%W_Comb, -1])
		
		index = np.where(permuted_approved[:,0] > (low%W_Comb) )[0]
		if(len(index) == 0):
			index = len(permuted_approved)
		else:
			index = index[0]
		
		
		j = index
		
		while(True):
			
			
			if(j == num_Comb):
				
				j -= num_Comb
				location = (location + W_Comb)%n
				locinv = (location * a) % n
				
				
			
			approved_loc = location + permuted_approved[j][0]
			if((low < high and (approved_loc >= high or approved_loc < low)) or (low > high and (approved_loc >= high and approved_loc < low))):
				break
			
			loc = (locinv + permuted_approved[j][1])%n
			score[loc] += 1
			
			if(score[loc]==loop_threshold):
				hits[hits_found]=loc
				hits_found +=1
			
			j += 1
			
		
		
		
		location = low - (low % W_Comb)
		locinv = (location * a) % n
		
		
		
		
	
	
	
	
	
	





def estimate_values(hits, hits_found, x_samp, loops, n, permute, B, B2, filter_loc, filter_est, location_loops):
	
	ans = {}
	values = np.zeros((2, loops), dtype=float)
	
	for i in range(hits_found):
		
		position = 0
		
		
		for j in range(loops):
			
			cur_B = 0
			current_filter = None
			
			if(j<location_loops):
				cur_B = B
				current_filter = filter_loc
			else:
				cur_B = B2
				current_filter = filter_est
			
			
			permuted_index= ((permute[j] * hits[i]) % n)
			hashed_to = permuted_index / (n / cur_B)
			dist = permuted_index % (n / cur_B)
			
			if (dist > ((n/cur_B)/2) ):
				hashed_to = (hashed_to + 1)%cur_B
				dist -= n/cur_B
			
			
			dist = (n - dist) % n
			
			filter_value = current_filter.sig_f[dist]
			values[0][position] = (x_samp[j][hashed_to] / filter_value).real
			values[1][position] = (x_samp[j][hashed_to] / filter_value).imag
			position += 1
			
			
		
		
		location = (loops - 1) / 2
		
		for a in range(2):
			
			values[a] = nth_element(values[a], location)
			
		
		
		realv = values[0][location]
		imagv = values[1][location]
		
		##ans.append( ( hits[i], complex(realv, imagv) ) )
		ans[hits[i]] = complex(realv, imagv)
		
	
	
	
	return ans





def outer_loop(origx,n,filter_loc,filter_est,B2,num,B,W_Comb,Comb_loops,loop_threshold,location_loops,loops, ALG_TYPE):
	# outer_loop(x, n, filter_loc, filter_est, B_est, B_thresh, B_loc, W_Comb, comb_loops, threshold_loops, loc_loops, loc_loops + est_loops)
	# outer_loop(x, n, filter, filter_est, B_est, B_thresh, B_loc, W_Comb, Comb_loops, loops_thresh, loops_loc, loops_loc + loops_est)
	permute  = np.zeros(loops, dtype=int)
	permuteb = np.zeros(loops, dtype=int)
	
	
	
	
	# x_sampt 2d array of loops x B or B2
	x_samp = []
	
	#for i in range(loops):
	#	
	#	tmp = None
	#	
	#	if(i < location_loops):
	#		tmp = np.zeros(B, dtype=complex)
	#	else:
	#		tmp  = np.zeros(B2, dtype=complex)
	#	
	#	
	#	x_samp.append(tmp)
	#	
	
	#x_samp = np.array(x_samp)
	
	
	hits_found = 0
	
	hits = np.zeros(n, dtype=int)
	score = np.zeros(n, dtype = int)
	num_Comb = num
	##Comb_Approved = np.zeros(Comb_loops*num, dtype = int)
	Comb_Approved = []
	
	if(ALG_TYPE == 2):
		# WITH_COMB
		for i in range(Comb_loops):
			
			num_largest = Comb_Filt(origx, n, num, W_Comb)
			Comb_Approved.append(num_largest)
			# add the num_largest indices to the Comb_Approved array of indices
			
		
		
		Comb_Approved = np.array(Comb_Approved)
		Comb_Approved = np.reshape(Comb_Approved, (-1,) )
		Comb_Approved = np.unique(Comb_Approved)
		num_Comb = Comb_Approved.shape[0]
		
		
		hits_found= Comb_Approved.shape[0] * (n/W_Comb)
		for j in range(n/W_Comb):
			for i in range(num_Comb):
				
				hits[j*num_Comb + i] = j*W_Comb + Comb_Approved[i]
				
				
		
		
	
	
	
	
	### inner loops
	
	for i in range(loops):
		
		a = 0
		b = 0
		
		while(gcd(a,n) != 1):
			
			a = np.random.randint(0, n, 1)[0]
			
		
		
		ai = mod_inverse(a, n)
		
		
		permute[i]=ai
		permuteb[i] = b
		
		perform_location = int(i < location_loops)
		
		current_filter = None
		current_B = 0
		
		if(perform_location):
			
			current_filter = filter_loc
			current_B = B
			
		else:
			
			current_filter = filter_est
			current_B  = B2
			
		
		###J = np.zeros(num, dtype = int)
		
		x_samp_i, J = inner_loop_locate(origx, n, current_filter, num, current_B, a, ai, b)
		x_samp.append(x_samp_i)
		
		
		
		if (perform_location):
			if (ALG_TYPE == 1):
				
				#inner_loop_filter_regular(J, n, num, cur_B, a, ai, b, loop_threshold, score, hits, hits_found, G_T)
				hits_found = inner_loop_filter_regular(J, n, num, current_B, a, ai, b, loop_threshold, score, hits, hits_found)
			else:
				hits_found = inner_loop_filter_Comb(J, n, num, current_B, a, ai, b, loop_threshold, score, hits, hits_found, G_T, Comb_Approved, num_Comb, W_Comb)
				
			
			
		
		
		
		
		
		
	
	x_samp = np.array(x_samp)
	
	
	ans = estimate_values(hits, hits_found, x_samp,  loops, n, permute, B, B2, filter_loc, filter_est, location_loops)
	
	return ans













































