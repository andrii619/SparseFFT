import math
import cmath
import numpy as np
from scipy import fftpack
from scipy.fftpack import fft, fftshift, ifft


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





def inner_loop_locate(origx, n, current_filter, num, B, a, ai, b, J):
	
	
	if (n % B):
		print("Warning: n is not divisible by B")
	
	x_sampt = np.zeros(n, dtype = complex)
	
	index=b
	
	for i in range(current_filter.sig_t.shape[0]):
		
		x_sampt[i%B] += origx[index] * current_filter.sig_t[i]
		index = (index+ai) %n
		
		
	
	
	x_samp = fftpack.fft(x_sampt, B)
	
	samples = np.power( np.absolute(x_samp), 2)
	
	J = find_largest_indices(num, samples)
	
	
	
	
	
	return x_samp, J




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
			
			if(score[loc]==loop_threshold)
				hits[hits_found]=loc
				hits_found +=1
			
			j += 1
			
		
		
		
		location = low - (low % W_Comb)
		locinv = (location * a) % n
		
		
		
		
	
	
	
	
	
	











def outer_loop(origx,n,filter_loc,filter_est,B2,num,B,W_Comb,Comb_loops,loop_threshold,location_loops,loops, ALG_TYPE):
	
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
		
		x_samp_i, J = inner_loop_locate(origx, n, cur_filter,num, cur_B,a, ai, b,x_samp[i], J,PF_T, BC_T)
		x_samp.append(x_samp_i)
		
		
		if (perform_location):
			if (ALG_TYPE == 1):
				
				#inner_loop_filter_regular(J, n, num, cur_B, a, ai, b, loop_threshold, score, hits, hits_found, G_T)
				hits_found = inner_loop_filter_regular(J, n, num, cur_B, a, ai, b, loop_threshold, score, hits, hits_found)
			else:
				#inner_loop_filter_Comb(J, n, num, cur_B, a, ai, b, loop_threshold, score, hits, hits_found, G_T, Comb_Approved, num_Comb, W_Comb)
				
			
			
		
		
		
		
		
		
	
	
	
	
	std::map<int, complex_t> ans = estimate_values(hits, hits_found, x_samp,  loops, n, permute, B, B2, filter, filter_Est, location_loops);
	
	













































