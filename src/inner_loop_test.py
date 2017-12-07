
import fractions
import math
import cmath
import numpy as np
from scipy.stats import binom
from scipy import fftpack
from scipy.fftpack import fft, fftshift, ifft

import matplotlib.pyplot as plt


from filters import *
from utils import *




def Comb_Filt(origx, n, num, W_Comb):
	# randomly samples W_Comb elements from input array. Then performs fft on it and finds the largest frequencies. 
	# returns Comb_Approved which contains the num largest frequencies.
	if(n%W_Comb):
		
		
		raise Exception('W_Comb is not divisible by N, which algorithm expects')
	
	
	x_sampt = np.zeros(W_Comb, dtype = complex)
	
	sigma = n/W_Comb;
	#offset = int(np.random.rand(1)[0] * sigma)
	offset = int(0.5 * sigma)
	
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
	# inner_loop_locate(origx, n, current_filter, num, current_B, a, ai, b)
	
	if (n % B):
		print("Warning: n is not divisible by B")
	
	x_sampt = np.zeros(B, dtype = complex)
	index=b
	##print("Curr filt t shape: %d" % current_filter.sig_t.shape[0])
	for i in range(current_filter.sig_t.shape[0]):
		
		
		##print("Index: %d" % index)
		
		x_sampt[i%B] += origx[index] * current_filter.sig_t[i]
		index = (index+ai) %n
		
		
	
	plt.plot(abs(x_sampt))
	plt.title("subs t")
	plt.show()
	
	
	x_samp_i = fftpack.fft(x_sampt, B)
	
	samples = np.power( np.absolute(x_samp_i), 2)
	
	#plt.plot(samples)
	#plt.show()
	
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
		
		#for a in range(2):
		#	
		#	values[a] = nth_element(values[a], location)
		#	
		
		
		realv = values[0][location]
		imagv = values[1][location]
		
		##ans.append( ( hits[i], complex(realv, imagv) ) )
		ans[hits[i]] = complex(realv, imagv)
		
	
	
	
	return ans




##def main():

n = 8192
k=10
Bcst_loc =2
Bcst_est =  2
Comb_cst =1
comb_loops =1
est_loops =16
loc_loops =7
threshold_loops =6
tolerance_loc =1e-8
tolerance_est =1e-8



BB_loc = Bcst_loc*np.sqrt( (n*k) /np.log2(n) ) 
BB_est = Bcst_est*np.sqrt( (n*k)/np.log2(n) ) 

lobefrac_loc = 0.5/BB_loc
lobefrac_est = 0.5/BB_est
b_loc = int(1.2*1.1*(n/BB_loc))
b_est = int(1.4*1.1*(n/BB_est))

B_loc = floor_to_pow2(BB_loc)
B_thresh = 2*k
B_est = floor_to_pow2(BB_est)

W_Comb = floor_to_pow2(Comb_cst*n/B_loc)

x_f = np.zeros(n)
large_freq = np.zeros(k, dtype = int)
large_freq[0] = 6711
large_freq[1] = 5050
large_freq[2] = 3456
large_freq[3] = 5493
large_freq[4] = 10
large_freq[5] = 11
large_freq[6] = 2347
large_freq[7] = 100
large_freq[8] = 523
large_freq[9] = 55
x_f[large_freq] = 1.0

plt.plot(x_f)
plt.show()



x = n*fftpack.ifft(x_f, n)
##x = fftpack.fft(x_f, n)



plt.plot(abs(x[:200]))
plt.show()


## add noise to x

x_f = fftpack.fft(x, n)/n



filter_t = chebyshev_window(lobefrac_loc, tolerance_loc)
filter_loc = make_multiple(filter_t, n, b_loc)

filter_t = chebyshev_window(lobefrac_est, tolerance_est)
filter_est = make_multiple(filter_t, n, b_est)

#plt.plot(filter_loc.sig_t)
#plt.title("Loc filter t %d" %n)
#plt.show()

#print("w_loc %d" % filter_loc.sig_t.shape[0])

#for l in range(filter_loc.sig_t.shape[0]):
#	
#	print("sig_t[%d] = %3.8f" %(l, (filter_loc.sig_t[l]*1000.0)))
#	



#plt.plot(np.abs(filter_loc.sig_f))
#plt.title("Loc filter f %d" % n)
#plt.show()


#plt.plot(filter_est.sig_t)
#plt.title("EST filter t %d" % n)
#plt.show()


#plt.plot(np.abs(filter_est.sig_f) )
#plt.title("Est filter f %d" % n)
#plt.show()



origx = x
B2 = B_est
num = B_thresh
B = B_loc
Comb_loops = comb_loops
loop_threshold = threshold_loops
location_loops = loc_loops
loops = loc_loops + est_loops


permute  = np.zeros(loops, dtype=int)
permuteb = np.zeros(loops, dtype=int)


x_samp = []

hits_found = 0

hits = np.zeros(n, dtype=int)
score = np.zeros(n, dtype = int)
num_Comb = num
##Comb_Approved = np.zeros(Comb_loops*num, dtype = int)
Comb_Approved = []



##num_largest = Comb_Filt(origx, n, num, W_Comb)






for i in range(loops):
	
	a = 7773
	b = 0
	
	#while(gcd(a,n) != 1):
	#	
	#	a = np.random.randint(0, n, 1)[0]
	#	
	
	
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
	##for l in range(len(J)):
	##	
	##	print("Loop: %d, J[%d] = %d" % (i, l, J[l] ))
	##	
	#	
	
	
	####x_samp.append(x_samp_i)
	
	
	
	
	if (perform_location):
		hits_found = inner_loop_filter_regular(J, n, num, current_B, a, ai, b, loop_threshold, score, hits, hits_found)
	
	
	
	
	
	

#plt.plot(score[:150])
#plt.title("score")
#plt.show()

#plt.plot(hits[:150])
#plt.title("hits")
#plt.show()

x_samp = np.array(x_samp)

ans = estimate_values(hits, hits_found, x_samp,  loops, n, permute, B, B2, filter_loc, filter_est, location_loops)




num_candidates = len(ans)
candidates = np.zeros((num_candidates, 2))
x_f_Large = np.zeros(n, dtype=complex)
ans_Large = np.zeros(n, dtype = complex)
counter = 0
ERROR = 0.0

for key in sorted(ans.iterkeys()):
	
	value = ans[key]
	
	#print("Key: %d" % key)
	#print("Val: %3.5f " % abs(ans[key]))
	
	candidates[counter][1] = int(key)
	candidates[counter][0] = abs(value)
	counter += 1


for l in range(counter):
	
	key = candidates[l][1]
	value = candidates[l][0]
	#print("candidates[%d][key] = %d" % (l,key))
	#print("candidates[%d][val] = %3.5f " % (l, value))
	


for i in range(k):
	x_f_Large[large_freq[i]]=x_f[large_freq[i]]
	


tmp = np.argpartition(candidates[:,0], num_candidates-k)
candidates = candidates[tmp]

for l in range(k):
	
	key = int(candidates[num_candidates - k + l][1])
	ans_Large[key] = ans[key]
	

for i in range(n):
	ERROR += abs(ans_Large[i] - x_f_Large[i])





#idx_largest = find_largest_indices(k, candidates[:,1])
#keys = np.array(candidates[idx_largest][:,0], dtype=int)
#values = candidates[idx_largest][:,1]

#for key in keys:
#	ans_Large[key] = ans[key]
#	












