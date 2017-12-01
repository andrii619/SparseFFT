
import argparse
import numpy as np
from utils import *

# Default parameter values that can be changed via command line

GRAPH_TYPE = 1
REPETITIONS = 1
VERBOSE = True
ALG_TYPE = 1



#BCST_LOC = 1
#BCST_EST = 1
#COMB_CST = 2
#LOC_LOOPS = 4
#EST_LOOPS = 16
#THRESHOLD_LOOPS = 3
#COMB_LOOPS = 1
#TOLERANCE_LOC = 1e-8
#TOLERANCE_EST = 1e-8
#SNR = 100




def run_experiment(x,x_f,large_freq,k,n,lobefrac_loc,tolerance_loc,b_loc,B_loc,B_thresh,loc_loops,threshold_loops,lobefrac_est,tolerance_est,
						b_est,B_est,est_loops,W_Comb,comb_loops):
	
	
	print("0: n: %d" % n)
	print("1:lobefrac_loc %3.3f" % lobefrac_loc)
	print("2: tolerance_loc %3.3f" % tolerance_loc)
	print("3: b_loc %d" % b_loc)
	print("4: B_loc %d" % B_loc)
	print("5: B_thresh %d" % B_thresh)
	print("6: loc_loops %d" % loc_loops)
	print("7: threshold_loops %d" % threshold_loops)
	print("8: lobefrac_est %3.3f" % lobefrac_est)
	print("9: tolerance_est %3.3f" % tolerance_est)
	print("10: b_est %d" % b_est)
	print("11: B_est %d" % B_est)
	print("12: est_loops %d" % est_loops)
	print("13: W_Comb %d" % W_Comb)
	print("14: Comb_loops %d" % comb_loops)
	print("15: k %d" % k)
	
	
	
	
	
	
	
	
	



def main():
	global GRAPH_TYPE, REPETITIONS, ALG_TYPE, VERBOSE
	
	
	parser = argparse.ArgumentParser(description='Sparse fft')
	##parser.add_argument('-h','--help', help='Print Help', required=False)
	parser.add_argument('-g','--graph_type', help='Type of the graph to plot', type=int, required=False)
	parser.add_argument('-r','--repetitions', help='Number of repetitions', type=int, required=False)
	parser.add_argument('-v','--verbose', help='Verbose mode', type=bool, required=False)
	parser.add_argument('-a','--alg_type', help='sFFT 1.0 or sFFT 2.0', type = int, required=False)
	
	args = vars(parser.parse_args())
	
	if(args['graph_type']):
		
		g = args['graph_type']
		if(g >= 1 and g <= 3):
			GRAPH_TYPE = args['graph_type']
			print("Program arg graph_type: %d" % GRAPH_TYPE)
		else:
			print("graph_type: Invalid option")
		
		
	
	if(args['repetitions']):
		
		r = args['repetitions']
		if(r>=1):
			REPETITIONS = args['repetitions']
			print("Program arg repetirions: %d" % REPETITIONS)
			
		else:
			print("repetitions: Invalid Option")
		
		
	
	
	if(args['verbose']):
		
		
		VERBOSE = args['verbose']
		print(VERBOSE)
		
		
	
	
	if(args['alg_type']):
		
		a = args['alg_type']
		if(a==1 or a== 2):
			
			ALG_TYPE = args['alg_type']
			
			
			
		else:
			print("alg_type: Invalid Option")
		
		
		
	
	
	N_vec = np.array([8192, 16384, 32768, 65536, 131072, 262144,  524288, 1048576, 2097152, 4194304, 8388608, 16777216], dtype=np.int)
	##N_vec = np.array([8192, 16384, 32768])
	K_vec = np.array([50, 100, 200, 500, 1000, 2000, 2500, 4000], dtype=np.int)
	SNR_vec = np.array([-20, -10, -7, -3, 0, 3, 7, 10, 20, 30, 40, 50, 60, 120])
	
	
	length = 1
	if(GRAPH_TYPE == 1):
		length = len(N_vec)
	elif(GRAPH_TYPE == 2):
		length = len(K_vec)
	else:
		length = len(SNR_vec)
	
	
	sFFT_times = np.zeros(length)
	fft_times = np.zeros(length)
	sFFT_errors = np.zeros(length)
	
	
	for i in range(length):
		
		
		if(GRAPH_TYPE == 1):
			
			n = N_vec[i]
			k = 50
			
			experiment_parameters = get_expermient_vs_N_parameters(n, ALG_TYPE)
			
			###Bcst_loc, Bcst_est, Comb_cst, comb_loops, est_loops, loc_loops, threshold_loops, tolerance_loc,tolerance_est =1e-8;
			
			
		elif(GRAPH_TYPE == 2):
			n = 4194304;
			k = K_vec[i]
			
			experiment_parameters = get_expermient_vs_K_parameters(k, ALG_TYPE)
			
			
			
		else:
			n = 4194304
			k = 50
			SNR = SNR_vec[i]; 
			experiment_parameters = get_expermient_vs_N_parameters(n, ALG_TYPE)
			
			
			
			
		
		
		if(experiment_parameters['Bcst_loc']):
			BCST_LOC = experiment_parameters['Bcst_loc']
		
		if(experiment_parameters['Bcst_est']):
			BCST_EST = experiment_parameters['Bcst_est']
		
		if(experiment_parameters['Comb_cst']):
			COMB_CST = experiment_parameters['Comb_cst']
		
		if(experiment_parameters['comb_loops']):
			COMB_LOOPS = experiment_parameters['comb_loops']
		
		if(experiment_parameters['est_loops']):
			EST_LOOPS = experiment_parameters['est_loops']
		
		if(experiment_parameters['loc_loops']):
			LOC_LOOPS = experiment_parameters['loc_loops']
		
		if(experiment_parameters['threshold_loops']):
			THRESHOLD_LOOPS = experiment_parameters['threshold_loops']
		
		if(experiment_parameters['tolerance_loc']):
			TOLERANCE_LOC = experiment_parameters['tolerance_loc']
		
		if(experiment_parameters['tolerance_est']):
			TOLERANCE_EST = experiment_parameters['tolerance_est']
		
		
		BB_loc = BCST_LOC*np.sqrt( (n*k) /np.log2(n) ) 
		BB_est = BCST_EST*np.sqrt( (n*k)/np.log2(n) ) 
		
		
		print("BB_loc: %4.5f" % BB_loc)
		print("BB_est: %4.5f" % BB_est)
		
		
		LOBEFRAC_LOC = 0.5/BB_loc
		LOBEFRAC_EST = 0.5/BB_est
		b_loc = int(1.2*1.1*(n/BB_loc))
		b_est = int(1.4*1.1*(n/BB_est))
		
		B_loc = floor_to_pow2(BB_loc)
		B_thresh = 2*k
		B_est = floor_to_pow2(BB_est)
		
		W_Comb = floor_to_pow2(COMB_CST*n/B_loc)
		
		for j in range(REPETITIONS):
			
			x, x_f, large_freq = generate_random_signal(n, k)
			
			run_experiment(x,x_f,large_freq,k,n,LOBEFRAC_LOC,TOLERANCE_LOC,b_loc,B_loc,B_thresh,LOC_LOOPS,THRESHOLD_LOOPS,LOBEFRAC_EST,TOLERANCE_EST,
						b_est,B_est,EST_LOOPS,W_Comb,COMB_LOOPS)
			
			
			
			
			
			
			
			
			
			
			
			
		
		
		
		
		
		
		
		
		
		
	
	
	



if __name__ == "__main__":
	
	
	main()











































































































































































































