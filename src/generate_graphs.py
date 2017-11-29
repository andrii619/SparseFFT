
import argparse


# Default parameter values that can be changed via command line

GRAPH_TYPE = 1
REPETITIONS = 10
VERBOSE = True
ALG_TYPE = 1



BCST_LOC = 1
BCST_EST = 1
COMB_CST = 2
LOC_LOOPS = 4
EST_LOOPS = 16
THRESHOLD_LOOPS = 3
COMB_LOOPS = 1
TOLERANCE_LOC = 1e-8
TOLERANCE_EST = 1e-8
SNR = 100





def main():
	
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
		
		
		
	
	
	##N_vec = np.array([8192, 16384, 32768, 65536, 131072, 262144,  524288, 1048576, 2097152, 4194304, 8388608, 16777216])
	N_vec = np.array([8192, 16384, 32768])
	K_vec = np.array([50, 100, 200, 500, 1000, 2000, 2500, 4000])
	SNR_vec[14] = np.array([-20, -10, -7, -3, 0, 3, 7, 10, 20, 30, 40, 50, 60, 120])
	
	
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
			
			experiment_parameters = get_expermient_vs_K_parameters(n, ALG_TYPE)
			
			
			
		else:
			n = 4194304
			k = 50
			SNR = SNR_vec[i]; 
			experiment_parameters = get_expermient_vs_N_parameters(n, WITH_COMB)
			
			
			
			
		
		
		if(experiment_parameters['Bcst_loc']):
			BCST_LOC = experiment_parameters['Bcst_loc']
		
		if(experiment_parameters['Bcst_est']):
			BCST_EST = experiment_parameters['Bcst_est']
		
		if(experiment_parameters['Comb_cst']):
			COMB_CST = experiment_parameters['Comb_cst']
		
		if(experiment_parameters['Comb_loops']):
			COMB_LOOPS = experiment_parameters['Comb_loops']
		
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
		
		
		
		for j in range(REPETITIONS):
			
			
			
			x, x_f = generate_random_signal(n, k)
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
		
		
		
		
		
		
		
		
		
		
	
	
	



if __name__ == "__main__":
	
	
	main()











































































































































































































