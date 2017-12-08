# SparseFFT
# Authors: Andrii Hlyvko, Wentao Zhu
# Fall 2017
A Python implementation of the sparse fourier transform. The code is implemented in Python 2.7. The code requires the numpy, scipy, matplotlib Python packages to be installed
to run. The main code file is the generate_graphs.py. This file will run the experiments and generate all plots.
The arguments to generate_graphs.py are: 
--graph_type     integer 1,2, or 3. Tells which graphs to generate. 1 will run the experiment for various ranges of N input signal sizes.
2 will run the experiment for various K keeping N constant. 3 will keep K and N constant and run the experiment for various SNR.

Example (assuming running from src folder):

python ./generate_graphs --graph_type 1




Note on runtime of the code: The Sparse fourier transform runtime is compared against the scipy.fftpack.fft which was implemented in C using vast
vectorized operations. The runtimes of the fft were faster than the fftw which the paper compared against. Since our implementation was not fully vectorized 
it runs slower than the one in the paper.

The plots folder contains the generated plots. The src folder contains all code. The report folder contains the pdf of the report and the presentation
is in the presentation folder.




