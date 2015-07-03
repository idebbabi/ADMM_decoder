// Copyright 2012 Xishuo Liu, Stark C. Draper, Benjamin Recht
//
// This program is distributed under the terms of the GNU General Public License.
//
// This file is part of ADMM Decoder.
//
//    ADMM Decoder is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ADMM Decoder is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ADMM Decoder.  If not, see <http://www.gnu.org/licenses/>.
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Project:				ADMM Decoder
// Files:				LDPC_Simulator_Example.cpp, LDPC_Class.cpp,
//						LDPC_Simulator_Data_Def.h, LDPC_Class.h,
//						MersenneTwister.h
// Date:				8.30.2012
//
// Author:				Xishuo Liu, xliu94@wisc.edu
// Thanks to:			S. Barman, S. Draper and B. Recht.
//
// Papers:				1. S. Barman, X. Liu, S. Draper and B. Recht,  
//						"Decomposition Methods for Large Scale LP Decoding"
//						http://arxiv.org/abs/1204.0556
//						2. X. Liu, S. Draper and B. Recht,
//						"Suppressing Pseudocodewords by Penalizing the Objective of LP Decoding"
//						IEEE Information Theory Workshop (ITW), 2012. Lausanne: Switzerland, 2012
///////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
This file contains examples for the ADMM decoder classes. 
*/
//#define DEBUG
#include "LDPC_Class.h"
using namespace std;

int main(){
	/*************************************************
		Example 1: Demonstrate different decoders
	**************************************************/
	Simulator ldpcsim(2640, 1320, "marguliscode.txt", "AWGN", 1.6); // simulate the [2640, 1320] "Margulis" code. 
																	// simulate the AWGN channel with Eb/N0 = 1.6dB
	ldpcsim.SetCodeword("codeword2640.txt");						// simulate the codeword stored in file "codeword2640.txt
	
	ldpcsim.SetCommandLineOutput(100, 2);							// create an command line output every 100 simulations. Output all details for simulation
	ldpcsim.SetDecoderParameters(1000, 1e-5, 5.5, 1.9);				// set max iteration = 1000,
																	// ADMM ending tolerance = 1e-5
																	// ADMM \mu = 5.5
																	// ADMM over-relaxation parameter = 1.9
	ldpcsim.SetTargets(0, 100, 0);									// set target number of simulation = 100
	int seed = rand();												// generate a seed
	
	// ADMM LP decoder
	ldpcsim.SetDecoder("ADMMLP");
	ldpcsim.SetChannelSeed(seed);
	ldpcsim.RunSim();
	// ADMM L1 penalized decoder
	ldpcsim.SetDecoder("ADMML1");
	ldpcsim.SetDecoderParameters(1000, 1e-5, 5.5, 1.9, 0.6); // \alpha = 0.6
	ldpcsim.SetChannelSeed(seed);
	ldpcsim.RunSim();
	// ADMM L2 penalized decoder
	ldpcsim.SetDecoder("ADMML2");
	ldpcsim.SetDecoderParameters(1000, 1e-5, 5.5, 1.9, 0.8);// \alpha = 0.8
	ldpcsim.SetChannelSeed(seed);
	ldpcsim.RunSim();
	// ADMM Entropy penalized decoder
	ldpcsim.SetDecoder("ADMMEntropy");
	ldpcsim.SetDecoderParameters(1000, 1e-5, 5.5, 1.9, 0.4, 0.2);// \alpha = 0.4, threshold = 0.2
	ldpcsim.SetChannelSeed(seed);
	ldpcsim.RunSim();
	// non-saturating BP decoder
	ldpcsim.SetDecoder("BP");
	ldpcsim.SetDecoderParameters(1000, 1e-5, 5.5, 1.9, true, 0);
	ldpcsim.SetChannelSeed(seed);
	ldpcsim.RunSim();
	// original BP, saturate at 7
	ldpcsim.SetDecoderParameters(1000, 1e-5, 5.5, 1.9, false, 7);
	ldpcsim.SetChannelSeed(seed);
	ldpcsim.RunSim();

	/*The output from Example 1 should be something like the following:
	TotalSims = 100  Time = 49.14*99%(s)  Total error = 45  Undetected error = 0  ML error = 0  Avg # iteration = 490.49  Avg # iteration correct = 73.6182  Avg # iteration wrong = 1000  Avg decode time = 0.4912(s)  Avg decode time correct = 0.0772727(s)  Avg decode time wrong = 0.997111(s)
	TotalSims = 100  Time = 5.97*98%(s)  Total error = 2  Undetected error = 0  ML error = 0  Avg # iteration = 55.31  Avg # iteration correct = 36.0306  Avg # iteration wrong = 1000  Avg decode time = 0.059(s)  Avg decode time correct = 0.0381633(s)  Avg decode time wrong = 1.08(s)
	TotalSims = 100  Time = 4.35*99%(s)  Total error = 1  Undetected error = 0  ML error = 0  Avg # iteration = 40.63  Avg # iteration correct = 30.9394  Avg # iteration wrong = 1000  Avg decode time = 0.0433(s)  Avg decode time correct = 0.0327273(s)  Avg decode time wrong = 1.09(s)
	TotalSims = 100  Time = 6.15*99%(s)  Total error = 1  Undetected error = 0  ML error = 0  Avg # iteration = 38.69  Avg # iteration correct = 28.9798  Avg # iteration wrong = 1000  Avg decode time = 0.0612(s)  Avg decode time correct = 0.0458586(s)  Avg decode time wrong = 1.58(s)
	TotalSims = 100  Time = 15.3*99%(s)  Total error = 1  Undetected error = 0  ML error = 0  Avg # iteration = 28.2  Avg # iteration correct = 18.3838  Avg # iteration wrong = 1000  Avg decode time = 0.1525(s)  Avg decode time correct = 0.0987879(s)  Avg decode time wrong = 5.47(s)
	TotalSims = 100  Time = 7.74*99%(s)  Total error = 1  Undetected error = 0  ML error = 0  Avg # iteration = 26.97  Avg # iteration correct = 17.1414  Avg # iteration wrong = 1000  Avg decode time = 0.0772(s)  Avg decode time correct = 0.0491919(s)  Avg decode time wrong = 2.85(s)
	*/



	/******************************************************************************
		Example 2: Demonstrate effects of penalty parameter and log it on files
	*******************************************************************************/
	ofstream ofile("results.txt");	// set output files
	ofstream ofile2("results_brief.txt");
	
	for(double alpha = 0.2; alpha < 1.6; alpha += 0.2)
	{
		ldpcsim.SetDecoder("ADMML1");
		ldpcsim.SetDecoderParameters(50000, 1e-5, 5.5, 1.9, alpha); // the maximum number of iterations is 50000 (which is really big number..)
		ldpcsim.SetTargets(200, 0, 0); // set target error = 200
		ldpcsim.RunSim();
		ofile<<"alpha = "<<alpha<<endl;
		ofile<<"TotalSims = "<<ldpcsim.GetTotalSims()
					<<"  Time = "<<(double)ldpcsim.GetTotalExeTime()/CLOCKS_PER_SEC<<"*"<<(ldpcsim.GetTotalDecodingTime()*100)/ldpcsim.GetTotalExeTime()<<"%(s)"
					<<"  Total error = "<<ldpcsim.GetTotalErrors()
					<<"  Undetected error = "<<ldpcsim.GetTotalUndetectedErrors()
					<<"  ML error = "<<ldpcsim.GetTotalMLErrors()
					<<"  Avg # iteration = "<<ldpcsim.GetTotalIterations()/(double)(ldpcsim.GetTotalSims())
					<<"  Avg # iteration correct = "<<ldpcsim.GetTotalCorrectIterations()/(double)(ldpcsim.GetTotalSims() - ldpcsim.GetTotalErrors())
					<<"  Avg # iteration wrong = "<<ldpcsim.GetTotalWrongIterations()/(double)(ldpcsim.GetTotalErrors())
					<<"  Avg decode time = "<<(double)ldpcsim.GetTotalDecodingTime()/CLOCKS_PER_SEC/ (double)(ldpcsim.GetTotalSims()) <<"(s)"
					<<"  Avg decode time correct = "<<(double)ldpcsim.GetTotalCorrectDecodingTime()/CLOCKS_PER_SEC/ (double)(ldpcsim.GetTotalSims() - ldpcsim.GetTotalErrors()) <<"(s)"
					<<"  Avg decode time wrong = "<<(double)ldpcsim.GetTotalWrongDecodingTime()/CLOCKS_PER_SEC/ (double)(ldpcsim.GetTotalErrors()) <<"(s)"
					<<endl<<endl;
		ofile2<<alpha;
		ofile2<<"  "<<ldpcsim.GetTotalSims()
					<<"  "<<(double)ldpcsim.GetTotalExeTime()/CLOCKS_PER_SEC<<"*"<<(ldpcsim.GetTotalDecodingTime()*100)/ldpcsim.GetTotalExeTime()<<"%"
					<<"  "<<ldpcsim.GetTotalErrors()
					<<"  "<<ldpcsim.GetTotalUndetectedErrors()
					<<"  "<<ldpcsim.GetTotalMLErrors()
					<<"  "<<ldpcsim.GetTotalIterations()/(double)(ldpcsim.GetTotalSims())
					<<"  "<<ldpcsim.GetTotalCorrectIterations()/(double)(ldpcsim.GetTotalSims() - ldpcsim.GetTotalErrors())
					<<"  "<<ldpcsim.GetTotalWrongIterations()/(double)(ldpcsim.GetTotalErrors())
					<<"  "<<(double)ldpcsim.GetTotalDecodingTime()/CLOCKS_PER_SEC/ (double)(ldpcsim.GetTotalSims())
					<<"  "<<(double)ldpcsim.GetTotalCorrectDecodingTime()/CLOCKS_PER_SEC/ (double)(ldpcsim.GetTotalSims() - ldpcsim.GetTotalErrors()) 
					<<"  "<<(double)ldpcsim.GetTotalWrongDecodingTime()/CLOCKS_PER_SEC/ (double)(ldpcsim.GetTotalErrors()) 
					<<endl;
	}


	/* There should be two files logging the results.
		==========|results.txt|==========
		
		alpha = 0
		TotalSims = 434  Time = 8902.05*99%(s)  Total error = 200  Undetected error = 0  ML error = 0  Avg # iteration = 23148.5  Avg # iteration correct = 198.517  Avg # iteration wrong = 50000  Avg decode time = 20.5112(s)  Avg decode time correct = 0.191154(s)  Avg decode time wrong = 44.2858(s)

		alpha = 0.2
		TotalSims = 2839  Time = 9575.68*99%(s)  Total error = 200  Undetected error = 0  ML error = 0  Avg # iteration = 3584.66  Avg # iteration correct = 67.0159  Avg # iteration wrong = 50000  Avg decode time = 3.3725(s)  Avg decode time correct = 0.0687192(s)  Avg decode time wrong = 46.9659(s)

		alpha = 0.4
		TotalSims = 7460  Time = 10025.6*99%(s)  Total error = 200  Undetected error = 0  ML error = 0  Avg # iteration = 1410.07  Avg # iteration correct = 71.5015  Avg # iteration wrong = 50000  Avg decode time = 1.34351(s)  Avg decode time correct = 0.0732631(s)  Avg decode time wrong = 47.4533(s)

		alpha = 0.6
		TotalSims = 12019  Time = 11227.8*99%(s)  Total error = 200  Undetected error = 0  ML error = 0  Avg # iteration = 958.432  Avg # iteration correct = 128.555  Avg # iteration wrong = 50000  Avg decode time = 0.933796(s)  Avg decode time correct = 0.131289(s)  Avg decode time wrong = 48.358(s)

		alpha = 0.8
		TotalSims = 26797  Time = 16406.2*99%(s)  Total error = 200  Undetected error = 0  ML error = 0  Avg # iteration = 601.754  Avg # iteration correct = 230.296  Avg # iteration wrong = 50000  Avg decode time = 0.611842(s)  Avg decode time correct = 0.235946(s)  Avg decode time wrong = 50.6003(s)

		alpha = 1
		TotalSims = 36616  Time = 20182.9*99%(s)  Total error = 200  Undetected error = 0  ML error = 0  Avg # iteration = 540.549  Avg # iteration correct = 268.913  Avg # iteration wrong = 50000  Avg decode time = 0.550811(s)  Avg decode time correct = 0.273554(s)  Avg decode time wrong = 51.0336(s)

		alpha = 1.2
		TotalSims = 38594  Time = 23427.9*99%(s)  Total error = 200  Undetected error = 0  ML error = 0  Avg # iteration = 598.451  Avg # iteration correct = 341.111  Avg # iteration wrong = 50000  Avg decode time = 0.606615(s)  Avg decode time correct = 0.345096(s)  Avg decode time wrong = 50.8104(s)

		alpha = 1.4
		TotalSims = 25835  Time = 22645.6*99%(s)  Total error = 200  Undetected error = 0  ML error = 0  Avg # iteration = 876.518  Avg # iteration correct = 493.265  Avg # iteration wrong = 50000  Avg decode time = 0.876115(s)  Avg decode time correct = 0.491814(s)  Avg decode time wrong = 50.1339(s)

		alpha = 1.6
		TotalSims = 18381  Time = 21088.5*99%(s)  Total error = 200  Undetected error = 0  ML error = 0  Avg # iteration = 1157.86  Avg # iteration correct = 620.575  Avg # iteration wrong = 50000  Avg decode time = 1.1469(s)  Avg decode time correct = 0.612875(s)  Avg decode time wrong = 49.6921(s)

		============|results_brief.txt|===========

		0  434  8902.05*99%  200  0  0  23148.5  198.517  50000  20.5112  0.191154  44.2858
		0.2  2839  9575.68*99%  200  0  0  3584.66  67.0159  50000  3.3725  0.0687192  46.9659
		0.4  7460  10025.6*99%  200  0  0  1410.07  71.5015  50000  1.34351  0.0732631  47.4533
		0.6  12019  11227.8*99%  200  0  0  958.432  128.555  50000  0.933796  0.131289  48.358
		0.8  26797  16406.2*99%  200  0  0  601.754  230.296  50000  0.611842  0.235946  50.6003
		1  36616  20182.9*99%  200  0  0  540.549  268.913  50000  0.550811  0.273554  51.0336
		1.2  38594  23427.9*99%  200  0  0  598.451  341.111  50000  0.606615  0.345096  50.8104
		1.4  25835  22645.6*99%  200  0  0  876.518  493.265  50000  0.876115  0.491814  50.1339
		1.6  18381  21088.5*99%  200  0  0  1157.86  620.575  50000  1.1469  0.612875  49.6921
	*/
	return 0;
}