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
// Name:				ADMM Decoder
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
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
This file contains methods for the ADMM decoder classes. 
*/

#include "LDPC_Class.h"

//Noisy sequence Class
NoisySeq::NoisySeq(int blocklength)
{
	Blocklength = blocklength;
	NoisySequence = new double[blocklength];
}
NoisySeq::~NoisySeq()
{
#ifdef DEBUG
cout<<"~NoisySeq()"<<endl;
#endif
	delete [] NoisySequence;
}

//Channel Class
void Channel::SetChannel(int blocklength, string channel, double parameter)
{
	if(channel == "AWGN")
	{
		mSTD = parameter;
	}
	else if (channel == "BSC")
	{
		mCrossoverProb = parameter;
	}
	else
	{
		cout<<"channel not supported"<<endl;
		exit(0);
	}
		
	mBlocklength = blocklength;
	mChannelName = channel;
}
void Channel::GenerateOutput(int input[], NoisySeq &ChannelOutput)
{
	
	if(mChannelName == "AWGN")
	{
		ChannelOutput.ChannelParameter = mSTD;
		mGenerateAWGN(input, ChannelOutput);
	}
	else if (mChannelName == "BSC")
	{
		ChannelOutput.ChannelParameter = mCrossoverProb;
		mGenerateBitFlip(input, ChannelOutput);
	}
	else
	{
		cout<<"channel not supported"<<endl;
		exit(0);
	}
}
void Channel::mGenerateAWGN(int input[], NoisySeq &ChannelOutput)// generate noisy sequence for AWGN
{
	double std_check = 0;
	for (int i = 0; i < mBlocklength; i++)
	{
		int transmittedbit = 0;
		transmittedbit = 2 * input[i] - 1;

		double receivedbit = channelrand.randNorm(transmittedbit, mSTD); // assume all zero codeword is sent
		ChannelOutput.NoisySequence[i] = receivedbit;
		//_LogLikelihoodRatio[i] = ((receivedbit - 1)*(receivedbit - 1) - (receivedbit + 1)*(receivedbit + 1))/2.0/std/std;
		std_check += (receivedbit - transmittedbit)*(receivedbit - transmittedbit);
	}
#ifdef DEBUG
	cout<<"std = "<<sqrt(std_check/(double)mBlocklength)<<" real std = "<<mSTD<<endl;
#endif
}
void Channel::mGenerateBitFlip(int input[], NoisySeq &ChannelOutput)// generate noisy sequence for BSC
{
	int total_bit_flips = 0;
	for (int i = 0; i < mBlocklength; i++)
	{
		double rand_num = channelrand.rand();
		if(rand_num > mCrossoverProb)
		{
			ChannelOutput.NoisySequence[i] = input[i];
		}
		else 
		{
			ChannelOutput.NoisySequence[i] = 1 - input[i];
			total_bit_flips++;
		}
	}
#ifdef DEBUG
	cout<<"p = "<<(double)total_bit_flips/mBlocklength<<" real p = "<<mCrossoverProb<<endl;
#endif
}
void Channel::SetSeed(int seed){channelrand.seed(seed);}


// Decoder Class
Decoder::Decoder(string FileName,int nChecks, int BlockLength)
{
	// Learn and read parity check matrix
	mPCMatrixFileName = FileName;
	mNChecks = nChecks;
	mBlocklength = BlockLength;
	mPCheckMapSize = 0;
	CheckDegree = new int[mNChecks];
	VariableDegree = new int[mBlocklength];
	OutputFromDecoder = new double[mBlocklength];
	_LogLikelihoodRatio = new double [mBlocklength];
	u0 = new double [mBlocklength];

	mSetDefaultParameters();
	// learn the structure of parity check matrix: degree, # edges
	mLearnParityCheckMatrix(); 
	mPCheckMap = new TRIPLE[mPCheckMapSize];
	mParityCheckMatrix = new TRIPLE[mPCheckMapSize];
	u = new MESSAGE<double> [mPCheckMapSize];
	v = new MESSAGE<double> [mPCheckMapSize];
	mReadParityCheckMatrix();
	warmstartchecks = new double[mPCheckMapSize];
	warmstartsave = new double[mPCheckMapSize];

	RateOfCode = (double)(mBlocklength -mNChecks)/mBlocklength;
}
Decoder::~Decoder()
{
	delete [] u;
	delete [] v;
	delete [] u0;
	delete [] CheckDegree;
	delete [] VariableDegree;
	delete [] warmstartsave;
	delete [] warmstartchecks;
	delete [] mPCheckMap;
	delete [] mParityCheckMatrix;
	delete [] OutputFromDecoder;
	delete [] _LogLikelihoodRatio;
}

//Decodersettings
void Decoder::SetParameters(int mxIt, double feas,  double p_mu, double p_rho)
{
	maxIteration = mxIt;
	para_end_feas = feas;
	para_mu = p_mu;
	para_rho = p_rho;
}
void Decoder::SetBSCnormalized(bool normal)
{
	if(normal)
		normalizeBSC = 1;
	else
		normalizeBSC = 0;
}
void Decoder::SetBPSaturation(bool nonsat, double value)
{
	BP_non_sat = nonsat;
	BP_sat_value = value;
}
void Decoder::SetDecoder(string decoder)
{
	mRealDecoderName = decoder;
	if(decoder == "BP")
	{
		DecoderName = decoder;
		mInitializeMessages();
	}
	else if(decoder == "ADMMLP")
	{
		DecoderName = "ADMM";
		enumADMM = LinearProgram;
	}
	else if(decoder == "ADMML1")
	{
		DecoderName = "ADMM";
		enumADMM = L1Penalty;
	}
	else if(decoder == "ADMML2")
	{
		DecoderName = "ADMM";
		enumADMM = L2Penalty;
	}
	else if(decoder == "ADMMEntropy")
	{
		DecoderName = "ADMM";
		enumADMM = EntropyPenalty;
	}
	else if(decoder == "MSA")
	{
		DecoderName = "MSA";
	}
	else
	{
		cout<<"decoder not supported"<<endl;
		exit(0);
	}
}

// Decoder initializations
void Decoder::mSetDefaultParameters()
{
	// decoder settings
	maxIteration = 1000;
	para_end_feas = 1e-6;
	para_mu = 5;
	para_rho= 1.8;

	alpha = 1;
	for(int i = 0; i < mBlocklength; i++)
	{
		VariableDegree[i] = 0;
	}
	for(int i = 0; i < mNChecks; i++)
	{
		CheckDegree[i] = 0;
	}
	normalizeBSC = 0.0;
}
void Decoder::mLearnParityCheckMatrix()// learn the degree of the parity check matrix
{
	ifstream myfile (mPCMatrixFileName.c_str());
	int tempvalue = 0;
	int i = 0;
	string line;
	if (myfile.is_open())
	{
		//cout<<"here 2.1"<<endl;
		while(getline(myfile, line)) 
		{
		   istringstream is(line);
		   int curr_degree = 0;
		   while( is >> tempvalue ) 
		   {
			   
			   VariableDegree[tempvalue]++;
			   curr_degree++;
			   mPCheckMapSize++;
		   }
		   CheckDegree[i] = curr_degree;
		   //cout<<"here 2.11-"<<i<<endl;
		   i++;
		}
		//cout<<"here 2.2"<<endl;
		if (myfile.is_open())
			myfile.close();
		//cout<<"here 2.3"<<endl;

	}
	else 
	{
		cout << "Unable to open file"; 
		exit(0);
	}
}
void Decoder::mReadParityCheckMatrix()// read the parity check matrix, initialize message passing
{
	ifstream myfile (mPCMatrixFileName.c_str());
	int tempvalue = 0;
	int i = 0;
	int count = 0;
	string line;
	//cout<<"here 3.1"<<endl;
	if (myfile.is_open())
	{
		while(getline(myfile, line)) 
		{
			istringstream is(line);
			while( is >> tempvalue ) 
			{
				mParityCheckMatrix[count].col = tempvalue;
				mParityCheckMatrix[count].row = i;

				mPCheckMap[count].col = tempvalue;
				mPCheckMap[count].row = count;
			
				u[count].SetDeg(CheckDegree[i]);

				u[count].Mcol = mParityCheckMatrix[count].col;
				u[count].Mrow = mParityCheckMatrix[count].row;
				u[count].Mprob = 0;

				v[count].SetDeg(VariableDegree[tempvalue]);
				v[count].Mcol = mParityCheckMatrix[count].col;
				v[count].Mrow = mParityCheckMatrix[count].row;
				v[count].Mprob = 0;
			
				count++;
			}
			i++;
		}
	}
	if (myfile.is_open())
		myfile.close();
}
void Decoder::mInitializeMessages() 
{
	for (int i = 0; i < mPCheckMapSize; i++)
	{
		int countu = 0, countv = 0;
		for (int j = 0; j < mPCheckMapSize; j++)
		{
			if (u[i].Mrow == u[j].Mrow && i!=j)
			{
				u[i].i[countu] = j;
				countu++;
			}
			if (v[i].Mcol == v[j].Mcol && i!=j)
			{
				v[i].i[countv] = j;
				countv++;
			}
			if(countu == u[i].degree - 1 && countv == v[i].degree - 1)
				break;
		}
		countu = 0; countv = 0;
			//cout<<(double)i/mPCheckMapSize *100<<endl;
	}
}

// Decode functions
void Decoder::mGenerateLLR(NoisySeq &channeloutput)
{
	if(ChannelName == "AWGN")
	{
		double std = channeloutput.ChannelParameter;
		for(int i = 0; i < mBlocklength; i++)
		{
			double receivedbit = channeloutput.NoisySequence[i];
			_LogLikelihoodRatio[i] = ((receivedbit - 1)*(receivedbit - 1) - (receivedbit + 1)*(receivedbit + 1))/2.0/std/std;
		}
	}
	else if(ChannelName == "BSC")
	{
		double CrossOverProb = channeloutput.ChannelParameter;
		for(int i = 0; i < mBlocklength; i++)
		{
			int receivedbit = floor(channeloutput.NoisySequence[i] + 0.5);
			if(receivedbit == 0)
				_LogLikelihoodRatio[i] = - log(CrossOverProb/(1 - CrossOverProb));
			else
				_LogLikelihoodRatio[i] = log(CrossOverProb/(1 - CrossOverProb));
		}
	}
	else
	{
		cout<<"channel not supported"<<endl;
	}
}
void Decoder::Decode(NoisySeq &channeloutput, DecodedSeq &decoderesult)
{
	mExeTime = 0;
	mAlgorithmConverge = false;
	mValidCodeword = false;

	mGenerateLLR(channeloutput);
	
	clock_t tStart = clock();
	if(DecoderName == "ADMM")
		ADMMDecoder();
	else if(DecoderName == "BP")
	{
# ifdef DEBUG
		cout<<"BP"<<endl;
#endif
		BPDecoder();
	}
	else if(DecoderName == "MSA")
	{
# ifdef DEBUG
		cout<<"MSA"<<endl;
#endif
		MSADecoder();
	}
	else
	{
		cout<<"Decoder not supported"<<endl;
		exit(0);
	}
	mExeTime = clock() - tStart;
	decoderesult.ExeTime = mExeTime;
	decoderesult.Iteration = mIteration;
	decoderesult.ValidCodeword = mValidCodeword;
	decoderesult.AlgorithmConverge = mAlgorithmConverge;
	decoderesult.Iteration = mIteration;

	for(int i = 0; i < mBlocklength; i++)
	{
		if(my_isnan(OutputFromDecoder[i]))
			decoderesult.IsNan = true;
		if(my_isinf(OutputFromDecoder[i]))
			decoderesult.IsInf = true;
		decoderesult.SoftInfo[i] = OutputFromDecoder[i];
		
		if(DecoderName == "ADMM")
		{
			decoderesult.HardDecision[i] = floor(OutputFromDecoder[i] + 0.5);
		}
		else // DecoderName == "BP" or "MSA"
		{
			if(OutputFromDecoder[i] < 0)
				decoderesult.HardDecision[i] = 1;
			else
				decoderesult.HardDecision[i] = 0;
		}
	}
	//cout<<"done?"<<endl;
}
//! comparison function used in sorting
inline bool sort_compare_de(const NODE & a, const NODE & b){return (a.value > b.value);}

double *Decoder::mProjectionPPd(NODE v[], int length)
{
	double zero_tol_offset = 1e-10;
	double *results = new double[length];
	bool AllZero = true, AllOne = true;
	for(int i = 0; i < length; i++)
	{
		if(v[i].value > 0)
			AllZero = false;
		if(v[i].value <= 1)
			AllOne = false;
	}
	if(AllZero) // exit if the vector has all negative values, which means that the projection is all zero vector
	{
		for(int i = 0;i < length; i++)
			results[i] = 0;
		return results;
	}
	if (AllOne && (length%2 == 0)) // exit if the vector should be all one vector.
	{
		for(int i = 0;i < length; i++)
			results[i] = 1;
		return results;
	}

	NODE *zSort = new NODE[length]; // vector that save the sorted input vector
	for(int i = 0;i < length; i++) // keep the original indices while sorting, 
	{
		zSort[i].index = v[i].index;
		zSort[i].value = v[i].value;
	}
	sort(zSort, zSort + length, sort_compare_de);//sort input vector, decreasing order

	
	int clip_idx = 0, zero_idx= 0;
	double constituent = 0;
	NODE *zClip = new NODE[length]; 
	for(int i = 0;i < length; i++)// project on the [0,1]^d cube
	{
		zClip[i].value = min(max(zSort[i].value, 0.0),1.0);
		zClip[i].index = zSort[i].index;
		constituent += zClip[i].value;
	} 
	int r = (int)floor(constituent);
	if (r & 1)
		r--;
	// calculate constituent parity

	// calculate sum_Clip = $f_r^T z$
	double sum_Clip = 0;
	for(int i = 0; i < r+1; i++)
		sum_Clip += zClip[i].value;
	for(int i = r + 1; i < length; i++)
		sum_Clip -= zClip[i].value;
	

	if (sum_Clip <= r) // then done, return projection, beta = 0
	{
		for(int i = 0; i < length; i++)
			results[zClip[i].index] = zClip[i].value;
		delete [] zSort; delete [] zClip;
		return results;
	}

	double beta = 0;
	double beta_max = 0;
	if (r + 2 <= length)
		beta_max = (zSort[r].value - zSort[r+1].value)/2; // assign beta_max
	else
		beta_max = zSort[r].value;

	NODE *zBetaRep = new NODE[length]; 

	// merge beta, save sort
	int left_idx = r, right_idx = r + 1;
	int count_idx = 0;
	while(count_idx < length)
	{
		double temp_a, temp_b;
		if(left_idx < 0)
		{
			while(count_idx < length)
			{
				zBetaRep[count_idx].index = right_idx;
				zBetaRep[count_idx].value = -zSort[right_idx].value;
				right_idx++;count_idx++;
			}
			break;
		}
		if(right_idx >= length)
		{
			while(count_idx < length)
			{
				zBetaRep[count_idx].index = left_idx;
				zBetaRep[count_idx].value = zSort[left_idx].value - 1;
				left_idx--;count_idx++;
			}
			break;
		}
		temp_a = zSort[left_idx].value - 1; temp_b = -zSort[right_idx].value;
		if(temp_a > temp_b || left_idx < 0)
		{
			zBetaRep[count_idx].index = right_idx;
			zBetaRep[count_idx].value = temp_b;
			right_idx++;count_idx++;
		}
		else
		{
			zBetaRep[count_idx].index = left_idx;
			zBetaRep[count_idx].value = temp_a;
			left_idx--;count_idx++;
		}
	}
	


	for(int i = 0; i < length; i++)
	{
		if (zSort[i].value > 1)
			clip_idx++;
		if (zSort[i].value >= 0 - zero_tol_offset)
			zero_idx++;
	}
	clip_idx--;
	bool first_positive = true, below_beta_max = true;
	int idx_start = 0, idx_end = 0;
	for(int i = 0;i < length; i++)
	{
		if(zBetaRep[i].value < 0 + zero_tol_offset && first_positive)
		{
			idx_start++;
		}
		if(zBetaRep[i].value < beta_max && below_beta_max)
		{
			idx_end++;
		}
	}
	idx_end--;
	double active_sum = 0;
	for(int i = 0;i < length; i++)
	{
		if(i > clip_idx && i <= r)
			active_sum += zSort[i].value;
		if(i > r && i < zero_idx)
			active_sum -= zSort[i].value;
	}
	double total_sum = 0;
	total_sum = active_sum + clip_idx + 1;

	int previous_clip_idx, previous_zero_idx;
	double previous_active_sum;
	bool change_pre = true;
	previous_clip_idx = clip_idx;
	previous_zero_idx = zero_idx;
	previous_active_sum = active_sum;
	for(int i = idx_start; i <= idx_end; i++)
	{
		if(change_pre)
		{
			// save previous things
			previous_clip_idx = clip_idx;
			previous_zero_idx = zero_idx;
			previous_active_sum = active_sum;
		}
		change_pre = false;

		beta = zBetaRep[i].value;
		if(zBetaRep[i].index <= r)
		{
			clip_idx--;
			active_sum += zSort[zBetaRep[i].index].value;
		}
		else
		{
			zero_idx++;
			active_sum -= zSort[zBetaRep[i].index].value;
		}
		if (i < length - 1)
		{
			if (beta != zBetaRep[i + 1].value)
			{
				total_sum = (clip_idx+1) + active_sum - beta * (zero_idx - clip_idx - 1);
				change_pre = true;
				if(total_sum < r)
					break;
			}
			else
			{
				if(0)
					cout<<"here here"<<endl;
			}
		}
		else if (i == length - 1)
		{
			total_sum = (clip_idx + 1)  + active_sum - beta * (zero_idx - clip_idx - 1);
			change_pre = true;
		}
	}
	if (total_sum > r)
	{
		beta = -(r - clip_idx - 1 - active_sum)/(zero_idx - clip_idx - 1);
	}
	else
	{
		beta = -(r - previous_clip_idx - 1 - previous_active_sum)/(previous_zero_idx - previous_clip_idx - 1);
	}
	for(int i = 0;i < length; i++)
	{
		if (i <= r)
		{
			results[zSort[i].index] = min(max(zSort[i].value - beta, 0.0),1.0);
		}
		else
		{
			results[zSort[i].index] = min(max(zSort[i].value + beta, 0.0),1.0);
		}

	}
	delete [] zSort; delete [] zClip; delete [] zBetaRep;
	return results;
}
double Decoder::xUpdateADMMLP(double degree, double mu, double llr, double t)
{
	double x = (double)1.0/degree * t;
	if (x > 1)
		x = 1;
	else if (x < 0)
		x = 0;

	return x;
}
double Decoder::xUpdateADMML1(double degree, double mu, double llr, double t, double alpha)
{
	double temp_threshold = degree/2.0;
	double x;

	if (t < temp_threshold)
		x = (double)1.0/degree * (t - alpha/mu);
	else
		x = (double)1.0/degree * (t +  alpha/mu);

	if (x > 1)
		x = 1;
	else if (x < 0)
		x = 0;
	return x;
}
double Decoder::xUpdateADMML2(double degree, double mu, double llr, double t, double alpha)
{
	double x = (t  -  alpha/mu)/(degree - 2*alpha/mu);
	if (x > 1)
		x = 1;
	else if (x < 0)
		x = 0;
	return x;
}
double Decoder::xUpdateADMMEntropy(double prevx, double degree, double mu, double llr, double t, double alpha)
{
	double x = newtonsol(0, prevx, 10, t, degree, mu, alpha);
	if (x > 1)
		x = 1;
	else if (x < 0)
		x = 0;
	return x;
}
void Decoder::ADMMDecoder()
{
	double feas_tol = para_end_feas;
	int maxIter = maxIteration;

	int numIter = 0;
	double mu = para_mu;
	double rho = para_rho;

	double *x = new double[mBlocklength];
	double *temp = new double[mBlocklength];
	double *Lambda = new double[mPCheckMapSize];
	double *zReplica = new double[mPCheckMapSize];
	double *zOld = new double[mPCheckMapSize];
	double *xPred = new double[mPCheckMapSize];


	//initialize values

	for (int i = 0;i < mPCheckMapSize; i++)
	{	
		Lambda[i] = 0;
		zReplica[i] = 0.5;	
	}

	for(int i = 0; i < maxIter; i++)
	{
		mIteration = i + 1;

		clock_t testStart = clock();
		int warmstart_activation_number = 0;
		
		//update x
		// calculate PCheckMap'*(Lambda + z)
		for (int j = 0; j < mBlocklength; j++)
			temp[j] = 0;
		
		for(int j = 0; j < mPCheckMapSize; j++)
		{
			if (mPCheckMap[j].col > mBlocklength) // warning, can be removed
				cout<<"warning: "<<mPCheckMap[j].col<<endl;
			temp[mPCheckMap[j].col] += (zReplica[mPCheckMap[j].row] + Lambda[mPCheckMap[j].row]);
		}

		// calculate new x and projection on to [0,1]
		for(int j = 0; j < mBlocklength; j++)
		{
			double llr = _LogLikelihoodRatio[j];
			double t = temp[j] - _LogLikelihoodRatio[j] / mu;
			double deg = (double)VariableDegree[j];

			switch(enumADMM)
			{
			case LinearProgram:
				x[j] = xUpdateADMMLP(deg, mu, llr, t);
				break;
			case L1Penalty:
				x[j] = xUpdateADMML1(deg, mu, llr, t, alpha);
				break;
			case L2Penalty:
				x[j] = xUpdateADMML2(deg, mu, llr, t, alpha);
				break;
			case EntropyPenalty:
				if(i == 0)
					x[j] = 0.5;
				x[j] = xUpdateADMMEntropy(x[j], deg, mu, llr, t, alpha);
				break;
			default:
				cout<<"ADMM not supported"<<endl;
			}
#ifdef DEBUG
			//cout<<x[j]<<" | "<<(double)1.0/VariableDegree[j] * (temp[j] - _LogLikelihoodRatio[j] / mu)<<endl;
#endif
		}

		//update z
		for (int j = 0; j < mPCheckMapSize; j++)
			xPred[j] = 0;

		for(int j = 0; j < mPCheckMapSize; j++)
		{
			xPred[mPCheckMap[j].row] += x[mPCheckMap[j].col];
		}
		// get old z
		for(int j = 0; j < mPCheckMapSize; j++)
		{
			zOld[j] = zReplica[j];
		}
		//do projection for each check node
		int CumSumCheckDegree = 0; // cumulative position of currect edge in factor graph
		for(int j = 0; j < mNChecks; j++)
		{
			double warmstartdiff = 0;
			NODE *vector_before_proj = new NODE[CheckDegree[j]]; // initialize projection vector
			for(int k = 0; k < CheckDegree[j]; k++) // fill out the vector
			{
				vector_before_proj[k].index = k; // set index, for reconstruction after sorting
				vector_before_proj[k].value = rho * xPred[CumSumCheckDegree + k] + (1.0 - rho) * zOld[CumSumCheckDegree + k] - Lambda[CumSumCheckDegree + k];
				if(i == 0)
				{
					warmstartchecks[CumSumCheckDegree + k] = vector_before_proj[k].value;
					warmstartdiff = 999;
				}
				else
				{
					warmstartdiff += abs(warmstartchecks[CumSumCheckDegree + k] - vector_before_proj[k].value);
					warmstartchecks[CumSumCheckDegree + k] = vector_before_proj[k].value;
				}
			}
			
			if(warmstartdiff > feas_tol)//(double)0.1/(i*i+1))
			{
				double *ztemp  = mProjectionPPd(vector_before_proj, CheckDegree[j]);
				for(int k = 0; k < CheckDegree[j]; k++)
				{
					zReplica[CumSumCheckDegree + k] = ztemp[k];
					warmstartsave[CumSumCheckDegree + k] = ztemp[k];

				}
				delete [] ztemp;
			}
			else
			{
				for(int k = 0; k < CheckDegree[j]; k++)
				{
					zReplica[CumSumCheckDegree + k] = warmstartsave[CumSumCheckDegree + k];
				}
			}

			CumSumCheckDegree += CheckDegree[j];
			
			delete [] vector_before_proj;
		}

		for(int j = 0; j < mBlocklength; j++)
		{
			OutputFromDecoder[j] = x[j];
		}

		/////////Important tweak/////////////////
		if(ValidateCodeword())
		{
			mAlgorithmConverge = true;
			mValidCodeword = true;
			break;
		}

		double feas = -9999;
		for(int j = 0; j < mPCheckMapSize; j++)
		{
			Lambda[j] = Lambda[j] + (rho * (zReplica[j] - xPred[j]) + (1 - rho) * (zReplica[j] - zOld[j]));
			double absolute = zReplica[j] > xPred[j] ? zReplica[j] - xPred[j] : xPred[j] - zReplica[j];
			if(absolute > feas)
			{
				feas = absolute;
			}
		}

		if (feas < feas_tol)
		{
			mAlgorithmConverge = true;
			break;
		}
		
	}
	for(int i = 0; i < mBlocklength; i++)
	{
		OutputFromDecoder[i] = x[i];
	}
#ifdef DEBUG
	if (numIter > 10000)
		cout<<"numIter = "<<numIter<<endl;
#endif
	delete [] x;
	delete [] Lambda;
	delete [] zReplica;
	delete [] zOld;
	delete [] xPred;
	delete [] temp;

}
double Decoder::NonSatUpd(double llr1, double llr2)	
{
	double pos_llr1, pos_llr2;
	int sign_llr1, sign_llr2;
	if(llr1 > 0)
	{
		pos_llr1 = llr1;
		sign_llr1 = 1;
	}
	else
	{
		pos_llr1 = -llr1;
		sign_llr1 = -1;
	}
	if(llr2 > 0)
	{
		pos_llr2 = llr2;
		sign_llr2 = 1;
	}
	else
	{
		pos_llr2 = -llr2;
		sign_llr2 = -1;
	}
	double llrout = sign_llr1 * sign_llr2 * min(pos_llr1, pos_llr2);
	double sumsum = llr1 + llr2;
	double minusminus = llr1 - llr2;
	llrout = llrout + log(1 + exp(- abs(sumsum))) - log(1 + exp(- abs(minusminus)));
	return llrout;
}
void Decoder::BPDecoder()
{
	bool error = false;
	//double feas_tol = para_end_feas;
	int maxIter = maxIteration;
	int numIter = 0;
	for (int i = 0; i < mPCheckMapSize; i++)
	{
		v[i].Mprob =  _LogLikelihoodRatio[v[i].Mcol];
	}
	for(int i = 0; i < mBlocklength; i++)
	{
		u0[i] = _LogLikelihoodRatio[i];
	}

	int itnum;
	for(itnum = 0; itnum < maxIter; itnum++)
	{
		mIteration = itnum + 1;
		//	cout<<"finished calculating u"<<endl;
		for (int i = 0; i < mPCheckMapSize; i++)
		{
			double temp;
			if(BP_non_sat)
			{
				temp = 1.0;
				temp = v[u[i].i[0]].Mprob;
				for(int j = 1; j < u[i].degree - 1; j++)
				{
					temp = NonSatUpd(temp, v[u[i].i[j]].Mprob);
				}
				if(u[i].degree % 2 != 0)
					temp = temp;
				u[i].Mprob = temp;
			}
			else
			{
				temp = 1.0;
				for(int j = 0; j < u[i].degree - 1; j++)
				{
					if(abs(v[u[i].i[j]].Mprob) < BP_sat_value)
						temp = temp * tanh(0.5 * v[u[i].i[j]].Mprob);
					else
					{
						if(v[u[i].i[j]].Mprob > 0)
							temp = temp * tanh(0.5 * BP_sat_value);
						else
							temp = temp * tanh(- 0.5 * BP_sat_value);
					}
				}
				u[i].Mprob = log((1 + temp)/(1 - temp));
				
			}
			if(my_isinf(u[i].Mprob) || my_isnan(u[i].Mprob))
			{
				cout<<"is inf or is nan"<<endl;
			}
		}
		//cout<<"one iter"<<endl;
		for (int i=0;i<mPCheckMapSize;i++)
		{
				v[i].Mprob = u0[v[i].Mcol];
				for(int j = 0; j < v[i].degree - 1; j++)
					v[i].Mprob += u[v[i].i[j]].Mprob;
					
		}

		// decision at each iteration
		for(int i = 0; i < mBlocklength; i++)
			OutputFromDecoder[i] = u0[i];
		for(int i = 0; i < mPCheckMapSize; i++)
		{
			OutputFromDecoder[u[i].Mcol] += u[i].Mprob;
		}
		for(int i = 0; i < mBlocklength; i++)
		{
			if(my_isinf(OutputFromDecoder[i]))
			{
				cout<<"isinf"<<endl;
				break;
			}
			if(my_isnan(OutputFromDecoder[i]))
			{
				cout<<"isnan"<<endl;
				break;
			}
		}
		if (ValidateCodeword())
		{
			mValidCodeword = true;
			mAlgorithmConverge = true;
			break;
		}
	}
	
}
void Decoder::MSADecoder()
{
	bool error = false;
	bool MSAConverge = false;
	//double feas_tol = para_end_feas;
	int maxIter = maxIteration;
	int numIter = 0;
	for(int i = 0; i < mPCheckMapSize; i++)
	{
		v[i].Mprob = - _LogLikelihoodRatio[v[i].Mcol];
	}
	for(int i = 0; i < mBlocklength; i++)
	{
		u0[i] = -_LogLikelihoodRatio[i];
	}


	int itnum;
	for(itnum = 0; itnum < maxIter; itnum++)
	{
		for (int i = 0; i < mPCheckMapSize; i++)
		{
			double temp = abs(v[u[i].i[0]].Mprob);
			int sign = 1;
			for(int j = 0; j < u[i].degree - 1; j++)
			{
				double magnitude;
				if(v[u[i].i[j]].Mprob < 0)
				{
					sign = (-1) * sign;
					magnitude = -v[u[i].i[j]].Mprob;
				}
				else
				{
					magnitude = v[u[i].i[j]].Mprob;
				}
				if(magnitude < temp)
					temp = magnitude;
			}
			temp = sign * temp;

			if(u[i].degree % 2 != 0)
				temp = -temp;
			u[i].Mprob = temp;

		}

		for (int i = 0; i < mPCheckMapSize; i++)
		{
				v[i].Mprob = u0[v[i].Mcol];
				for(int j = 0; j < v[i].degree - 1; j++)
					v[i].Mprob += u[v[i].i[j]].Mprob;		
		}

		// decision at each iteration
		for(int i = 0; i < mBlocklength; i++)
		{
			OutputFromDecoder[i] = u0[i]; 
		}
		for(int i = 0; i < mPCheckMapSize; i++)
		{
			OutputFromDecoder[u[i].Mcol] += u[i].Mprob;
		}
		int errorbits = 0;
		for(int i = 0; i < mBlocklength; i++)
		{
			if(my_isinf(OutputFromDecoder[i]))
			{
				cout<<"isinf"<<endl;
				//break;
			}
			if(my_isnan(OutputFromDecoder[i]))
			{
				cout<<"isnan"<<endl;
				//break;
			}
			if(abs(OutputFromDecoder[i]) > 10000000)
			{
				cout<<"big number"<<endl;
				//break;
			}
		}
		if(itnum > 5)
		{
			if (ValidateCodeword())
			{
				mValidCodeword = true;
				mAlgorithmConverge = true;
				break;
			}
		}
	}
	mIteration = itnum + 1;
}
bool Decoder::ValidateCodeword()
{
	int *syndrome = new int[mNChecks];
	for(int i = 0; i < mNChecks; i++)
		syndrome[i] = 0;
	for(int i = 0; i < mPCheckMapSize; i++)
	{
		if(DecoderName == "ADMM")
		{
			syndrome[mParityCheckMatrix[i].row] += (int) floor(OutputFromDecoder[mParityCheckMatrix[i].col] + 0.5);
		}
		else // DecoderName == "BP" or "MSA"
		{
			if(OutputFromDecoder[mParityCheckMatrix[i].col] < 0)
				syndrome[mParityCheckMatrix[i].row]++;
		}
	}
	bool validcodeword = true;
	for(int i = 0; i < mNChecks; i++)
	{
		if (syndrome[i]%2 == 1)
		{
			validcodeword = false;
			break;
		}
	}
	delete [] syndrome;
	return validcodeword;
}
bool Decoder::ValidateCodeword(int input[])
{
	int *syndrome = new int[mNChecks];
	for(int i = 0; i < mNChecks; i++)
		syndrome[i] = 0;
	for(int i = 0; i < mPCheckMapSize; i++)
	{
		syndrome[mParityCheckMatrix[i].row] += input[mParityCheckMatrix[i].col];
	}
	bool validcodeword = true;
	for(int i = 0; i < mNChecks; i++)
	{
		if (syndrome[i]%2 == 1)
		{
			validcodeword = false;
			break;
		}
	}
	delete [] syndrome;
	return validcodeword;
}



void Simulator::SetCodeword(string codewordfile)// use user defined codeword
{
	ifstream cwfile(codewordfile.c_str());
	int temp_input = 0;
	int i = 0;
	if (cwfile.is_open())
	{
		while (cwfile >> temp_input) 
		{
			mInput[i] = temp_input;
			i++;
			if(i >= mBlocklength)
				break;
		}
	}
	if(!mDecoder.ValidateCodeword(mInput))
	{
		cout<<"invalid c.w., use all zero instead"<<endl;
		for(int i = 0; i < mBlocklength; i++)
			mInput[i] = 0;
	}
}
void Simulator::mUpdateStats()
{
	mErrorFlag = false;
	for(int i = 0; i < mBlocklength; i++)
	{
		if(mInput[i] != mDecodedSequence.HardDecision[i]){mErrorFlag = true; break;} // decoded seq does not agree with the input
		if(mDecodedSequence.IsNan){mErrorFlag = true; break;}// decoded seq is nan
		if(mDecodedSequence.IsInf){mErrorFlag = true; break;}// decoded seq is inf
	}
	if(mErrorFlag)
	{
		mTotalErrors++;
		mTotalWrongIterations += mDecodedSequence.Iteration;
		mTotalWrongDecodingTime += mDecodedSequence.ExeTime;
		if(mDecodedSequence.ValidCodeword)
		{
			mTotalUndetectedErrors++; 
			if(mChannelName == "AWGN")
			{
				double dis_orig = 0;
				double dis_decode = 0;
				for(int i = 0; i < mBlocklength; i++)
				{
					dis_orig += ((double)mInput[i]*2.0 - 1 - mNoisySequence.NoisySequence[i])*((double)mInput[i]*2.0 - 1 - mNoisySequence.NoisySequence[i]);
					dis_decode += ((double)floor(mDecodedSequence.HardDecision[i] + 0.5)*2.0 - 1 - mNoisySequence.NoisySequence[i]) * ((double)floor(mDecodedSequence.HardDecision[i] + 0.5)*2.0 - 1 - mNoisySequence.NoisySequence[i]);
				}
				if (dis_orig > dis_decode)
				{
					mTotalMLErrors++;
				}
			}
			else if(mChannelName == "BSC")
			{
				double hmdis_orig = 0;
				double hmdis_decoder = 0;
				for(int i = 0; i < mBlocklength; i++)
				{
					hmdis_orig += abs(mInput[i] -  mNoisySequence.NoisySequence[i]);
					hmdis_decoder += abs(mDecodedSequence.HardDecision[i] -  mNoisySequence.NoisySequence[i]);
				}
				if (hmdis_orig > hmdis_decoder)
				{
					mTotalMLErrors++;
				}
			}
		}
	}
	else
	{
		mTotalCorrectIterations += mDecodedSequence.Iteration;
		mTotalCorrectDecodingTime += mDecodedSequence.ExeTime;
	}
	mTotalSims++;
	mTotalDecodingTime += mDecodedSequence.ExeTime;
	mTotalIterations += mDecodedSequence.Iteration;


}
void Simulator::RunSim()
{
	clock_t tStart = clock();
	mClear(); // clear stats
	while(true)
	{
			
		mChannel.GenerateOutput(mInput, mNoisySequence);
		mDecoder.Decode(mNoisySequence, mDecodedSequence);
		mUpdateStats(); // update stats
		mTotalExeTime = clock() - tStart;
		mOutputCommand();
			
		if(mTotalErrors >= mTargetErrors && mTargetErrors > 0)
			break;
		if(mTotalSims >= mTargetSims && mTargetSims > 0)
			break;
		if(mTotalExeTime >= mTargetExeTime && mTargetExeTime > 0)
			break;
			
	}

    
    mTotalExeTime = clock() - tStart;
}