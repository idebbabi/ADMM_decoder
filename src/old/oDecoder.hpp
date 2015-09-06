/*
 * oDecoder.hpp
 *
 *  Created on: 17 avr. 2015
 *      Author: legal
 */

#ifndef _oDecoder_HPP_
#define _oDecoder_HPP_

#include "oMatriceProjection.hpp"

class oDecoder{
public:
	NODE   zSort   [32]; // DECLARATION STATIQUE POUR GAGNER DU TEMPS !
	NODE   zClip   [32]; // DECLARATION STATIQUE POUR GAGNER DU TEMPS !
	NODE   zBetaRep[32]; // DECLARATION STATIQUE POUR GAGNER DU TEMPS !
	double results [32]; // DECLARATION STATIQUE POUR GAGNER DU TEMPS !
	float* vector_before_proj; // BLG initialize projection vector

	#define NEW_PROJECTION
//	_MatriceProjection<double,  6> mp_6;
//	_MatriceProjection<double,  7> mp_7;
	oMatriceProjection mp;

	void SetDecoder(string foDecoder);
private:

	string oDecoderName;
	string ChannelName;
	string mRealoDecoderName;

	//! Projection Algorithm. Use two-slice property.
    /*!
      \param v Input vector
	  \param length Length of the vector
	  \return Projection onto the parity polytope
    */
//	double* mProjectionPPd   (NODE v[], int length);

	// parameters for oDecoders
	string mPCMatrixFileName;  /*!< parity check matrix file */
	int mNChecks, mBlocklength, mPCheckMapSize;
	int *CheckDegree;
	int *VariableDegree;
	double RateOfCode; /*!< code rate, k/n */

	TRIPLE *mPCheckMap; /*!< parity check mapping */
	TRIPLE *mParityCheckMatrix; /*!< parity check matrix */

	MESSAGE<double> *u; /*!< messages in message passing */
	MESSAGE<double> *v;
	double *u0;


	// algorithm parameter
	double para_end_feas;   /*!<  ADMM alg end tolerance, typically 1e-5 */
	double para_mu;         /*!< ADMM alg \mu, typically 5.5 */
	double para_rho;        /*!< ADMM over relaxation para, typically 1.8, 1.9 */
	int maxIteration;       /*!< max iteration */

	bool BP_non_sat; /*!< true if nonsaturating BP is used */
	double BP_sat_value;  /*!< saturation value if saturating BP is used*/

	//! Initialize BP bipartite graph. Note this function is NOT efficient!!!
	/*!
	  Link u message with v messages. Simply search for all possible edges and set pointer.
	  The BP functions provided are not fully optimized. Please use with care.
	*/
	void mInitializeMessages();

	//! Learn the degree of the parity check matrix
	/*!
	  The parity check matrix should be in correct format. This function is useful for irregular LDPC codes.
	*/
	void mLearnParityCheckMatrix();
	//! Read and store the parity check matrix
	void mReadParityCheckMatrix();
	void mSetDefaultParameters();

	//! Calculate log likelihood ratio using the output from the channel
	/*!
	  \param channeloutput Output from channel.
	  \sa NoisySeq
	*/
	void mGenerateLLR(NoisySeq &channeloutput);

	//oDecoders:
	void ADMMoDecoder();
	void BPoDecoder();
	void MSAoDecoder();
	// data
	double *_LogLikelihoodRatio; /*!< log-likelihood ratio from received vector */
	double *OutputFromoDecoder; /*!<soft information from oDecoder. could be pseudocodeword or message from BP  */


	//
	// TO AVOID INTENPESTIVE ALLOCATIONS !
	//
	int* syndrome; // BLG

	double alpha; /*!< the constant for penalty term */

	double normalizeBSC;

	double *warmstartchecks;	/*!< warm start ADMM */
	double *warmstartsave;	/*!< save results for warm start*/

	//oDecoder stats
	unsigned long int mExeTime; /*!< exevution time */
	bool mAlgorithmConverge; /*!< true if the oDecoder converges*/
	int mIteration; /*!< number of iterations used for decoding */
	bool mValidCodeword; /*!< true if the output is a valid codeword */
	// ADMM oDecoder select
	enum ADMMdecoding {
                 LinearProgram, /*!< Enum value LP. */
                 L1Penalty, /*!< Enum value L1. */
                 L2Penalty,  /*!< Enum value L2. */
				 EntropyPenalty /*!< Enum value Entropy. */
               }
		enumADMM;
	// function for different x-updates
	double xUpdateADMMLP(double degree, double mu, double llr, double t);
	double xUpdateADMML1(double degree, double mu, double llr, double t, double alpha);
	double xUpdateADMML2(double degree, double mu, double llr, double t, double alpha);
	double xUpdateADMMEntropy(double prevx, double degree, double mu, double llr, double t, double alpha);
	// function for nonsaturating BP
	double NonSatUpd(double llr1, double llr2);


	// functions for entropy penalty
	double entropy_th; /*!< Threshold for entropy penalty */
	double funcWithEntropy(double x, double t, double d, double mu, double alpha)
	{
		double upth = 1 - entropy_th, lowth = entropy_th;
        double y = 0.0;
		if(x < upth && x > lowth)
		{
			y = t/d + alpha/d/mu*log(x) - alpha/d/mu*log(1-x) - x;
		}
		if(x >= upth)
		{
			y = t/d + alpha/d/mu*log(upth) - alpha/d/mu*log(1 - upth) - upth
				+ (alpha/d/mu/upth/(1 - upth) - 1)*(x - upth);
		}
		if(x <= lowth)
		{
			y = t/d + alpha/d/mu*log(lowth) - alpha/d/mu*log(1 - lowth) - lowth
				+ (alpha/d/mu/lowth/(1 - lowth) - 1)*(x - lowth);
		}
		return y;
	}

    double funcWithEntropyDrv(double x, double t, double d, double mu, double alpha)
	{
		double upth = 1 - entropy_th, lowth = entropy_th;
		double y = 0.0;
		if(x < upth && x > lowth)
		{
			y = alpha/d/mu/x/(1 - x) - 1;
		}
		if(x >= upth)
		{
			y = alpha/d/mu/upth/(1 - upth) - 1;
		}
		if(x <= lowth)
		{
			y = alpha/d/mu/lowth/(1 - lowth) - 1;
		}
		return y;
	}

    double newtonsol(double y, double xinit, int iter, double t, double d, double mu, double alpha)
	{
		if(iter < 0)
			iter = 5;
		double x = xinit;
		for(int i = 0; i < iter; i++)
		{
			double a = funcWithEntropy(x,t,d,mu,alpha);
			double b = funcWithEntropyDrv(x,t,d,mu,alpha);
			double xnew = x - (a - y)/b;
			if(xnew > 1) xnew = 1;
			if(xnew < 0) xnew = 0;
			if(abs (x - xnew) < 0.0000001)
			{
				x = xnew;
				break;
			}
			x = xnew;
		}
		return x;
	}

public:
	//! Set oDecoder parameters.
	/*!
	 \param mxIt Maximum number of iterations. Default: 1000
	 \param feas ADMM stopping parameter. Default: 1e-6
	 \param p_mu ADMM step parameter mu. Default: 5
	 \param p_rho ADMM over relaxation parameter. Default: 1.8
	*/
	void SetParameters(int mxIt, double feas,  double p_mu, double p_rho);
	//! Set BSC normalized log likelihood ratios.
	/*!
	 \param normal if true, then the LLRs are 1 or -1. if false, then it is log(p/(1-p)) or log((1-p)/p)
	*/
	void SetBSCnormalized(bool normal);
	//! Set non-saturating BP oDecoder.
	/*!
	 \param nonsat True if the oDecoder should be NON-Saturated oDecoder
	 \param value Saturation level for saturated oDecoder. nonsat should be false for this option to work.
	*/
	void SetBPSaturation(bool nonsat, double value);
	//! Set oDecoder type.
	/*!
	 \param oDecoder oDecoder name. Currently support "BP", "MSA", "ADMMLP", "ADMML1", "ADMML2" and "ADMMEntropy"
	*/
	void SetoDecoder(string oDecoder);
	//! Set channel type.
	/*!
	 \param channel Channel name. Currently support "BSC" and "AWGN"
	*/
	void SetChannel(string channel){ChannelName = channel;}
	//! Set entropy penalty constant
	/*!
	 \param alp Penalty constant
	*/
	void SetPenaltyConstant(double alp){alpha = alp;}
	//! Set entropy penalty threshold
	/*!
	 \param th Threshold
	*/
	void SetEntropyParameter(double th){entropy_th = th;}
	//! See if the output from oDecoder is valid codeword
	/*!
		Note that the output soft information from ADMM and message passing algorithms(BP and MSA) are different.
	*/
	bool ValidateCodeword();
	//! See if the input sequence is valid codeword
	/*!
		\param input Input sequence array.
	*/
	bool ValidateCodeword(int input[]);
	//! Invoke decode function
	/*!
		\param channeloutput Noisy sequence from channel.
		\param oDecoderesult Decode results. Contains decoded sequence and other information.
		\sa DecodedSeq
		\sa NoisySeq
	*/
	double Decode(NoisySeq &channeloutput, DecodedSeq &oDecoderesult);

	//! Constructor.
	/*!
		\param FileName File for parity check matrix.
		\param blocklength Blocklength of code (i.e. n)
		\param nChecks Number of checks (i.e. row number of parity check matrix or (n-k))
	*/
	oDecoder(string FileName,int nChecks, int BlockLength);
	//! Destructor
	~oDecoder();
};


// oDecoder Class
oDecoder::oDecoder(string FileName,int nChecks, int BlockLength) : mp()
{
	// Learn and read parity check matrix
	mPCMatrixFileName = FileName;
	mNChecks = nChecks;
	mBlocklength = BlockLength;
	mPCheckMapSize = 0;
	CheckDegree = new int[mNChecks];
	VariableDegree = new int[mBlocklength];
	OutputFromoDecoder = new double[mBlocklength];
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

	syndrome           = new int[mNChecks]; // BLG
	vector_before_proj = new float[ 64 ];    // BLG initialize projection vector

}
oDecoder::~oDecoder()
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
	delete [] OutputFromoDecoder;
	delete [] _LogLikelihoodRatio;
	delete [] vector_before_proj; // BLG initialize projection vector

	delete [] syndrome; // BLG

}

//oDecodersettings
void oDecoder::SetParameters(int mxIt, double feas,  double p_mu, double p_rho)
{
	maxIteration = mxIt;
	para_end_feas = feas;
	para_mu = p_mu;
	para_rho = p_rho;
}
void oDecoder::SetBSCnormalized(bool normal)
{
	if(normal)
		normalizeBSC = 1;
	else
		normalizeBSC = 0;
}
void oDecoder::SetBPSaturation(bool nonsat, double value)
{
	BP_non_sat = nonsat;
	BP_sat_value = value;
}
void oDecoder::SetDecoder(string foDecoder)
{
	mRealoDecoderName = foDecoder;
	if(foDecoder == "BP")
	{
		oDecoderName = foDecoder;
		mInitializeMessages();
	}
	else if(foDecoder == "ADMMLP")
	{
		oDecoderName = "ADMM";
		enumADMM = LinearProgram;
	}
	else if(foDecoder == "ADMML1")
	{
		oDecoderName = "ADMM";
		enumADMM = L1Penalty;
	}
	else if(foDecoder == "ADMML2")
	{
		oDecoderName = "ADMM";
		enumADMM = L2Penalty;
	}
	else if(foDecoder == "ADMMEntropy")
	{
		oDecoderName = "ADMM";
		enumADMM = EntropyPenalty;
	}
	else if(foDecoder == "MSA")
	{
		oDecoderName = "MSA";
        mInitializeMessages();
	}
	else
	{
		cout<<"(oDecoder) choice is not supported : " << foDecoder <<endl;
		exit(0);
	}
}

// oDecoder initializations
void oDecoder::mSetDefaultParameters()
{
	// oDecoder settings
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
void oDecoder::mLearnParityCheckMatrix()// learn the degree of the parity check matrix
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

void oDecoder::mReadParityCheckMatrix()// read the parity check matrix, initialize message passing
{
	ifstream myfile (mPCMatrixFileName.c_str());
	int tempvalue = 0;
	int i = 0;
	int count = 0;
	string line;

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
				u[count].Mcol  = mParityCheckMatrix[count].col;
				u[count].Mrow  = mParityCheckMatrix[count].row;
				u[count].Mprob = 0;

				v[count].SetDeg(VariableDegree[tempvalue]);
				v[count].Mcol  = mParityCheckMatrix[count].col;
				v[count].Mrow  = mParityCheckMatrix[count].row;
				v[count].Mprob = 0;

				count++;
			}
			i++;
		}
	}
	if (myfile.is_open())
		myfile.close();
}

void oDecoder::mInitializeMessages()
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
void oDecoder::mGenerateLLR(NoisySeq &channeloutput)
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
double oDecoder::Decode(NoisySeq &channeloutput, DecodedSeq &oDecoderesult)
{
	mExeTime           = 0;
	mAlgorithmConverge = false;
	mValidCodeword     = false;

	mGenerateLLR(channeloutput);

    auto start     = chrono::steady_clock::now();
	clock_t tStart = clock();
	if(oDecoderName == "ADMM")
		ADMMoDecoder();
	else if(oDecoderName == "BP")
	{
# ifdef DEBUG
		cout<<"BP"<<endl;
#endif
		BPoDecoder();
	}
	else if(oDecoderName == "MSA")
	{
# ifdef DEBUG
		cout<<"MSA"<<endl;
#endif
		MSAoDecoder();
	}
	else
	{
		cout<<"oDecoder not supported"<<endl;
		exit(0);
	}
    auto end       = chrono::steady_clock::now();
	mExeTime                       = clock() - tStart;
	oDecoderesult.ExeTime           = mExeTime;
	oDecoderesult.Iteration         = mIteration;
	oDecoderesult.ValidCodeword     = mValidCodeword;
	oDecoderesult.AlgorithmConverge = mAlgorithmConverge;
	oDecoderesult.Iteration         = mIteration;

	for(int i = 0; i < mBlocklength; i++)
	{
		if(my_isnan(OutputFromoDecoder[i]))
			oDecoderesult.IsNan = true;
		if(my_isinf(OutputFromoDecoder[i]))
			oDecoderesult.IsInf = true;
		oDecoderesult.SoftInfo[i] = OutputFromoDecoder[i];

		if(oDecoderName == "ADMM")
		{
			oDecoderesult.HardDecision[i] = floor(OutputFromoDecoder[i] + 0.5);
		}
		else // oDecoderName == "BP" or "MSA"
		{
			if(OutputFromoDecoder[i] < 0)
				oDecoderesult.HardDecision[i] = 1;
			else
				oDecoderesult.HardDecision[i] = 0;
		}
	}
    auto diff      = end - start;
    return chrono::duration <double, milli> (diff).count();
}

double oDecoder::xUpdateADMMLP(double degree, double mu, double llr, double t)
{
	double x = (double)1.0/degree * t;
	if (x > 1)
		x = 1;
	else if (x < 0)
		x = 0;

	return x;
}

double oDecoder::xUpdateADMML1(double degree, double mu, double llr, double t, double alpha)
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
double oDecoder::xUpdateADMML2(double degree, double mu, double llr, double t, double alpha)
{
	double x = (t  -  alpha/mu)/(degree - 2*alpha/mu);
	if (x > 1)
		x = 1;
	else if (x < 0)
		x = 0;
	return x;
}
double oDecoder::xUpdateADMMEntropy(double prevx, double degree, double mu, double llr, double t, double alpha)
{
	double x = newtonsol(0, prevx, 10, t, degree, mu, alpha);
	if (x > 1)
		x = 1;
	else if (x < 0)
		x = 0;
	return x;
}


void oDecoder::ADMMoDecoder()
{
	double feas_tol = para_end_feas;
	int maxIter     = maxIteration;

	int numIter = 0;
	double mu   = para_mu;
	double rho  = para_rho;

	double *x        = new double[mBlocklength];
	double *temp     = new double[mBlocklength];
	double *Lambda   = new double[mPCheckMapSize];
	double *zReplica = new double[mPCheckMapSize];
	double *zOld     = new double[mPCheckMapSize];
	double *xPred    = new double[mPCheckMapSize];

	long long int min = 1000000;
	long long int max = 0;
	long long int sum = 0;
	long long int exe = 0;

	//initialize values

	for (int i = 0;i < mPCheckMapSize; i++)
	{
		Lambda[i]   = 0;
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
			double t   = temp[j] - _LogLikelihoodRatio[j] / mu;
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

		//
		// get old z
		//
		for(int j = 0; j < mPCheckMapSize; j++)
		{
			zOld[j] = zReplica[j];
		}

		//
		//do projection for each check node
		//
		int CumSumCheckDegree = 0; // cumulative position of currect edge in factor graph
		for(int j = 0; j < mNChecks; j++)
		{
			double warmstartdiff = 0;
			NODE vector_before_proj[32]; // initialize projection vector
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
				double *ztemp;

				const auto start  = timer();
				ztemp = mp.mProjectionPPd(vector_before_proj, CheckDegree[j]);
				const auto temps  = timer() - start;
				min  = (min < temps) ? min : temps;
				max  = (max > temps) ? max : temps;
				sum += temps;
				exe += 1;

				for(int k = 0; k < CheckDegree[j]; k++)
				{
					zReplica     [CumSumCheckDegree + k] = ztemp[k];
					warmstartsave[CumSumCheckDegree + k] = ztemp[k];

				}
				//delete [] ztemp;
			}
			else
			{
				for(int k = 0; k < CheckDegree[j]; k++)
				{
					zReplica[CumSumCheckDegree + k] = warmstartsave[CumSumCheckDegree + k];
				}
			}

			CumSumCheckDegree += CheckDegree[j];

// BLG			delete [] vector_before_proj;
		}

		//memcpy( );

		for(int j = 0; j < mBlocklength; j++)
		{
			OutputFromoDecoder[j] = x[j];
		}

		/////////Important tweak/////////////////
		if( ValidateCodeword() )
		{
			mAlgorithmConverge = true;
			mValidCodeword     = true;
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
		OutputFromoDecoder[i] = x[i];
	}
#ifdef DEBUG
	if (numIter > 10000)
		cout<<"numIter = "<<numIter<<endl;
#endif
	const float avg = sum / exe;
	const auto start  = timer();
	const float fake  = tanh(avg);
	const auto t_tanh = timer() - start;
	printf("mProjectionPPd : min = %ld, max = %ld, avg = %ld [tanh = %d - fake: %1.3f]\n", min, max, (int)avg, t_tanh, fake);

	delete [] x;
	delete [] Lambda;
	delete [] zReplica;
	delete [] zOld;
	delete [] xPred;
	delete [] temp;

}

double oDecoder::NonSatUpd(double llr1, double llr2)
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


void oDecoder::BPoDecoder()
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
				temp = v[ u[i].i[0] ].Mprob;
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
			OutputFromoDecoder[i] = u0[i];
		for(int i = 0; i < mPCheckMapSize; i++)
		{
			OutputFromoDecoder[u[i].Mcol] += u[i].Mprob;
		}
		for(int i = 0; i < mBlocklength; i++)
		{
			if(my_isinf(OutputFromoDecoder[i]))
			{
				cout<<"isinf"<<endl;
				break;
			}
			if(my_isnan(OutputFromoDecoder[i]))
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


void oDecoder::MSAoDecoder()
{
	bool error       = false;
	bool MSAConverge = false;
	int maxIter      = maxIteration;
	int numIter      = 0;

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
		//
		//
		//
		for (int i = 0; i < mPCheckMapSize; i++)
		{
			int sign    = 1;
			double temp = abs(v[u[i].i[0]].Mprob);
			for(int j = 0; j < u[i].degree - 1; j++)
			{
				double magnitude;
				if( v[ u[i].i[j] ].Mprob < 0 )
				{
					sign      = (-1) * sign;
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

		//
		//
		//
		for (int i = 0; i < mPCheckMapSize; i++)
		{
				v[i].Mprob = u0[v[i].Mcol];
				for(int j = 0; j < v[i].degree - 1; j++)
					v[i].Mprob += u[v[i].i[j]].Mprob;
		}

		//
		// decision at each iteration
		//
		for(int i = 0; i < mBlocklength; i++)
		{
			OutputFromoDecoder[i] = u0[i];
		}

		//
		//
		//
		for(int i = 0; i < mPCheckMapSize; i++)
		{
			OutputFromoDecoder[u[i].Mcol] += u[i].Mprob;
		}

		//
		// DES SANITY CHECKS
		//
		//int errorbits = 0;
		for(int i = 0; i < mBlocklength; i++)
		{
			if(my_isinf(OutputFromoDecoder[i]))
			{
				cout<<"isinf"<<endl;
				//break;
			}
			if(my_isnan(OutputFromoDecoder[i]))
			{
				cout<<"isnan"<<endl;
				//break;
			}
			if(abs(OutputFromoDecoder[i]) > 10000000)
			{
				cout << "big number OutputFromoDecoder[" << i << "] = " << OutputFromoDecoder[i] <<endl;
				//break;
			}
		}


		//
		//
		//
		if(itnum > 5)
		{
			if (ValidateCodeword())
			{
				mValidCodeword     = true;
				mAlgorithmConverge = true;
				break;
			}
		}
	}

	mIteration = itnum + 1;
}


bool oDecoder::ValidateCodeword()
{
	// BLG	int *syndrome = new int[mNChecks];
	for(int i = 0; i < mNChecks; i++){
		syndrome[i] = 0;
	}

	if(oDecoderName == "ADMM")
	{
		for(int i = 0; i < mPCheckMapSize; i++)
		{
			syndrome[mParityCheckMatrix[i].row] += (int) floor(OutputFromoDecoder[mParityCheckMatrix[i].col] + 0.5);
		}
	}
	else // oDecoderName == "BP" or "MSA"
	{
		for(int i = 0; i < mPCheckMapSize; i++)
		{
			if(OutputFromoDecoder[mParityCheckMatrix[i].col] < 0)
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
	// BLG 	delete [] syndrome;
	return validcodeword;
}


bool oDecoder::ValidateCodeword(int input[])
{
// BLG	int *syndrome = new int[mNChecks];
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
// BLG	delete [] syndrome;
	return validcodeword;
}

#endif /* oDecoder_HPP_ */
