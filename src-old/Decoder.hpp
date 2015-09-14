/*
 * Decoder.hpp
 *
 *  Created on: 17 avr. 2015
 *      Author: legal
 */
//#define PROFILE_ON

#ifndef DECODER_HPP_
#define DECODER_HPP_

#include "MatriceProjection.hpp"
#include "./intel/RDTSC.hpp"
#define type float // BLG

class Decoder{
public:
	NODE   zSort   [32]; // DECLARATION STATIQUE POUR GAGNER DU TEMPS !
	NODE   zClip   [32]; // DECLARATION STATIQUE POUR GAGNER DU TEMPS !
	NODE   zBetaRep[32]; // DECLARATION STATIQUE POUR GAGNER DU TEMPS !
	double results [32]; // DECLARATION STATIQUE POUR GAGNER DU TEMPS !
	type* vector_before_proj; // BLG initialize projection vector


    long long int t_vn = 0;
    long long int t_cn = 0;
    long long int t_pj = 0;
    long long int t_ex = 0;


	#define NEW_PROJECTION
//	MatriceProjection<double,  6> mp_6;
//	MatriceProjection<double,  7> mp_7;
	MatriceProjection<type, 32> mp;
private:

	string DecoderName;
	string ChannelName;
	string mRealDecoderName;

	//! Projection Algorithm. Use two-slice property.
    /*!
      \param v Input vector
	  \param length Length of the vector
	  \return Projection onto the parity polytope
    */
	double* mProjectionPPd   (NODE v[], int length);
	double* mProjectionPPd_d6(NODE v[]);
	double* mProjectionPPd_d7(NODE v[]);

	// parameters for decoders
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

	//Decoders:
	void ADMMDecoder();
	void ADMMDecoder_float();
	void BPDecoder();
	void MSADecoder();
	// data
	double *_LogLikelihoodRatio; /*!< log-likelihood ratio from received vector */
	double *OutputFromDecoder; /*!<soft information from decoder. could be pseudocodeword or message from BP  */


	//
	// TO AVOID INTENPESTIVE ALLOCATIONS !
	//
	int* syndrome; // BLG

	double alpha; /*!< the constant for penalty term */

	double normalizeBSC;

	double *warmstartchecks;	/*!< warm start ADMM */
	double *warmstartsave;	/*!< save results for warm start*/

	//decoder stats
	unsigned long int mExeTime; /*!< exevution time */
	bool mAlgorithmConverge; /*!< true if the decoder converges*/
	int mIteration; /*!< number of iterations used for decoding */
	bool mValidCodeword; /*!< true if the output is a valid codeword */
	// ADMM decoder select
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
	//! Set decoder parameters.
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
	//! Set non-saturating BP decoder.
	/*!
	 \param nonsat True if the decoder should be NON-Saturated decoder
	 \param value Saturation level for saturated decoder. nonsat should be false for this option to work.
	*/
	void SetBPSaturation(bool nonsat, double value);
	//! Set decoder type.
	/*!
	 \param decoder Decoder name. Currently support "BP", "MSA", "ADMMLP", "ADMML1", "ADMML2" and "ADMMEntropy"
	*/
	void SetDecoder(string decoder);
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
	//! See if the output from decoder is valid codeword
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
		\param decoderesult Decode results. Contains decoded sequence and other information.
		\sa DecodedSeq
		\sa NoisySeq
	*/
	void Decode(NoisySeq &channeloutput, DecodedSeq &decoderesult);



	//! Constructor.
	/*!
		\param FileName File for parity check matrix.
		\param blocklength Blocklength of code (i.e. n)
		\param nChecks Number of checks (i.e. row number of parity check matrix or (n-k))
	*/
	Decoder(string FileName,int nChecks, int BlockLength);
	//! Destructor
	~Decoder();
};


// Decoder Class
Decoder::Decoder(string FileName,int nChecks, int BlockLength) : mp()
{
#ifdef PROFILE_ON
    t_vn = 0;
    t_cn = 0;
    t_pj = 0;
    t_ex = 0;
#endif

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

	syndrome           = new int[mNChecks]; // BLG
	vector_before_proj = new float[ 64 ];    // BLG initialize projection vector

}
Decoder::~Decoder()
{

#ifdef PROFILE_ON
	//File *f=fopen("time.txt","a");//imen

    if( t_ex != 0 ){
        t_cn /= t_ex;
        t_vn /= t_ex;
        t_pj /= t_ex;
        long long int sum = t_cn + t_vn;
        double CN = (100.0 * t_cn) / (double)sum;
        double VN = (100.0 * t_vn) / (double)sum;
        double PJ = (100.0 * t_pj) / (double)t_cn;
        printf("VN : %1.3f - CN = %1.3f (PJ = %1.3f)\n", VN, CN, PJ);
printf("cycles VN : %ld - cycles CN = %ld \n", t_vn, t_cn);
printf("cy_per_node VN: %1.3f - cy_per_node CN: %1.3f - cy_per_node PJ = %1.3f \n", (double) t_vn/mBlocklength, (double) t_cn/mNChecks, (double) t_pj/mNChecks);

         // fprintf("%1.3f %1.3f %1.3f\n", VN, CN, PJ);//imen
       // fclose (f);//imen
    }else{
    	printf("(WW) PROFILE WAS SWITCHED ON AND NO DATA ARE AVAILABLE !\n");
    }
#endif

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
	delete [] vector_before_proj; // BLG initialize projection vector

	delete [] syndrome; // BLG

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
	mExeTime           = 0;
	mAlgorithmConverge = false;
	mValidCodeword     = false;

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
	mExeTime                       = clock() - tStart;
	decoderesult.ExeTime           = mExeTime;
	decoderesult.Iteration         = mIteration;
	decoderesult.ValidCodeword     = mValidCodeword;
	decoderesult.AlgorithmConverge = mAlgorithmConverge;
	decoderesult.Iteration         = mIteration;

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
//===========================================================================================
double Decoder::xUpdateADMML2(double degree, double mu, double llr, double t, double alpha)
{
    float tableau[7]    = { 0.0f };
    double x ; 
    if ((mBlocklength == 576) && (mNChecks == 288))
    {
        tableau[2]  = 0.00001f;
	tableau[3]  = 2.00928f;
	tableau[6]  = 4.69438f;
	x = (double)(t  -  tableau[(int)degree]/mu)/(degree - 2*tableau[(int)degree]/mu);

    }
    else if((mBlocklength == 2304) && (mNChecks == 1152) )
    {
        tableau[2]  = 0.29669288f; 
	tableau[3]  = 0.46964023f;
	tableau[6]  = 3.19548154f;
	x = (double)(t  -  tableau[(int)degree]/mu)/(degree - 2*tableau[(int)degree]/mu);
    }
    else
    {
	x = (t  - alpha/mu)/(degree - 2*alpha/mu);
    }

	if (x > 1)
		x = 1;
	else if (x < 0)
		x = 0;
	return x;
   
}
//===========================================================================================
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
	int maxIter     = maxIteration;

	int numIter = 0;
	double mu   = para_mu;
	double rho  = para_rho;
    if ((mBlocklength == 576) && (mNChecks == 288))
    {
     	mu        = 3.37309f;//penalty
    }
    else if((mBlocklength == 2304) && (mNChecks == 1152) )
    {
    	mu        = 3.81398683f;//penalty
    }
    else{
        	mu   = para_mu;
	}
	double *x        = new double[mBlocklength];
	double *temp     = new double[mBlocklength];
	double *Lambda   = new double[mPCheckMapSize];
	double *zReplica = new double[mPCheckMapSize];
	double *zOld     = new double[mPCheckMapSize];
	double *xPred    = new double[mPCheckMapSize];

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
		//    	//
	//    	
    	// MEASURE OF THE VN EXECUTION TIME
    	//

		#ifdef PROFILE_ON
				const auto start_vn = timer();
		#endif


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
   	//
    	// MEASURE OF THE VN EXECUTION TIME
    	//
		#ifdef PROFILE_ON
				t_vn   += (timer() - start_vn);
		#endif

   	//
    	// MEASURE OF THE CN EXECUTION TIME
    	//
		#ifdef PROFILE_ON
				const auto start_cn = timer();
		#endif



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
//			NODE *vector_before_proj = new NODE[CheckDegree[j]]; // initialize projection vector
			for(int k = 0; k < CheckDegree[j]; k++) // fill out the vector
			{
				//vector_before_proj[k].index = k; // set index, for reconstruction after sorting
				vector_before_proj[k] = rho * xPred[CumSumCheckDegree + k] + (1.0 - rho) * zOld[CumSumCheckDegree + k] - Lambda[CumSumCheckDegree + k];
				if(i == 0)
				{
					warmstartchecks[CumSumCheckDegree + k] = vector_before_proj[k];
					warmstartdiff = 999;
				}
				else
				{
					warmstartdiff += abs(warmstartchecks[CumSumCheckDegree + k] - vector_before_proj[k]);
					warmstartchecks[CumSumCheckDegree + k] = vector_before_proj[k];
				}
			}

			if(warmstartdiff > feas_tol)//(double)0.1/(i*i+1))
			{
				type *ztemp;
	    	//
	    	// MEASURE OF THE PROJECTION EXECUTION TIME
	    	//
            #ifdef PROFILE_ON
    			const auto start_pj = timer();
	    #endif
				
#ifdef NEW_PROJECTION
				ztemp = mp.mProjectionPPd(vector_before_proj, CheckDegree[j]);
#else
				if( CheckDegree[j] == 6){
					ztemp = mProjectionPPd_d6(vector_before_proj);
				}else if( CheckDegree[j] == 7){
					ztemp = mProjectionPPd_d7(vector_before_proj);
				}else{
					ztemp = mProjectionPPd(vector_before_proj, CheckDegree[j]);
				}
#endif
   	//
    	// MEASURE OF THE CN EXECUTION TIME
    	//
		#ifdef PROFILE_ON
				t_pj   += (timer() - start_pj);
		#endif

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
			OutputFromDecoder[j] = x[j];
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
   	//
    	// MEASURE OF THE CN EXECUTION TIME
    	//
		#ifdef PROFILE_ON
				t_cn   += (timer() - start_cn);
		#endif

		#ifdef PROFILE_ON
		t_ex += 1;
	        #endif

		if (feas < feas_tol)
		{
			mAlgorithmConverge = true;
			break;
		}

	}
		#ifdef PROFILE_ON
		//t_ex += 1;
	        #endif
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


void Decoder::ADMMDecoder_float()
{
	float feas_tol = para_end_feas;
	int maxIter     = maxIteration;

	int numIter = 0;
	float mu    = para_mu;
	float rho   = para_rho;

	float *x        = new float[mBlocklength];
	float *temp     = new float[mBlocklength];
	float *Lambda   = new float[mPCheckMapSize];
	float *zReplica = new float[mPCheckMapSize];
	float *zOld     = new float[mPCheckMapSize];
	float *xPred    = new float[mPCheckMapSize];

	//initialize values

	for (int i = 0;i < mPCheckMapSize; i++)
	{
		Lambda[i]   = 0;
		zReplica[i] = 0.5;
	}

//	NODE *vector_before_proj = new NODE[ 64 ]; // BLG initialize projection vector

	for(int i = 0; i < maxIter; i++)
	{
		mIteration = i + 1;

//		clock_t testStart = clock();
//		int warmstart_activation_number = 0;

		//update x
		// calculate PCheckMap'*(Lambda + z)
		for (int j = 0; j < mBlocklength; j++)
			temp[j] = 0;

		for(int j = 0; j < mPCheckMapSize; j++)
		{
			const auto pos = mPCheckMap[j];
			temp[pos.col] += (zReplica[pos.row] + Lambda[pos.row]);
//			temp[mPCheckMap[j].col] += (zReplica[mPCheckMap[j].row] + Lambda[mPCheckMap[j].row]);
		}

		// calculate new x and projection on to [0,1]
		for(int j = 0; j < mBlocklength; j++)
		{
			float llr = _LogLikelihoodRatio[j];
			float t   = temp[j] - _LogLikelihoodRatio[j] / mu;
			float deg = (double)VariableDegree[j];

			switch(enumADMM)
			{
			case LinearProgram:
				x[j] = xUpdateADMMLP(deg, mu, llr, t);
				break;
			case L1Penalty:
//				float tab[6] = {0.0f, 0.0f, 0.2147f, 0.8618f, 0.0f, 0.0f, 3.68f};
//				mu = 1.0f;
//				alpha = tab[ VariableDegree[j] ];
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

		//
		//do projection for each check node
		//
		int CumSumCheckDegree = 0; // cumulative position of currect edge in factor graph
		for(int j = 0; j < mNChecks; j++)
		{
			float warmstartdiff = 0;
// INITIALIZED BEFORE			NODE *vector_before_proj = new NODE[CheckDegree[j]]; // initialize projection vector
			for(int k = 0; k < CheckDegree[j]; k++) // fill out the vector
			{
				//vector_before_proj[k].index = k; // set index, for reconstruction after sorting
				vector_before_proj[k] = rho * xPred[CumSumCheckDegree + k] + (1.0 - rho) * zOld[CumSumCheckDegree + k] - Lambda[CumSumCheckDegree + k];
				if(i == 0)
				{
					warmstartchecks[CumSumCheckDegree + k] = vector_before_proj[k];
					warmstartdiff = 999;
				}
				else
				{
					warmstartdiff += abs(warmstartchecks[CumSumCheckDegree + k] - vector_before_proj[k]);
					warmstartchecks[CumSumCheckDegree + k] = vector_before_proj[k];
				}
			}

			if(warmstartdiff > feas_tol)//(double)0.1/(i*i+1))
			{
//				double *ztemp  = mProjectionPPd(vector_before_proj, CheckDegree[j]);
				type *ztemp;
#ifdef NEW_PROJECTION
				ztemp = mp.mProjectionPPd(vector_before_proj, CheckDegree[j]);
#else
				if( CheckDegree[j] == 6){
					ztemp = mProjectionPPd_d6(vector_before_proj);
				}else if( CheckDegree[j] == 7){
					ztemp = mProjectionPPd_d7(vector_before_proj);
				}else{
					ztemp = mProjectionPPd(vector_before_proj, CheckDegree[j]);
				}
#endif
				for(int k = 0; k < CheckDegree[j]; k++)
				{
					zReplica[CumSumCheckDegree + k] = ztemp[k];
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

// BLG BLG			delete [] vector_before_proj;
		}

		for(int j = 0; j < mBlocklength; j++)
		{
			OutputFromDecoder[j] = x[j];
		}

		/////////Important tweak/////////////////
		if(ValidateCodeword())
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
			cout<<"We gone there : " << feas << " < " << feas_tol <<endl;
			mAlgorithmConverge = true;
			break;
		}else{
			cout<<"We gone here : " << feas << " < " << feas_tol <<endl;
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

	delete [] vector_before_proj; // BLG initialize projection vector

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
				if(v[u[i].i[j]].Mprob < 0)
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
			OutputFromDecoder[i] = u0[i];
		}

		//
		//
		//
		for(int i = 0; i < mPCheckMapSize; i++)
		{
			OutputFromDecoder[u[i].Mcol] += u[i].Mprob;
		}

		//
		// DES SANITY CHECKS
		//
		//int errorbits = 0;
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


bool Decoder::ValidateCodeword()
{
	// BLG	int *syndrome = new int[mNChecks];
	for(int i = 0; i < mNChecks; i++){
		syndrome[i] = 0;
	}

	if(DecoderName == "ADMM")
	{
		for(int i = 0; i < mPCheckMapSize; i++)
		{
			syndrome[mParityCheckMatrix[i].row] += (int) floor(OutputFromDecoder[mParityCheckMatrix[i].col] + 0.5);
		}
	}
	else // DecoderName == "BP" or "MSA"
	{
		for(int i = 0; i < mPCheckMapSize; i++)
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
	// BLG 	delete [] syndrome;
	return validcodeword;
}


bool Decoder::ValidateCodeword(int input[])
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

#endif /* DECODER_HPP_ */
