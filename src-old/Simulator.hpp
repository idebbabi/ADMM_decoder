/*
 * Simulator.hpp
 *
 *  Created on: 17 avr. 2015
 *      Author: legal
 */

#ifndef SIMULATOR_HPP_
#define SIMULATOR_HPP_

#include "Decoder.hpp"

//!  Class for simulation
/*!
	Run simulations.
	Take input sequence, pass it through channel to get noisy sequence. Pass noisy sequence to decoder and then compare decoded sequence with original input.
*/
class Simulator{
private:
	typedef unsigned long uint32;
	int *mInput;
	NoisySeq mNoisySequence;
	Channel mChannel;
	DecodedSeq mDecodedSequence;
	Decoder mDecoder;

	string mChannelName;
	string mDecoderName;

	string mParityCheckFile;
	string mCodewordFile;

	int mBlocklength;
	int mNChecks;
	double mRateOfCode;
	double EbN0;
	double mSTD;
	double mCrossOverProb;
	double mChannelParameter;
	bool mErrorFlag;

	// Simulation control
	uint32 mOutputFileInterval;
	string mOutputFileName;
	ofstream mOutputFileStream;

	uint32 mOutputCommandInterval;
	uint32 mOutputCommandDetailLevel;

	// simulation stats
	uint32 mTargetErrors;
	uint32 mTargetSims;
	uint32 mTargetExeTime;

	double mTotalExeTime;

	uint32 mTotalErrors;
	uint32 mTotalSims;
	uint32 mTotalDecodingTime;
	uint32 mTotalIterations;

	uint32 mTotalCorrectIterations;
	uint32 mTotalWrongIterations;

	uint32 mTotalCorrectDecodingTime;
	uint32 mTotalWrongDecodingTime;

	uint32 mTotalMLErrors;
	uint32 mTotalUndetectedErrors;

	//! The input is Eb/N0 (dB). This should be converted to std in AWGN or crossover probability in BSC
	void mTranslateEbN0()
	{
		mSTD = (double)sqrt(1.0 / pow(10.0, EbN0/10) / 2.0 / mRateOfCode);
		double temp = sqrt(2 * mRateOfCode * pow(10.0 , EbN0/10));
		mCrossOverProb = 0.5 - 0.5 * myerf(temp/sqrt(2.0));
		if(mChannelName == "AWGN")
			mChannelParameter = mSTD;
		else
			mChannelParameter = mCrossOverProb;
	}
	void mUpdateStats();
	void mClear()
	{
		mTotalErrors = 0;
		mTotalSims = 0;
		mTotalExeTime = 0;
		mTotalDecodingTime = 0;
		mTotalIterations = 0;
		mTotalCorrectIterations = 0;
		mTotalWrongIterations = 0;
		mTotalCorrectDecodingTime = 0;
		mTotalWrongDecodingTime = 0;
		mTotalMLErrors = 0;
		mTotalUndetectedErrors = 0;
	}

	void ShowTime(unsigned long secondes)
	{
	    int ss = secondes % 60;
	    int mn = (secondes / 60) % 60;
	    int hh = (secondes / 3600);
	    printf("%2.2dh%2.2d'%2.2d", hh, mn, ss);
	}

	void mOutputCommandFlush()
	{
//		if (mTotalSims % mOutputCommandInterval == 0)
//		{
			if(mOutputCommandDetailLevel == 1)
			{
				cout<<"TotalSims = "<<mTotalSims
					<<"  Time = "<<(double)mTotalExeTime/CLOCKS_PER_SEC<<"s  Total error = "<<mTotalErrors
					<<"  Avg # iteration = "<<mTotalIterations/(double)(mTotalSims)
					<<"  Undetected error = "<<mTotalUndetectedErrors
					<<"  ML error = "<<mTotalMLErrors
					<<"  Avg decode time = "<<(double)mTotalDecodingTime/CLOCKS_PER_SEC/ (double)(mTotalSims) <<"(s)"
					<<endl;
			}
			else if(mOutputCommandDetailLevel == 2)
			{
				cout<<"TotalSims = " << mTotalSims
					<<"  Time = "<<(double)mTotalExeTime/CLOCKS_PER_SEC<<"s"
					<<"  Total error = "<<mTotalErrors
					<<"  Undetected error = "<< mTotalUndetectedErrors
					<<"  ML error = "<<mTotalMLErrors
					<<"  Avg # iteration = "<<mTotalIterations/(double)(mTotalSims)
					<<"  Avg # iteration correct = "<<mTotalCorrectIterations/(double)(mTotalSims - mTotalErrors)
					<<"  Avg # iteration wrong = "<<mTotalWrongIterations/(double)(mTotalErrors)
					<<"  Avg decode time = "<<(double)mTotalDecodingTime/CLOCKS_PER_SEC/ (double)(mTotalSims) <<"(s)"
					<<"  Avg decode time correct = "<<(double)mTotalCorrectDecodingTime/CLOCKS_PER_SEC/ (double)(mTotalSims - mTotalErrors) <<"(s)"
					<<"  Avg decode time wrong = "<<(double)mTotalWrongDecodingTime/CLOCKS_PER_SEC/ (double)(mTotalErrors) <<"(s)"
					<<endl;
			}
//		}

			double tBER  = (mBlocklength * (double)mTotalErrors) / (mBlocklength * (double)mTotalSims);
		    double tFER  = (double)mTotalErrors / (double)mTotalSims;
		    double temps = (double)mTotalExeTime/(double)CLOCKS_PER_SEC;
		    double fpmn  = (60 * mTotalSims) / temps;
		    double bps   = ((double)fpmn * (double)mBlocklength) / 60.0 / 1000.0 / 1000.0;
		    double be_par_fe = mTotalUndetectedErrors;
		    printf("SNR = %.2f | BER =  %2.3e | FER =  %2.3e | BPS =  %2.2f | MATRICES = %10ld| FE = %d | ERR = %d | MlErr = %d | RUNTIME = ", EbN0, tBER, tFER, bps, mTotalSims, (int) mTotalErrors, (int) mTotalUndetectedErrors, mTotalUndetectedErrors);
		//    printf("SNR = %.2f | BER =  %2.3e | FER =  %2.3e | BPS =  %2.2f | MATRICES = %10ld| FE = %d | BE = %d | RUNTIME = ", Eb_N0, tBER, tFER, bps, counter->nb_processed_frames(), (int)counter->nb_fe(), (int)counter->nb_be());
		    ShowTime( temps );
		    printf("\n");
		    fflush(stdout);
	}

	void mOutputCommand()
	{
		if (mTotalSims % mOutputCommandInterval == 0)
		{
			if(mOutputCommandDetailLevel == 1)
			{
				cout<<"TotalSims = "<<mTotalSims
					<<"  Time = "<<(double)mTotalExeTime/CLOCKS_PER_SEC<< "s"
					<<"  Avg # iteration = "<<mTotalIterations/(double)(mTotalSims)
					<<"  Undetected error = "<<mTotalUndetectedErrors
					<<"  ML error = "<<mTotalMLErrors
					<<"  Avg decode time = "<<(double)mTotalDecodingTime/CLOCKS_PER_SEC/ (double)(mTotalSims) <<"(s)"
					<<endl;
			}
			else if(mOutputCommandDetailLevel == 2)
			{
				cout<<"TotalSims = "<<mTotalSims
					<<"  Time = "<<(double)mTotalExeTime/CLOCKS_PER_SEC<< "s"
					<<"  Total error = "<<mTotalErrors
					<<"  Undetected error = "<<mTotalUndetectedErrors
					<<"  ML error = "<<mTotalMLErrors
					<<"  Avg # iteration = "<<mTotalIterations/(double)(mTotalSims)
					<<"  Avg # iteration correct = "<<mTotalCorrectIterations/(double)(mTotalSims - mTotalErrors)
					<<"  Avg # iteration wrong = "<<mTotalWrongIterations/(double)(mTotalErrors)
					<<"  Avg decode time = "<<(double)mTotalDecodingTime/CLOCKS_PER_SEC/ (double)(mTotalSims) <<"(s)"
					<<"  Avg decode time correct = "<<(double)mTotalCorrectDecodingTime/CLOCKS_PER_SEC/ (double)(mTotalSims - mTotalErrors) <<"(s)"
					<<"  Avg decode time wrong = "<<(double)mTotalWrongDecodingTime/CLOCKS_PER_SEC/ (double)(mTotalErrors) <<"(s)"
					<<endl;
			}
		}
	}


	void mWriteFile()
	{
		int RoundsCount = 1;
		if (mTotalSims == mOutputFileInterval)
		{

				mOutputFileStream<<"TotalSims = "<<mOutputFileInterval<<"*"<<RoundsCount
					<<"  Time = "<<(double)mTotalExeTime/CLOCKS_PER_SEC<<"s"
					<<"  Avg # iteration = "<<mTotalIterations/(double)(mTotalSims)
					<<"  Undetected error = "<<mTotalUndetectedErrors
					<<"  ML error = "<<mTotalMLErrors
					<<"  Avg decode time = "<<(double)mTotalDecodingTime/CLOCKS_PER_SEC/ (double)(mTotalSims) <<"(s)"
					<<endl;
				RoundsCount++;
				mClear();
		}
	}
public:

	//! Return total number of errors
	uint32 GetTotalErrors(){return mTotalErrors;}
	//! Return total number of simulations
	uint32 GetTotalSims(){return mTotalSims;}
	//! Return time used for decoding
	uint32 GetTotalDecodingTime(){return mTotalDecodingTime;}
	//! Return total number of iterations
	uint32 GetTotalIterations(){return mTotalIterations;}
	//! Return total number of iterations for correct decoding events
	uint32 GetTotalCorrectIterations(){return mTotalCorrectIterations;}
	//! Return total number of iterations for erroneous decoding events
	uint32 GetTotalWrongIterations(){return mTotalWrongIterations;}
	//! Return time used for correct decoding events
	uint32 GetTotalCorrectDecodingTime(){return mTotalCorrectDecodingTime;}
	//! Return time used for erroneous decoding events
	uint32 GetTotalWrongDecodingTime(){return mTotalWrongDecodingTime;}
	//! Return total number of errors that are sure to be ML errors.
	/*!
		The ML error is calculated whenever there is a decoded codeword but is not the ML error. This is a lower bound for number of ML errors.
	*/
	uint32 GetTotalMLErrors(){return mTotalMLErrors;}
	//! Return total number of undetected errors. These are valid codeword but not the transmitted codeword.
	uint32 GetTotalUndetectedErrors(){return mTotalUndetectedErrors;}
	//! Return time for simulation. This include decoding time and time for channels, etc.
	uint32 GetTotalExeTime(){return mTotalExeTime;}
	//! Set channel's random number seed
	/*!
		\param seed Seed for channel.
	*/
	void SetChannelSeed(int seed){mChannel.SetSeed(seed);}
	//! Set simulation targets
	/*!
		\param tError Target errors. Stop simulation if more than this number of errors are collected.
		\param tSim Target simulations. Stop simulation if more than this number of simulations are performed.
		\param tExeTime Target time. Stop simulation if more than this amount of time is consumed.
	*/
	void SetTargets(uint32 tError, uint32 tSim, uint32 tExeTime){mTargetErrors = tError; mTargetSims = tSim; mTargetExeTime = tExeTime;}
	//! Set codeword.
	/*!
		Currently there is no encoding function. But you can choose a valid codeword to simulation if you have one other than the all zero codeword.
		\param filename Codeword file. The file should contain codeword. It should be 0,1 sequence, broken by space. e.g. 0 1 0 0 1 0
	*/
	void SetCodeword(string filename);
	//! Set command line output appearance.
	/*!
		\param simInterval Interval for each command line output.
		\param detaillevel Detail level for output. Only takes values 1 and 2.
	*/
	void SetCommandLineOutput(int simInterval, int detaillevel){mOutputCommandInterval = simInterval; mOutputCommandDetailLevel = detaillevel;}

	//! Set decoder name.
	/*!
	\param decoder Decoder name.
	 \sa Decoder
	*/
	void SetDecoder(string decoder){mDecoderName = decoder; mDecoder.SetDecoder(mDecoderName);}

	//! Set decoder parameters. Sufficient for ADMM LP.
	/*!
	 \param mxIt Maximum number of iterations. Default: 1000
	 \param feas ADMM stopping parameter. Default: 1e-6
	 \param p_mu ADMM step parameter mu. Default: 5
	 \param p_rho ADMM over relaxation parameter. Default: 1.8
	 \sa Decoder
	*/
	void SetDecoderParameters(int mxIt, double feas,  double p_mu, double p_rho)
	{
		mDecoder.SetParameters(mxIt,feas,p_mu,p_rho);
	}
	//! Set decoder parameters. Sufficient for ADMM L1 and ADMM L2
	/*!
	 \param mxIt Maximum number of iterations. Default: 1000
	 \param feas ADMM stopping parameter. Default: 1e-6
	 \param p_mu ADMM step parameter mu. Default: 5
	 \param p_rho ADMM over relaxation parameter. Default: 1.8
	 \param alp Penalty constant.
	 \sa Decoder
	*/
	void SetDecoderParameters(int mxIt, double feas,  double p_mu, double p_rho, double alp)
	{
		mDecoder.SetParameters(mxIt,feas,p_mu,p_rho); mDecoder.SetPenaltyConstant(alp);
	}
	//! Set decoder parameters. Sufficient for all ADMM algorithms
	/*!
	 \param mxIt Maximum number of iterations. Default: 1000
	 \param feas ADMM stopping parameter. Default: 1e-6
	 \param p_mu ADMM step parameter mu. Default: 5
	 \param p_rho ADMM over relaxation parameter. Default: 1.8
	 \param alp Penalty constant.
	 \param th Penalty thershold for entropy penalty
	 \sa Decoder
	*/
	void SetDecoderParameters(int mxIt, double feas,  double p_mu, double p_rho, double alp, double th)
	{
		mDecoder.SetParameters(mxIt,feas,p_mu,p_rho); mDecoder.SetPenaltyConstant(alp); mDecoder.SetEntropyParameter(th);
	}
	//! Set decoder parameters. Sufficient for BP decoding
	/*!
	 \param mxIt Maximum number of iterations. Default: 1000
	 \param feas ADMM stopping parameter. Default: 1e-6
	 \param p_mu ADMM step parameter mu. Default: 5
	 \param p_rho ADMM over relaxation parameter. Default: 1.8
	 \param nansat True if non-saturating BP is used
	 \param value Saturation value is original BP is used
	 \sa Decoder
	*/
	void SetDecoderParameters(int mxIt, double feas,  double p_mu, double p_rho, bool nonsat, double value)
	{
		mDecoder.SetParameters(mxIt,feas,p_mu,p_rho); mDecoder.SetBPSaturation(nonsat, value);
	}
	//! Set decoder parameters. Setting for normalized LLR for BSC
	/*!
	 \param mxIt Maximum number of iterations. Default: 1000
	 \param feas ADMM stopping parameter. Default: 1e-6
	 \param p_mu ADMM step parameter mu. Default: 5
	 \param p_rho ADMM over relaxation parameter. Default: 1.8
	 \param normal True if the LLRs from BSC are normalized to +1 or -1.
	 \sa Decoder
	*/
	void SetDecoderParameters(int mxIt, double feas,  double p_mu, double p_rho, bool normal)
	{
		mDecoder.SetParameters(mxIt,feas,p_mu,p_rho); mDecoder.SetBSCnormalized(normal);
	}
	//! Run simulations
	void RunSim();

	//! Constructor
	/*!
		\param blocklength Blocklength of the code
		\param nchecks Number of checks
		\param pcfilename Parity check matrix file name.
		\param channelname Channel name.
		\param eb_over_nzero Eb/N0 for the channel.
	*/
	Simulator(int blocklength, int nchecks, string pcfilename, string channelname, double eb_over_nzero)
		:mNoisySequence(blocklength)
		, mDecodedSequence(blocklength)
		, mDecoder(pcfilename, nchecks, blocklength)
		, mChannel(blocklength)
	{
		mChannelName = channelname;
		mBlocklength = blocklength;
		mNChecks = nchecks;
		mParityCheckFile = pcfilename;
		mRateOfCode = (double)(mBlocklength - mNChecks)/mBlocklength;
		EbN0 = eb_over_nzero;
		mTranslateEbN0();

		mNoisySequence.SetProperties(mChannelName, mChannelParameter);
		mChannel.SetChannel(mBlocklength, mChannelName, mChannelParameter);
		mDecoder.SetChannel(mChannelName);
		mInput = new int[mBlocklength];
		for(int i = 0; i < mBlocklength; i++)
			mInput[i] = 0;
	}
	~Simulator()
	{
		delete [] mInput;
	}
};


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
	mTotalIterations   += mDecodedSequence.Iteration;


}


void Simulator::RunSim()
{
	mTotalExeTime = 0;
	mClear(); // clear stats
	while(true)
	{

		mChannel.GenerateOutput(mInput, mNoisySequence);
		clock_t tStart = clock();
		mDecoder.Decode(mNoisySequence, mDecodedSequence);
		mTotalExeTime += clock() - tStart;
		mUpdateStats(); // update stats
		mOutputCommand();

		if(mTotalErrors >= mTargetErrors && mTargetErrors > 0)
			break;
		if(mTotalSims >= mTargetSims && mTargetSims > 0)
			break;
		if(mTotalExeTime >= mTargetExeTime && mTargetExeTime > 0)
			break;

	}
	mOutputCommandFlush( );
//	mTotalExeTime = clock() - tStart;
}



#endif /* SIMULATOR_HPP_ */
