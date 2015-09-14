/*
 *
 *  Created on: 17 avr. 2015
 *      Author: legal
 * 
 * stable version: june 24th, 2015
 * still under modification
 */
#include <x86intrin.h>
#include <xmmintrin.h>

#ifndef DECODER_HPP_
#define DECODER_HPP_

#define PROFILE_ON

#include "MatriceProjection.hpp"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <chrono>
#include "./intel/IACA.hpp"
using namespace std;


#define _type_ float 

class Decoder{
public:
	//_type_* vector_before_proj;
	#define NEW_PROJECTION
	MatriceProjection<_type_, 32> mp;
private:

	string DecoderName;
	string ChannelName;
	string mRealDecoderName;

	double* mProjectionPPd   (NODE v[], int length);
	double* mProjectionPPd_d6(NODE v[]);
	double* mProjectionPPd_d7(NODE v[]);

    long long int t_vn = 0;
    long long int t_cn = 0;
    long long int t_pj = 0;
    long long int t_ex = 0;

	// parameters for decoders
	string mPCMatrixFileName;  /*!< parity check matrix file */
	int mNChecks, mBlocklength, mPCheckMapSize;
	int *CheckDegree;
	int *VariableDegree;
	double RateOfCode; /*!< code rate, k/n */

	TRIPLE *mPCheckMap; /*!< parity check mapping */
	TRIPLE *mParityCheckMatrix; /*!< parity check matrix */

	unsigned int* t_col;
	unsigned int* t_row;
	unsigned int* t_col1;
	unsigned int* t_row1;

	float *Lambda;
	float *zReplica;
	float *latestProjVector;

	// algorithm parameter
	double para_end_feas;   /*!< ADMM alg end tolerance, typically 1e-5 */
	double para_mu;         /*!< ADMM alg \mu, typically 5.5 */
	double para_rho;        /*!< ADMM over relaxation para, typically 1.8, 1.9 */
	int maxIteration;       /*!< max iteration */

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

	//
	//
	//

	template <int _mBlocklength = 2640, int _mNChecks = 1320, int _mPCheckMapSize = 7920>
    void ADMMDecoder_deg_6_3();

	//
	//
	//

	void ADMMDecoder_576x288();

	template <int _mBlocklength = 576, int _mNChecks = 288, int _mPCheckMapSize = 1824>
	void ADMMDecoder_deg_6_7_2_3_6();
    
	// data
	float *_LogLikelihoodRatio; /*!< log-likelihood ratio from received vector */
	float *OutputFromDecoder; /*!<soft information from decoder. could be pseudocodeword */

	double alpha; /*!< the constant for penalty term */

	//decoder stats
	unsigned long int mExeTime; /*!< exevution time */
	bool mAlgorithmConverge; /*!< true if the decoder converges*/
	int mIteration; /*!< number of iterations used for decoding */
	bool mValidCodeword; /*!< true if the output is a valid codeword */

public:
	//! Set decoder parameters.
	/*!
	 \param mxIt Maximum number of iterations. Default: 1000
	 \param feas ADMM stopping parameter. Default: 1e-6
	 \param p_mu ADMM step parameter mu. Default: 5
	 \param p_rho ADMM over relaxation parameter. Default: 1.8
	*/
	void SetParameters(int mxIt, double feas,  double p_mu, double p_rho);
	void SetDecoder(string decoder);
	void SetChannel(string channel){ChannelName = channel;}
	void SetPenaltyConstant(double alp){alpha = alp;}
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
	double Decode(NoisySeq &channeloutput, DecodedSeq &decoderesult);

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
	CheckDegree         = (int*)_mm_malloc(mNChecks * sizeof(int), 64);
	VariableDegree      = (int*)_mm_malloc(mBlocklength * sizeof(int), 64);
	OutputFromDecoder   = (float*)_mm_malloc(mBlocklength * sizeof(float), 64);
	_LogLikelihoodRatio = (float*)_mm_malloc(mBlocklength * sizeof(float), 64);

	mSetDefaultParameters();
	mLearnParityCheckMatrix();
	mPCheckMap = new TRIPLE[mPCheckMapSize];
	mParityCheckMatrix = new TRIPLE[mPCheckMapSize];
	mReadParityCheckMatrix();
	RateOfCode = (double)(mBlocklength -mNChecks)/mBlocklength;
	//vector_before_proj = (float*)_mm_malloc(64 * sizeof(float), 64);
        
	Lambda   = (float*)_mm_malloc(mPCheckMapSize * sizeof(float), 64);
	zReplica = (float*)_mm_malloc(mPCheckMapSize * sizeof(float), 64);
	latestProjVector = (float*)_mm_malloc(mPCheckMapSize * sizeof(float), 64);

	t_col              = (unsigned int*)_mm_malloc(mPCheckMapSize * sizeof(unsigned int), 64);
	int ptr = 0;
	for(int k = 0; k < mNChecks; k++)
	{
		for(int i = 0; i < mPCheckMapSize; i++)
		{
			if( mParityCheckMatrix[i].row == k )
			{
				t_col[ptr] = mParityCheckMatrix[i].col;
				ptr += 1;
			}

		}
	}

	t_col1             = (unsigned int*)_mm_malloc(mPCheckMapSize * sizeof(unsigned int), 64);
        ptr = 0;
	for(int k = 0; k < mPCheckMapSize; k++)
	{
		for(int i = 0; i < mPCheckMapSize; i++)
		{
			if( mPCheckMap[i].row == k )
			{
				t_col1[ptr] = mPCheckMap[i].col;
				ptr += 1;
                                break; // on en a fini avec le msg_k
			}

		}
	}

	t_row             = (unsigned int*)_mm_malloc(mPCheckMapSize * sizeof(unsigned int), 64);
        ptr = 0;
	for(int j = 0; j < mBlocklength; j++)
	{
		for(int i = 0; i < mPCheckMapSize; i++)
		{
			if( mPCheckMap[i].col == j )
			{
				t_row[ptr] = mPCheckMap[i].row;
				ptr += 1;
			}
		}
	}
	t_row1             = (unsigned int*)_mm_malloc(mPCheckMapSize * sizeof(unsigned int), 64);
        ptr = 0;
	for(int j = 0; j < mBlocklength; j++)
	{
		for(int i = 0; i < mPCheckMapSize; i++)
		{
			if( mPCheckMap[i].col == j )
			{
				t_row1[ptr] = mPCheckMap[i].row;
 				ptr += 1;
			}
		}
	}
}

Decoder::~Decoder()
{
#ifdef PROFILE_ON

    if( t_ex != 0 ){
        t_cn /= t_ex;
        t_vn /= t_ex;
        t_pj /= t_ex;
        long long int sum = t_cn + t_vn;
        double CN = (100.0 * t_cn) / (double)sum;
        double VN = (100.0 * t_vn) / (double)sum;
        double PJ = (100.0 * t_pj) / (double)t_cn;

        printf("cy_VN : %ld - cy_CN = %ld (cy_PJ = %ld) - t_ex: %ld \n", t_vn, t_cn, t_pj, t_ex);
	printf("VN : %1.3f - CN = %1.3f (PJ = %1.3f) \n", VN, CN, PJ);
	printf("cyVNpernode : %1.3f - cyCNpernode = %1.3f (cyPJpernode = %1.3f) \n", double(t_vn)/double(mBlocklength), double(t_cn)/double(mNChecks), double(t_pj)/double(mNChecks));

    }else{
    	printf("(WW) PROFILE WAS SWITCHED ON AND NO DATA ARE AVAILABLE !\n");
    }
#endif

	_mm_free ( CheckDegree );
	_mm_free ( VariableDegree );
	delete [] mPCheckMap;
	delete [] mParityCheckMatrix;
	_mm_free ( OutputFromDecoder );
	_mm_free (_LogLikelihoodRatio );
	//_mm_free ( vector_before_proj );
	_mm_free ( t_col );
	_mm_free ( t_row );
	_mm_free ( t_col1 );
	_mm_free ( t_row1 );
	_mm_free ( Lambda );
	_mm_free ( zReplica );
	_mm_free( latestProjVector );
}

//Decodersettings
void Decoder::SetParameters(int mxIt, double feas,  double p_mu, double p_rho)
{
	maxIteration = mxIt;
	para_end_feas = feas;
	para_mu = p_mu;
	para_rho = p_rho;
}
void Decoder::SetDecoder(string decoder)
{
	mRealDecoderName = decoder;
        DecoderName  = "ADMM";
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



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
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void Decoder::mLearnParityCheckMatrix()// learn the degree of the parity check matrix
{
	ifstream myfile (mPCMatrixFileName.c_str());
	int tempvalue = 0;
	int i = 0;
	string line;
	if (myfile.is_open())
	{
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
		   i++;
		}
		if (myfile.is_open())
			myfile.close();
	}
	else
	{
		cout << "Unable to open file";
		exit(0);
	}
	cout <<"mPCheckMapSize= "<< mPCheckMapSize<<endl;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void Decoder::mReadParityCheckMatrix()// read the parity check matrix, initialize message passing
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
				count++;
			}
			i++;
		}
	}
	if (myfile.is_open())
		myfile.close();
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



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



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



double Decoder::Decode(NoisySeq &channeloutput, DecodedSeq &decoderesult)
{
	mExeTime           = 0;
	mAlgorithmConverge = false;
	mValidCodeword     = false;

	mGenerateLLR(channeloutput);

	clock_t tStart = clock();
    auto start     = chrono::steady_clock::now();
    //IACA_START
    //
    //
    //
    if( (mBlocklength == 2640) && (mNChecks == 1320) ){
        ADMMDecoder_deg_6_3<2640, 1320, 7920>( );

    //
    //
    //
//
    } else if( (mBlocklength == 9216) && (mNChecks == 4608) ){
        ADMMDecoder_deg_6_3<9216, 4608, 27648>( );
    } else if( (mBlocklength == 4000) && (mNChecks == 2000) ){
        ADMMDecoder_deg_6_3<4000, 2000, 12000>( );

    //
    //
    //
    } else if( (mBlocklength == 576) && (mNChecks == 288) ){
    	//ADMMDecoder_576x288();

    	ADMMDecoder_deg_6_7_2_3_6<576, 288, 1824>( );

    } else if( (mBlocklength == 2304) && (mNChecks == 1152) ){
    	//ADMMDecoder_576x288();
    	ADMMDecoder_deg_6_7_2_3_6<2304, 1152, 7296>( );
        
    } else {
        ADMMDecoder_float();
    }
    auto end       = chrono::steady_clock::now();

    //IACA_END

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
		decoderesult.HardDecision[i] = floor(OutputFromDecoder[i] + 0.5);
	}
    auto diff      = end - start;
    return chrono::duration <double, milli> (diff).count();
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void show6(float* a)
{
	printf("__m256 : %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f\n", a[5], a[4], a[3], a[2], a[1], a[0]);
}

void show8(float* a)
{
	printf("__m256 : %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f\n", a[7], a[6], a[5], a[4], a[3], a[2], a[1], a[0]);
}

void show6(int* a)
{
	printf("__m256 : %d %d %d %d %d %d\n", a[5], a[4], a[3], a[2], a[1], a[0]);
}

void show8(int* a)
{
	printf("__m256 : %d %d %d %d %d %d %d %d\n", a[7], a[6], a[5], a[4], a[3], a[2], a[1], a[0]);
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "./decoders/decoder_576x288.hpp"

#include "./decoders/decoder_deg_6_3.hpp"

//template <int _mBlocklength, int _mNChecks, int _mPCheckMapSize>
//void Decoder::ADMMDecoder_deg_6_3()
//{
////	const auto ;
////	const auto ;
////	const auto ;
//	const auto _VariableDegree = 3;
//	const auto _CheckDegree    = 6;
//
////	float feas_tol = para_end_feas;
//	int maxIter    = maxIteration;
//
////	int numIter = 0;
//	float mu    = para_mu;
//	float rho   = para_rho;
//    const float un_m_rho = 1.0 - rho;
//
//    //////////////////////////////////////////////////////////////////////////////////////
//
//    const __m256 a = _mm256_set1_ps ( 0.0f );
//    const __m256 b = _mm256_set1_ps ( 0.5f );
//	#pragma  unroll
//	for( int j = 0; j < _mPCheckMapSize; j+=8 )
//    {
//		_mm256_store_ps(&Lambda  [j],         a);
//        _mm256_store_ps(&zReplica[j],         b);
//        _mm256_store_ps(&latestProjVector[j], b);
//    }
//    //////////////////////////////////////////////////////////////////////////////////////
//
//	for(int i = 0; i < maxIter; i++)
//	{
//		const float _amu_   = alpha/mu;
//		const float _2_amu_ = _amu_+ _amu_;
//		mIteration          = i + 1;
//
//	    //////////////////////////////////////////////////////////////////////////////////////
//
//		const auto t_mu   = _mm256_set1_ps ( mu      );
//		const auto t_amu  = _mm256_set1_ps ( _amu_   );
//		const auto t_2amu = _mm256_set1_ps ( _2_amu_ );
//		const auto un     = _mm256_set1_ps ( 1.0f   );
//		const auto zero   = _mm256_set1_ps ( 0.0f    );
//		const auto t_deg  = _mm256_set1_ps ( _VariableDegree );
//
//#ifdef PROFILE_ON
//		const auto start = timer();
//#endif
//
//		for(int j = 0; j < _mBlocklength; j+=8 )
//		{
//            auto ptr = j * _VariableDegree;
//            float M[8]      __attribute__((aligned(64)));
//            #pragma  unroll
//			for(int qq = 0; qq < 8; qq++)
//			{
//				M[qq] = 0.0f;
//				#pragma  unroll
//				for(int k = 0; k < _VariableDegree; k++)
//				{
//					M[qq] += (zReplica[ t_row[ptr] ] + Lambda[ t_row[ptr] ]);
//					ptr   += 1;
//				}
//			}
//			const auto m      = _mm256_loadu_ps( M );
//			const auto llr    = _mm256_loadu_ps( &_LogLikelihoodRatio[j] );
//			const auto t      = _mm256_div_ps(llr, t_mu);
//			const auto t1     = _mm256_sub_ps(m, t);
//			const auto x1     = _mm256_sub_ps(t1    , t_amu  );
//			const auto x2     = _mm256_sub_ps(t_deg , t_2amu );
//			const auto xx     = _mm256_div_ps(x1    , x2     );
//			const auto vMax   = _mm256_min_ps(  xx , un    );
//			const auto vMin   = _mm256_max_ps( vMax , zero);
//			_mm256_storeu_ps(&OutputFromDecoder[j], vMin);
//		}
//#ifdef PROFILE_ON
//        t_vn   += (timer() - start);
//#endif
//
//	    //////////////////////////////////////////////////////////////////////////////////////
////        int ptr     = 0;
//		//do projection for each check node
//		//
////		int CumSumCheckDegree = 0; // cumulative position of currect edge in factor graph
//        int allVerified = 0;
//
//        const auto  dot5      = _mm256_set1_ps(     0.5f );
//	    const auto  _rho      = _mm256_set1_ps(      rho );
//		const auto  _un_m_rho = _mm256_set1_ps( un_m_rho );
//		const auto mask_6     = _mm256_set_epi32(0x00000000, 0x00000000, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF);
//
//#ifdef PROFILE_ON
//		const auto starT = timer();
//#endif
//
//        for(int j = 0; j < _mNChecks; j++)
//		{
////            float vector_before_proj[8] __attribute__((aligned(64)));
//            int CumSumCheckDegree = j * _CheckDegree;
//
//            const auto offsets    = _mm256_loadu_si256  ((const __m256i*)&t_col1  [CumSumCheckDegree]);
//            const auto xpred      = _mm256_mask_i32gather_ps (zero, OutputFromDecoder, offsets, _mm256_castsi256_ps(mask_6), 4);
//			const auto synd       = _mm256_cmp_ps( xpred, dot5,   _CMP_GT_OS );
//			int test              = (_mm256_movemask_ps( synd ) & 0x3F);
//			const auto syndrom    = _mm_popcnt_u32( test );
//			const auto _Replica   = _mm256_loadu_ps(                &zReplica[CumSumCheckDegree]);
//			const auto _ambda     = _mm256_loadu_ps(                &Lambda  [CumSumCheckDegree]);
//			const auto v1         = _mm256_mul_ps  (   xpred,      _rho );
//			const auto v2         = _mm256_mul_ps  ( _Replica, _un_m_rho );
//			const auto v3         = _mm256_add_ps  ( v1, v2 );
//			const auto vect_proj  = _mm256_sub_ps  ( v3, _ambda );
//
//
//            //
//            // ON REALISE LA PROJECTION !!!
//            //
//            allVerified       += ( syndrom & 0x01 );
//#ifdef PROFILE_ON
//    		const auto START = timer();
//#endif
//			const auto latest    = _mm256_loadu_ps(&latestProjVector[CumSumCheckDegree]);
//			const auto different = _mm256_sub_ps ( vect_proj, latest );
//			const auto maskAbsol = _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF));
//			const auto absolute  = _mm256_and_ps ( different, maskAbsol );
//	        const auto seuilProj = _mm256_set1_ps( 1e-5f );
//	        const auto despass   = _mm256_cmp_ps( absolute, seuilProj, _CMP_GT_OS );
//	        int skip = (_mm256_movemask_ps( despass ) & 0x3F) == 0x00; // degree 6
//
//	        if( skip == false )
//	        {
//	            const auto _ztemp  = mp.projection_deg6( vect_proj );
////	        	printf("On skippe un calcul !!!\n");
//	    		const auto _ztemp1 = _mm256_sub_ps(_ztemp,    xpred );
//	    		const auto _ztemp2 = _mm256_sub_ps(_ztemp, _Replica );
//	    		const auto _ztemp3 = _mm256_mul_ps(_ztemp1, _rho);
//	    		const auto _ztemp4 = _mm256_mul_ps(_ztemp2, _un_m_rho);
//				const auto nLambda = _mm256_add_ps( _ambda,  _ztemp3 );
//				const auto mLambda = _mm256_add_ps( nLambda, _ztemp4 );
//	    		_mm256_maskstore_ps(&  Lambda[CumSumCheckDegree], mask_6, mLambda);
//	    		_mm256_maskstore_ps(&zReplica[CumSumCheckDegree], mask_6,  _ztemp);
//				_mm256_storeu_ps(&latestProjVector[CumSumCheckDegree], vect_proj);
//	        }
//#ifdef PROFILE_ON
//            t_pj   += (timer() - START);
//#endif
//		}
//
//#ifdef PROFILE_ON
//        t_cn   += (timer() - starT);
//#endif
//
//		if(allVerified == 0)
//		{
//			mAlgorithmConverge = true;
//			mValidCodeword     = true;
//			break;
//		}
//	}
//	t_ex += 1;
//}
//
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//deg_6_7_2_3_6
//void Decoder::ADMMDecoder_576x288()
//template <int _mBlocklength, int _mNChecks, int _mPCheckMapSize>
//void Decoder::ADMMDecoder_deg_6_7_2_3_6()
//{
//	int maxIter          = maxIteration;
//    const float mu       = 3.37309f;//para_mu;
//    const float rho      = 1.9f;    //para_rho;
//    const float un_m_rho = 1.0 - rho;
//    float tableau[]      = { 0.0f, 0.0f, 0.00001f, 2.00928f, 0.0f, 0.0f, 4.69438f };
//    float tableaX[7];
//
//    //
//    // ON PRECALCULE LES CONSTANTES
//    //
//    for (int i = 0; i < 7; i++)
//    {
//        tableaX[i] = tableau[ i ] / mu;
//    }
//
//    for (int i = 0;i < mPCheckMapSize; i++)
//	{
//		Lambda  [i] = 0;
//		zReplica[i] = 0.5;
//	}
//
//	for(int i = 0; i < maxIter; i++)
//	{
//        int ptr    = 0;
//		mIteration = i + 1;
//
//        //
//		// VN processing kernel
//		//
//		for (int j = 0; j < mBlocklength; j++)
//        {
//            const int degVn = VariableDegree[j];
//            float temp = 0;
//            for(int k = 0; k < degVn; k++)
//            {
//            	temp += (zReplica[ t_row[ptr] ] + Lambda[ t_row[ptr] ]);
//                ptr  +=1;
//            }
//            const float _amu_    = tableaX[ degVn ];
//            const float _2_amu_  = _amu_+ _amu_;
//            float llr  = _LogLikelihoodRatio[j];
//			float t    = temp - llr / mu;
//	        float xx   = (t  -  _amu_)/(degVn - _2_amu_);
//			float vMax = std::min(xx,   1.0f);
//			float vMin = std::max(vMax, 0.0f);
//			OutputFromDecoder[j] = vMin;
//        }
//
//		//
//		// CN processing kernel
//		//
//		int CumSumCheckDegree = 0; // cumulative position of currect edge in factor graph
//        int allVerified       = 0;
//		float vector_before_proj[8] __attribute__((aligned(64)));
//		for(int j = 0; j < mNChecks; j++)
//		{
//            if( CheckDegree[j] == 6 ){
//                int syndrom = 0;
//                #pragma unroll
//                for(int k = 0; k < 6; k++) // fill out the vector
//                {
//                    const auto ind         = CumSumCheckDegree + k;
//                    const auto xpred       = OutputFromDecoder[ t_col1[ind] ];
//                    syndrom               += (xpred - 0.5) > 0.0f;
//                    vector_before_proj[k]  = rho * xpred + un_m_rho * zReplica[ind] - Lambda[ind];
//                }
//                allVerified   += ( syndrom & 0x01 );
//                _type_ *ztemp  = mp.mProjectionPPd(vector_before_proj, 6);
//
//                for(int k = 0; k < 6; k++)
//                {
//                    const auto l   = CumSumCheckDegree + k;
//                    Lambda[l]      = Lambda[l] + (rho * (ztemp[k] - OutputFromDecoder[ t_col1[l] ]) + un_m_rho * (ztemp[k] - zReplica[l]));
//                }
//                #pragma unroll
//                for(int k = 0; k < 6; k++)      { zReplica[CumSumCheckDegree + k] = ztemp[k]; }
//                CumSumCheckDegree += CheckDegree[j];
//
//            }else if( CheckDegree[j] == 7 ){
//                int syndrom = 0;
//                #pragma unroll
//                for(int k = 0; k < 7; k++) // fill out the vector
//                {
//                    const auto ind         = CumSumCheckDegree + k;
//                    const auto xpred       = OutputFromDecoder[ t_col1[ind] ];
//                    syndrom               += (xpred - 0.5) > 0.0f;
//                    vector_before_proj[k]  = rho * xpred + un_m_rho * zReplica[ind] - Lambda[ind];
//                }
//                allVerified   += ( syndrom & 0x01 );
//                _type_ *ztemp  = mp.mProjectionPPd(vector_before_proj, 7);
//
//                for(int k = 0; k < CheckDegree[j]; k++)
//                {
//                    const auto l   = CumSumCheckDegree + k;
//                    Lambda[l]      = Lambda[l] + (rho * (ztemp[k] - OutputFromDecoder[ t_col1[l] ]) + un_m_rho * (ztemp[k] - zReplica[l]));
//                }
//                #pragma unroll
//                for(int k = 0; k < 7; k++) { zReplica[CumSumCheckDegree + k] = ztemp[k]; }
//                CumSumCheckDegree += CheckDegree[j];
//
//            }else{
//                exit( 0 );
//            }
//
//
//        }
//
//		if(allVerified == 0)
//		{
//			mAlgorithmConverge = true;
//			mValidCodeword     = true;
//			break;
//		}
//	}
//}
//
//
//

#include "decoders/decoder_deg_6_7_2_3_6.hpp"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "decoders/decoder_float.hpp"

//void Decoder::ADMMDecoder_float()
//{
//    int maxIter          = maxIteration;
//    const float mu       = para_mu;
//    const float rho      = para_rho;
//    const float un_m_rho = 1.0 - rho;
//    const float _amu_    = alpha/mu;
//    const float _2_amu_  = _amu_+ _amu_;
//
//    for (int i = 0;i < mPCheckMapSize; i++)
//    {
//        Lambda  [i] = 0;
//        zReplica[i] = 0.5;
//    }
//
//    for(int i = 0; i < maxIter; i++)
//    {
//        int ptr    = 0;
//        mIteration = i + 1;
//        //
//        // VN processing kernel
//        //
//        for (int j = 0; j < mBlocklength; j++)
//        {
//            float temp = 0;
//            for(int k = 0; k < VariableDegree[j]; k++)
//            {
//                temp += (zReplica[ t_row[ptr] ] + Lambda[ t_row[ptr] ]);
//                ptr  +=1;
//            }
//
//            float llr  = _LogLikelihoodRatio[j];
//            float t    = temp - llr / mu;
//            float deg  = (float)VariableDegree[j];
//            float xx   = (t  -  _amu_)/(deg - _2_amu_);
//            float vMax = std::min(xx,   1.0f);
//            float vMin = std::max(vMax, 0.0f);
//            OutputFromDecoder[j] = vMin;
//        }
//
//
//
//        //
//        // CN processing kernel
//        //
//        int CumSumCheckDegree = 0; // cumulative position of currect edge in factor graph
//        int allVerified       = 0;
//        float vector_before_proj[32] __attribute__((aligned(64)));
//        for(int j = 0; j < mNChecks; j++)
//        {
//            int syndrom = 0;
//#pragma unroll
//            for(int k = 0; k < CheckDegree[j]; k++) // fill out the vector
//            {
//                const auto ind         = CumSumCheckDegree + k;
//
//                //
//                // CALCUL TON VN...
//                //
//
//                const auto xpred       = OutputFromDecoder[ t_col1[ind] ];
//                syndrom               += (xpred - 0.5) > 0.0f;
//                vector_before_proj[k]  = rho * xpred + un_m_rho * zReplica[ind] - Lambda[ind];
//            }
//
//            allVerified += ( syndrom & 0x01 );
//            _type_ *ztemp  = mp.mProjectionPPd(vector_before_proj, CheckDegree[j]);
//
//            for(int k = 0; k < CheckDegree[j]; k++)
//            {
//                const auto l   = CumSumCheckDegree + k;
//                Lambda[l]      = Lambda[l] + (rho * (ztemp[k] - OutputFromDecoder[ t_col1[l] ]) + un_m_rho * (ztemp[k] - zReplica[l]));
//            }
//
//#pragma unroll
//            for(int k = 0; k < CheckDegree[j]; k++)
//            {
//                zReplica[CumSumCheckDegree + k] = ztemp[k];
//            }
//            CumSumCheckDegree += CheckDegree[j];
//
//            //
//            // MISE A JOUR DE VN
//            //
//        }
//
//        if(allVerified == 0)
//        {
//            mAlgorithmConverge = true;
//            mValidCodeword     = true;
//            break;
//        }
//    }
//}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



bool Decoder::ValidateCodeword()
{
	int ptr = 0;
	for(int k = 0; k < mNChecks; k++)
	{
		int syndrom = 0;
		int deg     = CheckDegree[k];
		for(int i = 0; i < deg; i++)
		{
			syndrom += (int) floor(OutputFromDecoder[ t_col[ptr] ] + 0.5);
			ptr += 1;
		}
		if( syndrom & 0x01 ) return false;
	}
	return true;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



bool Decoder::ValidateCodeword(int input[])
{
	int ptr = 0;
	for(int k = 0; k < mNChecks; k++)
	{
		int syndrom = 0;
		int deg     = CheckDegree[k];
		for(int i = 0; i < deg; i++)
		{
			syndrom += (int) floor(input[ t_col[ptr] ]);
			ptr += 1;
		}
		if( syndrom & 0x01 ) return false;
	}
	return true;
}

#endif /* DECODER_HPP_ */
