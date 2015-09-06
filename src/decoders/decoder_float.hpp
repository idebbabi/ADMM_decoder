

void Decoder::ADMMDecoder_float()
{
	float feas_tol       = para_end_feas;
	int maxIter          = maxIteration;
	int numIter          = 0;
	const float rho      = 1.9f;//over relaxation 
	const float un_m_rho = 1.0 - rho;
	float mu             = 5.5f; 
	float tableau[12]     = { 0.0f};

    if ((mBlocklength == 1152) && (mNChecks == 288))
    {
	mu          = 3.43369433f;//penalty
        tableau[2]  = 0.00001f; 
	tableau[3]  = 1.01134229f; 
	tableau[6]  = 4.43497372f;
    } 
    else if((mBlocklength == 1944) && (mNChecks == 972) )//dv 2,3,6,11
    {
	mu          = 3.43369433f;//penalty
        tableau[2]  = 0.8f;
	tableau[3]  = 0.8f;
	tableau[6]  = 0.8f;
	tableau[11] = 0.8f;
    }
   else 
    {
	mu          = 5.5;//penalty// works for dv 2, 3, 6
        tableau[2]  = 0.8f;
	tableau[3]  = 0.8f;
	tableau[6]  = 0.8f;
	tableau[11] = 0.8f;
    }
   	 for (int i = 0;i < mPCheckMapSize; i++)
	{
		Lambda  [i] = 0;
		zReplica[i] = 0.5;
	}
	for(int i = 0; i < maxIter; i++)
	{
        	int ptr    = 0;
		mIteration = i + 1;
    	//
    	// MEASURE OF THE VN EXECUTION TIME
    	//
	#ifdef PROFILE_ON
		const auto start = timer();
	#endif
       
	//
	// VN processing kernel
	//
	for (int j = 0; j < mBlocklength; j++)
        {
            float temp = 0;
            for(int k = 0; k < VariableDegree[j]; k++)
            {
            	temp += (zReplica[ t_row[ptr] ] + Lambda[ t_row[ptr] ]);
                ptr  +=1;
            }

	const auto val = tableau[ VariableDegree[j] ];
	const float _amu_    = val / mu;
	const float _2_amu_  = _amu_+ _amu_;

	float llr  = _LogLikelihoodRatio[j];
	float t    = temp - llr / mu;
	float deg  = (double)VariableDegree[j];
	float xx   = (t  -  _amu_)/(deg - _2_amu_);
	float vMax = std::min(xx,   1.0f);
	float vMin = std::max(vMax, 0.0f);
	OutputFromDecoder[j] = vMin;
        }

        //
    	// MEASURE OF THE VN EXECUTION TIME
    	//
	#ifdef PROFILE_ON
		t_vn   += (timer() - start);
	#endif
	//
	// CN processing kernel
	//
	// MEASURE OF THE CN EXECUTION TIME
    	//
	#ifdef PROFILE_ON
		const auto starT = timer();
	#endif

	int CumSumCheckDegree = 0; // cumulative position of currect edge in factor graph
        int allVerified = 0;
	float vector_before_proj[32] __attribute__((aligned(64)));
	for(int j = 0; j < mNChecks; j++)
	{
		int syndrom = 0;
		#pragma unroll
		    for(int k = 0; k < CheckDegree[j]; k++) // fill out the vector
		    {
		    	const auto ind         = CumSumCheckDegree + k;
 		    	const auto xpred       = OutputFromDecoder[ t_col1[ind] ];
		    	syndrom               += (xpred - 0.5) > 0.0f;
		    	vector_before_proj[k]  = rho * xpred + un_m_rho * zReplica[ind] - Lambda[ind];
		    }
		    allVerified += ( syndrom & 0x01 );

		#ifdef PROFILE_ON
				const auto START = timer();
		#endif
			_type_ *ztemp  = mp.mProjectionPPd(vector_before_proj, CheckDegree[j]);

		#ifdef PROFILE_ON
			t_pj   += (timer() - START);
		#endif

		for(int k = 0; k < CheckDegree[j]; k++) 
	    	{
		const auto l   = CumSumCheckDegree + k;
		Lambda[l]      = Lambda[l] + (rho * (ztemp[k] - OutputFromDecoder[ t_col1[l] ]) + un_m_rho * (ztemp[k] - zReplica[l]));
		}
		#pragma unroll
		for(int k = 0; k < CheckDegree[j]; k++)
		{
                	zReplica[CumSumCheckDegree + k] = ztemp[k];
		}
		CumSumCheckDegree += CheckDegree[j];
 	}
	        #ifdef PROFILE_ON
				t_cn   += (timer() - starT);
		#endif

		#ifdef PROFILE_ON
			t_ex += 1;
		#endif

		if(allVerified == 0)
		{
			mAlgorithmConverge = true;
			mValidCodeword     = true;
			break;
		}
	}
}

