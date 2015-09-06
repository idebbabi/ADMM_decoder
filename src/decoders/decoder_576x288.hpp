void Decoder::ADMMDecoder_576x288()
{

	int maxIter          = maxIteration;
    const float mu       = 3.37309f;//para_mu;
    const float rho      = 1.9f;    //para_rho;
    const float un_m_rho = 1.0 - rho;
    float tableau[]      = { 0.0f, 0.0f, 0.00001f, 2.00928f, 0.0f, 0.0f, 4.69438f };
    float tableaX[7];

    //
    // ON PRECALCULE LES CONSTANTES
    //
    for (int i = 0; i < 7; i++)
    {
        tableaX[i] = tableau[ i ] / mu;
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

		#ifdef PROFILE_ON
				const auto start = timer();
		#endif

        //
		// VN processing kernel
		//
           
		for (int j = 0; j < mBlocklength; j++)
        {
            const int degVn = VariableDegree[j];
            float temp = 0;
            for(int k = 0; k < degVn; k++)
            {
            	temp += (zReplica[ t_row[ptr] ] + Lambda[ t_row[ptr] ]);
                ptr  +=1;
            }
            const float _amu_    = tableaX[ degVn ];
            const float _2_amu_  = _amu_+ _amu_;
            float llr  = _LogLikelihoodRatio[j];
			float t    = temp - llr / mu;
	        float xx   = (t  -  _amu_)/(degVn - _2_amu_);
			float vMax = std::min(xx,   1.0f);
			float vMin = std::max(vMax, 0.0f);
			OutputFromDecoder[j] = vMin;
        }

		#ifdef PROFILE_ON
				t_vn   += (timer() - start);
		#endif
		

		//
		// CN processing kernel
		//
		#ifdef PROFILE_ON
				const auto starT = timer();
		#endif

		int CumSumCheckDegree = 0; // cumulative position of currect edge in factor graph
        int allVerified       = 0;
		float vector_before_proj[8] __attribute__((aligned(64)));
		for(int j = 0; j < mNChecks; j++)
		{
            if( CheckDegree[j] == 6 ){
                int syndrom = 0;
                #pragma unroll
                for(int k = 0; k < 6; k++) // fill out the vector
                {
                    const auto ind         = CumSumCheckDegree + k;
                    const auto xpred       = OutputFromDecoder[ t_col1[ind] ];
                    syndrom               += (xpred - 0.5) > 0.0f;
                    vector_before_proj[k]  = (rho * xpred) + (un_m_rho * zReplica[ind]) - Lambda[ind];
                }

                allVerified   += ( syndrom & 0x01 );

    	    	//
    	    	// MEASURE OF THE PROJECTION EXECUTION TIME
    	    	//
                #ifdef PROFILE_ON
							const auto START = timer();
				#endif

				_type_ *ztemp  = mp.mProjectionPPd(vector_before_proj, 6);

    	    	//
    	    	// MEASURE OF THE PROJECTION EXECUTION TIME
    	    	//
				#ifdef PROFILE_ON
							t_pj   += (timer() - START);
				#endif
                
                for(int k = 0; k < 6; k++)
                {
                    const auto l   = CumSumCheckDegree + k;
                    Lambda[l]      = Lambda[l] + (rho * (ztemp[k] - OutputFromDecoder[ t_col1[l] ]) + un_m_rho * (ztemp[k] - zReplica[l]));
                }
                #pragma unroll
                for(int k = 0; k < 6; k++)      { zReplica[CumSumCheckDegree + k] = ztemp[k]; }
                CumSumCheckDegree += CheckDegree[j];

            }else if( CheckDegree[j] == 7 ){
                int syndrom = 0;
                #pragma unroll
                for(int k = 0; k < 7; k++) // fill out the vector
                {
                    const auto ind         = CumSumCheckDegree + k;
                    const auto xpred       = OutputFromDecoder[ t_col1[ind] ];
                    syndrom               += (xpred - 0.5) > 0.0f;
                    vector_before_proj[k]  = rho * xpred + un_m_rho * zReplica[ind] - Lambda[ind];
                }
                allVerified   += ( syndrom & 0x01 );

    	    	//
    	    	// MEASURE OF THE PROJECTION EXECUTION TIME
    	    	//
                #ifdef PROFILE_ON
							const auto START = timer();
				#endif

				_type_ *ztemp  = mp.mProjectionPPd(vector_before_proj, 7);

    	    	//
    	    	// MEASURE OF THE PROJECTION EXECUTION TIME
    	    	//
                #ifdef PROFILE_ON
							t_pj   += (timer() - START);
				#endif
                
                for(int k = 0; k < CheckDegree[j]; k++)
                {
                    const auto l   = CumSumCheckDegree + k;
                    Lambda[l]      = Lambda[l] + (rho * (ztemp[k] - OutputFromDecoder[ t_col1[l] ]) + un_m_rho * (ztemp[k] - zReplica[l]));
                }
                #pragma unroll
                for(int k = 0; k < 7; k++) { zReplica[CumSumCheckDegree + k] = ztemp[k]; }
                CumSumCheckDegree += CheckDegree[j];

            }else{
                exit( 0 );
            }
        }

    	//
    	// MEASURE OF THE CN LOOP EXECUTION TIME
    	//
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

	//
	// MEASURE OF THE NUMBER OF EXECUTION
	//


}

