/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//deg_6_7_2_3_6
//void Decoder::ADMMDecoder_576x288()
template <int _mBlocklength, int _mNChecks, int _mPCheckMapSize>
void Decoder::ADMMDecoder_deg_6_7_2_3_6()
{
	int maxIter          = maxIteration;
	float mu             = 5.5f; 
	float tableau[12]    = { 0.0f };

    if ((mBlocklength == 576) && (mNChecks == 288))
    {
     	mu          = 3.37309f;//penalty
        tableau[2]  = 0.00001f;
	tableau[3]  = 2.00928f;
	tableau[6]  = 4.69438f;

    }
    else if((mBlocklength == 2304) && (mNChecks == 1152) )
    {
    	mu          = 3.81398683f;//penalty
        tableau[2]  = 0.29669288f; 
	tableau[3]  = 0.46964023f;
	tableau[6]  = 3.19548154f;
    }
    else
    {
    	mu          = 5.5;//penalty
        tableau[2]  = 0.8f;
	tableau[3]  = 0.8f;
	tableau[6]  = 0.8f;
    }

    const float rho      = 1.9f;    //over relaxation parameter;
    const float un_m_rho = 1.0 - rho;
    const auto  _rho      = _mm256_set1_ps(      rho );
    const auto  _un_m_rho = _mm256_set1_ps( un_m_rho );
    float tableaX[12];

    //
    // ON PRECALCULE LES CONSTANTES
    //
	#pragma  unroll
    for (int i = 0; i < 7; i++)
    {
        tableaX[i] = tableau[ i ] / mu;
    }
	const auto t_mu    = _mm256_set1_ps ( mu );

	const auto t2_amu  = _mm256_set1_ps (        tableau[ 2 ] / mu   );
	const auto t3_amu  = _mm256_set1_ps (        tableau[ 3 ] / mu   );
	const auto t6_amu  = _mm256_set1_ps (        tableau[ 6 ] / mu   );

	const auto t2_2amu = _mm256_set1_ps ( 2.0f * tableau[ 2 ] / mu );
	const auto t3_2amu = _mm256_set1_ps ( 2.0f * tableau[ 3 ] / mu );
	const auto t6_2amu = _mm256_set1_ps ( 2.0f * tableau[ 6 ] / mu );

	const auto t2_deg  = _mm256_set1_ps ( 2.0f );
	const auto t3_deg  = _mm256_set1_ps ( 3.0f );
	const auto t6_deg  = _mm256_set1_ps ( 6.0f );

	const auto zero    = _mm256_set1_ps ( 0.0f );
	const auto un      = _mm256_set1_ps ( 1.0f );
    const __m256 a     = _mm256_set1_ps ( 0.0f );
    const __m256 b     = _mm256_set1_ps ( 0.5f );

    //////////////////////////////////////////////////////////////////////////////////////
	#pragma  unroll
	for( int j = 0; j < _mPCheckMapSize; j+=8 )
    {
	_mm256_store_ps(&Lambda  [j],         a);
        _mm256_store_ps(&zReplica[j],         b);
        _mm256_store_ps(&latestProjVector[j], b);
    }
    //////////////////////////////////////////////////////////////////////////////////////

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
		#pragma  unroll
		for (int j = 0; j < _mBlocklength; j++)
        {
            const int degVn = VariableDegree[j];
            float M[8] __attribute__((aligned(64)));

            if( degVn == 2 ){
#if 1
            	const int dVN = 2;
            	for(int qq = 0; qq < 8; qq++) 
		{
    				M[qq] = (zReplica[ t_row[ptr] ] + Lambda[ t_row[ptr] ]);      ptr   += 1;
    				#pragma  unroll
    				for(int k = 1; k < dVN; k++) 
				{
    					M[qq] += (zReplica[ t_row[ptr] ] + Lambda[ t_row[ptr] ]); ptr   += 1;
    				}
    		}
    		const auto m      = _mm256_loadu_ps( M );
    		const auto llr    = _mm256_loadu_ps( &_LogLikelihoodRatio[j] );
    		const auto t1     = _mm256_sub_ps(m, _mm256_div_ps(llr, t_mu));
    		const auto xx     = _mm256_div_ps(_mm256_sub_ps(t1, t2_amu), _mm256_sub_ps(t2_deg, t2_2amu));
    		const auto vMin   = _mm256_max_ps(_mm256_min_ps(xx, un) , zero);
    		_mm256_storeu_ps(&OutputFromDecoder[j], vMin);
    		j += 7;
#else
            	const int degVN = 2;
                float temp = (zReplica[ t_row[ptr] ] + Lambda[ t_row[ptr] ]);
		#pragma unroll
		for(int k = 1; k < degVN; k++)
			temp += (zReplica[ t_row[ptr + k] ] + Lambda[ t_row[ptr + k] ]);
		ptr  += degVN;
	        const float _amu_    = tableaX[ degVN ];
	        const float _2_amu_  = _amu_+ _amu_;
	        const float llr  = _LogLikelihoodRatio[j];
	        const float t    = temp - llr / mu;
	        const float xx   = (t  -  _amu_)/(degVn - _2_amu_);
	        const float vMax = std::min(xx,   1.0f);
	        const float vMin = std::max(vMax, 0.0f);
		OutputFromDecoder[j] = vMin;
#endif
            }else if( degVn == 3 ){
#if 1
            	const int dVN = 3;
            	for(int qq = 0; qq < 8; qq++) 
		{
    			M[qq] = (zReplica[ t_row[ptr] ] + Lambda[ t_row[ptr] ]);      ptr   += 1;
    			#pragma  unroll
    			for(int k = 1; k < dVN; k++) 
			{
    				M[qq] += (zReplica[ t_row[ptr] ] + Lambda[ t_row[ptr] ]); ptr   += 1;
    			}
    		}
    		const auto m      = _mm256_loadu_ps( M );
    		const auto llr    = _mm256_loadu_ps( &_LogLikelihoodRatio[j] );
    		const auto t1     = _mm256_sub_ps(m, _mm256_div_ps(llr, t_mu));
    		const auto xx     = _mm256_div_ps(_mm256_sub_ps(t1, t3_amu), _mm256_sub_ps(t3_deg, t3_2amu));
    		const auto vMin   = _mm256_max_ps(_mm256_min_ps(xx, un) , zero);
    		_mm256_storeu_ps(&OutputFromDecoder[j], vMin);
    		j += 7;
#else
    		const int degVN = 3;
                float temp = (zReplica[ t_row[ptr] ] + Lambda[ t_row[ptr] ]);
		#pragma unroll
		for(int k = 1; k < degVN; k++)
			temp += (zReplica[ t_row[ptr + k] ] + Lambda[ t_row[ptr + k] ]);
		ptr  += degVN;
	        const float _amu_    = tableaX[ degVN ];
	        const float _2_amu_  = _amu_+ _amu_;
	        const float llr  = _LogLikelihoodRatio[j];
	        const float t    = temp - llr / mu;
	        const float xx   = (t  -  _amu_)/(degVn - _2_amu_);
	        const float vMax = std::min(xx,   1.0f);
	        const float vMin = std::max(vMax, 0.0f);
		OutputFromDecoder[j] = vMin;
#endif
		}else if( degVn == 6 ){
#if 1
            	const int dVN = 6;
            	for(int qq = 0; qq < 8; qq++) 
		{
    			M[qq] = (zReplica[ t_row[ptr] ] + Lambda[ t_row[ptr] ]);      ptr   += 1;
    			#pragma  unroll
    			for(int k = 1; k < dVN; k++) 
			{
    				M[qq] += (zReplica[ t_row[ptr] ] + Lambda[ t_row[ptr] ]); ptr   += 1;
    			}
    		}
    		const auto m      = _mm256_loadu_ps( M );
    		const auto llr    = _mm256_loadu_ps( &_LogLikelihoodRatio[j] );
    		const auto t1     = _mm256_sub_ps(m, _mm256_div_ps(llr, t_mu));
    		const auto xx     = _mm256_div_ps(_mm256_sub_ps(t1, t6_amu), _mm256_sub_ps(t6_deg, t6_2amu));
    		const auto vMin   = _mm256_max_ps(_mm256_min_ps(xx, un) , zero);
    		_mm256_storeu_ps(&OutputFromDecoder[j], vMin);
    		j += 7;
#else
    		const int degVN = 6;
                float temp = (zReplica[ t_row[ptr] ] + Lambda[ t_row[ptr] ]);
		#pragma unroll
		for(int k = 1; k < degVN; k++)
			temp += (zReplica[ t_row[ptr + k] ] + Lambda[ t_row[ptr + k] ]);
		ptr  += degVN;
	        const float _amu_    = tableaX[ degVN ];
	        const float _2_amu_  = _amu_+ _amu_;
	        const float llr  = _LogLikelihoodRatio[j];
	        const float t    = temp - llr / mu;
	        const float xx   = (t  -  _amu_)/(degVn - _2_amu_);
	        const float vMax = std::min(xx,   1.0f);
	        const float vMin = std::max(vMax, 0.0f);
		OutputFromDecoder[j] = vMin;
#endif
            }
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
	int CumSumCheckDegree = 0; // cumulative position of currect edge in factor graph
        int allVerified       = 0;
	float vector_before_proj[8] __attribute__((aligned(64)));

        const auto zero    = _mm256_set1_ps ( 0.0f    );
        const auto mask_6  = _mm256_set_epi32(0x00000000, 0x00000000, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF);
        const auto mask_7  = _mm256_set_epi32(0x00000000, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF);
        const auto  dot5   = _mm256_set1_ps(     0.5f );

    	//
    	// MEASURE OF THE CN EXECUTION TIME
    	//
		#ifdef PROFILE_ON
				const auto starT = timer();
		#endif

    	const auto seuilProj = _mm256_set1_ps( 1e-5f );
        for(int j = 0; j < _mNChecks; j++)
		{
            if( CheckDegree[j] == 6 ){
            	const int  cDeg6       = 0x3F;
                const auto offsets    = _mm256_loadu_si256  ((const __m256i*)&t_col1  [CumSumCheckDegree]);
                const auto xpred      = _mm256_mask_i32gather_ps (zero, OutputFromDecoder, offsets, _mm256_castsi256_ps(mask_6), 4);
    		const auto synd       = _mm256_cmp_ps( xpred, dot5,   _CMP_GT_OS );
    		int test              = (_mm256_movemask_ps( synd ) & cDeg6);  // deg 6
    		const auto syndrom    = _mm_popcnt_u32( test );
    		const auto _Replica   = _mm256_loadu_ps( &zReplica[CumSumCheckDegree]);
    		const auto _ambda     = _mm256_loadu_ps( &Lambda  [CumSumCheckDegree]);
    		const auto v1         = _mm256_mul_ps  (xpred,      _rho );
    		const auto v2         = _mm256_mul_ps  ( _Replica, _un_m_rho );
    		const auto v3         = _mm256_add_ps  ( v1, v2 );
    		const auto vect_proj  = _mm256_sub_ps  ( v3, _ambda );

                //
                // ON REALISE LA PROJECTION !!!
                //
                allVerified       += ( syndrom & 0x01 );

    	    	//
    	    	// MEASURE OF THE PROJECTION EXECUTION TIME
    	    	//
                #ifdef PROFILE_ON
					const auto START = timer();
		#endif
    		const auto latest    = _mm256_loadu_ps(&latestProjVector[CumSumCheckDegree]);
    		const auto different = _mm256_sub_ps ( vect_proj, latest );
    		const auto maskAbsol = _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF));
    		const auto absolute  = _mm256_and_ps ( different, maskAbsol );
    	        const auto despass   = _mm256_cmp_ps( absolute, seuilProj, _CMP_GT_OS );
    	        int skip = (_mm256_movemask_ps( despass ) & cDeg6) == 0x00; // degree 6

    	        if( skip == false )
    	        {
    	        	const auto _ztemp  = mp.projection_deg6( vect_proj );
    	    		const auto _ztemp1 = _mm256_sub_ps(_ztemp,    xpred );
    	    		const auto _ztemp2 = _mm256_sub_ps(_ztemp, _Replica );
	    	    	const auto _ztemp3 = _mm256_mul_ps(_ztemp1, _rho);
	    	    	const auto _ztemp4 = _mm256_mul_ps(_ztemp2, _un_m_rho);
    			const auto nLambda = _mm256_add_ps( _ambda,  _ztemp3 );
    			const auto mLambda = _mm256_add_ps( nLambda, _ztemp4 );
    	    		_mm256_maskstore_ps(&  Lambda[CumSumCheckDegree],         mask_6,   mLambda);
    	    		_mm256_maskstore_ps(&zReplica[CumSumCheckDegree],         mask_6,    _ztemp);
    	        }
	    	_mm256_maskstore_ps(&latestProjVector[CumSumCheckDegree], mask_6, vect_proj);

    	    	//
    	    	// MEASURE OF THE PROJECTION EXECUTION TIME
    	    	//
    	        #ifdef PROFILE_ON
					t_pj   += (timer() - START);
		#endif
                CumSumCheckDegree += 6;

            }else if( CheckDegree[j] == 7 )
	    {
            	const int  cDeg7       = 0x7F;
                const auto offsets    = _mm256_loadu_si256  ((const __m256i*)&t_col1  [CumSumCheckDegree]);
                const auto xpred      = _mm256_mask_i32gather_ps (zero, OutputFromDecoder, offsets, _mm256_castsi256_ps(mask_7), 4);
    		const auto synd       = _mm256_cmp_ps( xpred, dot5,   _CMP_GT_OS );
    		const int  test       = (_mm256_movemask_ps( synd ) & cDeg7); // deg 7
    		const auto syndrom    = _mm_popcnt_u32( test );
    		const auto _Replica   = _mm256_loadu_ps( &zReplica[CumSumCheckDegree]);
    		const auto _ambda     = _mm256_loadu_ps( &Lambda  [CumSumCheckDegree]);
    		const auto v1         = _mm256_mul_ps  ( xpred,    _rho );
    		const auto v2         = _mm256_mul_ps  ( _Replica, _un_m_rho );
    		const auto v3         = _mm256_add_ps  ( v1, v2 );
    		const auto vect_proj  = _mm256_sub_ps  ( v3, _ambda );

                //
                // ON REALISE LA PROJECTION !!!
                //
                allVerified         += ( syndrom & 0x01 );

    	    	//
    	    	// MEASURE OF THE PROJECTION EXECUTION TIME
    	    	//
                #ifdef PROFILE_ON
					const auto START = timer();
		#endif
    		const auto latest    = _mm256_loadu_ps(&latestProjVector[CumSumCheckDegree]);
    		const auto different = _mm256_sub_ps ( vect_proj, latest );
    		const auto maskAbsol = _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF));
    		const auto absolute  = _mm256_and_ps ( different, maskAbsol );
    	        const auto despass   = _mm256_cmp_ps( absolute, seuilProj, _CMP_GT_OS );
    	        int skip = (_mm256_movemask_ps( despass ) & cDeg7) == 0x00; // degree 7

    	        if( skip == false )
    	        {
			const auto _ztemp  = mp.projection_deg7( vect_proj );
    	    		const auto _ztemp1 = _mm256_sub_ps(_ztemp,    xpred );
    	    		const auto _ztemp2 = _mm256_sub_ps(_ztemp, _Replica );
    	    		const auto _ztemp3 = _mm256_mul_ps(_ztemp1, _rho);
    	    		const auto _ztemp4 = _mm256_mul_ps(_ztemp2, _un_m_rho);
    			const auto nLambda = _mm256_add_ps( _ambda,  _ztemp3 );
    			const auto mLambda = _mm256_add_ps( nLambda, _ztemp4 );
    	    		_mm256_maskstore_ps(&  Lambda        [CumSumCheckDegree], mask_7,   mLambda);
    	    		_mm256_maskstore_ps(&zReplica        [CumSumCheckDegree], mask_7,    _ztemp);
    	        }
	    	_mm256_maskstore_ps(&latestProjVector[CumSumCheckDegree], mask_7, vect_proj);

    	    	//
    	    	// MEASURE OF THE PROJECTION EXECUTION TIME
    	    	//
    	        #ifdef PROFILE_ON
							t_pj   += (timer() - START);
		#endif

                CumSumCheckDegree += 7;

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
		//FILE *ft=fopen("time.txt","a");
		//fprintf(ft,"%d \n", t_cn/t_ex);
		//fprintf(ft,"%d %d %d \n", t_cn, t_vn, t_pj);
		//fclose(ft);
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
//	#ifdef PROFILE_ON
//		t_ex += 1;
//	#endif

}

