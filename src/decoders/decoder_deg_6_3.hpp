template <int _mBlocklength, int _mNChecks, int _mPCheckMapSize>
void Decoder::ADMMDecoder_deg_6_3( )
{
	const auto _VariableDegree = 3;
	const auto _CheckDegree    = 6;
	int maxIter  = maxIteration;
	float mu     = 3.0f;
	float rho    = 1.9f;
    const float un_m_rho = 1.0 - rho;

    //////////////////////////////////////////////////////////////////////////////////////

    const __m256 a = _mm256_set1_ps ( 0.0f );
    const __m256 b = _mm256_set1_ps ( 0.5f );
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
		const float _amu_   = alpha/mu;
		const float _2_amu_ = _amu_+ _amu_;
		mIteration          = i + 1;

	    //////////////////////////////////////////////////////////////////////////////////////

		const auto t_mu   = _mm256_set1_ps ( mu      );
		const auto t_amu  = _mm256_set1_ps ( _amu_   );
		const auto t_2amu = _mm256_set1_ps ( _2_amu_ );
		const auto un     = _mm256_set1_ps ( 1.0f   );
		const auto zero   = _mm256_set1_ps ( 0.0f    );
		const auto t_deg  = _mm256_set1_ps ( _VariableDegree );

    	//
    	// MEASURE OF THE VN EXECUTION TIME
    	//
		#ifdef PROFILE_ON
				const auto start = timer();
		#endif

		for(int j = 0; j < _mBlocklength; j+=8 )
		{

            auto ptr = j * _VariableDegree;
            float M[8]      __attribute__((aligned(64)));
            #pragma  unroll
			for(int qq = 0; qq < 8; qq++)
			{
				M[qq] = 0.0f;
				#pragma  unroll
				for(int k = 0; k < _VariableDegree; k++)
				{
					M[qq] += (zReplica[ t_row[ptr] ] + Lambda[ t_row[ptr] ]);
					ptr   += 1;
				}
			}
			const auto m      = _mm256_loadu_ps( M );
			const auto llr    = _mm256_loadu_ps( &_LogLikelihoodRatio[j] );
			const auto t      = _mm256_div_ps(llr, t_mu);
			const auto t1     = _mm256_sub_ps(m, t);
			const auto x1     = _mm256_sub_ps(t1    , t_amu  );
			const auto x2     = _mm256_sub_ps(t_deg , t_2amu );
			const auto xx     = _mm256_div_ps(x1    , x2     );
			const auto vMax   = _mm256_min_ps(  xx , un    );
			const auto vMin   = _mm256_max_ps( vMax , zero);
			_mm256_storeu_ps(&OutputFromDecoder[j], vMin);

		}

    	//
    	// MEASURE OF THE VN EXECUTION TIME
    	//

		#ifdef PROFILE_ON
		t_vn   += (timer() - start);
		#endif

	    //////////////////////////////////////////////////////////////////////////////////////
        int allVerified       = 0;
        const auto  dot5      = _mm256_set1_ps(     0.5f );
	    const auto  _rho      = _mm256_set1_ps(      rho );
		const auto  _un_m_rho = _mm256_set1_ps( un_m_rho );
		const auto mask_6     = _mm256_set_epi32(0x00000000, 0x00000000, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF);

    	//
    	// MEASURE OF THE CN EXECUTION TIME
    	//
		#ifdef PROFILE_ON
				const auto starT = timer();
		#endif

        for(int j = 0; j < _mNChecks; j++)
		{
            int CumSumCheckDegree = j * _CheckDegree;
            const auto offsets    = _mm256_loadu_si256  ((const __m256i*)&t_col1  [CumSumCheckDegree]);
            const auto xpred      = _mm256_mask_i32gather_ps (zero, OutputFromDecoder, offsets, _mm256_castsi256_ps(mask_6), 4);
			const auto synd       = _mm256_cmp_ps( xpred, dot5,   _CMP_GT_OS );
			int test              = (_mm256_movemask_ps( synd ) & 0x3F);
			const auto syndrom    = _mm_popcnt_u32( test );
			const auto _Replica   = _mm256_loadu_ps(                &zReplica[CumSumCheckDegree]);
			const auto _ambda     = _mm256_loadu_ps(                &Lambda  [CumSumCheckDegree]);
			const auto v1         = _mm256_mul_ps  (   xpred,      _rho );
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
/*
            #ifdef PROFILE_ON
    			const auto START = timer();
			#endif
*/
			const auto latest    = _mm256_loadu_ps(&latestProjVector[CumSumCheckDegree]);
			const auto different = _mm256_sub_ps ( vect_proj, latest );
			const auto maskAbsol = _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF));
			const auto absolute  = _mm256_and_ps ( different, maskAbsol );
	        const auto seuilProj = _mm256_set1_ps( 1e-6f );
	        const auto despass   = _mm256_cmp_ps( absolute, seuilProj, _CMP_GT_OS );
	        int skip = (_mm256_movemask_ps( despass ) & 0x3F) == 0x00; // degree 6

	        if( skip == false )
	        {
            #ifdef PROFILE_ON
    			const auto START = timer();
			#endif
	            const auto _ztemp  = mp.projection_deg6( vect_proj );
  	#ifdef PROFILE_ON
				t_pj   += (timer() - START);
			#endif
	    		const auto _ztemp1 = _mm256_sub_ps(_ztemp,    xpred );
	    		const auto _ztemp2 = _mm256_sub_ps(_ztemp, _Replica );
	    		const auto _ztemp3 = _mm256_mul_ps(_ztemp1, _rho);
	    		const auto _ztemp4 = _mm256_mul_ps(_ztemp2, _un_m_rho);
				const auto nLambda = _mm256_add_ps( _ambda,  _ztemp3 );
				const auto mLambda = _mm256_add_ps( nLambda, _ztemp4 );
	    		_mm256_maskstore_ps(&  Lambda[CumSumCheckDegree], mask_6, mLambda);
	    		_mm256_maskstore_ps(&zReplica[CumSumCheckDegree], mask_6,  _ztemp);
				_mm256_storeu_ps(&latestProjVector[CumSumCheckDegree], vect_proj);
	        }

	        //
	    	// MEASURE OF THE PROJECTION EXECUTION TIME
	    	//
/*
  	#ifdef PROFILE_ON
				t_pj   += (timer() - START);
			#endif
*/
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
/*
	#ifdef PROFILE_ON
		t_ex += 1;

	#endif
*/
}
