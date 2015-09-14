/*
 * Decoder.hpp
 *
 *  Created on: 17 avr. 2015
 *      Author: legal
 */

#ifndef MatriceProjection_HPP_
#define MatriceProjection_HPP_
#include "./intel/RDTSC.hpp"
//!  Class for decoder
/*!
  Main program for ADMM decoders. Also included a sum-product decoder and a min-sum decoder. These two are not well optimized.
*/

inline bool the_compare_fx(const NODE & a, const NODE & b){return (a.value > b.value);}


template <class T = float, int N=32>
class MatriceProjection{
private:
	NODE  zSort   [32]; // DECLARATION STATIQUE POUR GAGNER DU TEMPS !
	NODE  zClip   [32]; // DECLARATION STATIQUE POUR GAGNER DU TEMPS !
	NODE  zBetaRep[32]; // DECLARATION STATIQUE POUR GAGNER DU TEMPS !
	float results [32]; // DECLARATION STATIQUE POUR GAGNER DU TEMPS !


	void function_sort_6VNs(NODE* v){
		const int nE = 6;
		int i, j;
		for (i = 1; i < nE; i++)
		{
			NODE tmp = v[i]; // ON RECUPERE LE NOEUD
			double value = tmp.value;
			for (j = i; (j >= 1) && (value > v[j - 1].value); j--)
			{
				v[j] = v[j - 1];
			}
			v[j] = tmp;
		}
	}

private:
        long long int t_sort1, t_sort2;
        int counter_sort1, counter_sort2;
	int counter_1;
	int counter_2;
	int counter_3;
	int counter_4;

public:
	MatriceProjection()
	{
                t_sort1       = 0;
                counter_sort1 = 0;
		t_sort2       = 0;
                counter_sort2 = 0;

		counter_1 = 0;
		counter_2 = 0;
		counter_3 = 0;
		counter_4 = 0;
	}

	~MatriceProjection()
	{
#ifdef PROFILE_ON

		printf("(PP) Counter sort value1 = %ld - time sort1 = %1.3f\n", counter_sort1, double(t_sort1)/counter_sort1);
		printf("(PP) Counter sort value2 = %ld - time sort2 = %1.3f \n", counter_sort2, double(t_sort2)/counter_sort2);
		printf("(PP) Counter 1 value = %d percent = %1.3f \n", counter_1, (float) counter_1 * 100/(counter_1 + counter_2 + counter_3 + counter_4));
		printf("(PP) Counter 2 value = %d percent = %1.3f \n", counter_2 , (float) counter_2 * 100/(counter_1 + counter_2 + counter_3 + counter_4));
		printf("(PP) Counter 3 value = %d percent = %1.3f \n", counter_3 , (float) counter_3 * 100/(counter_1 + counter_2 + counter_3 + counter_4));
		printf("(PP) Counter 4 value = %d percent = %1.3f \n", counter_4 , (float) counter_4 * 100/(counter_1 + counter_2 + counter_3 + counter_4));
#endif
	}

	T* mProjectionPPd(float v[], int length)
	{
		switch(length)
		{
		case 6:
			return mProjectionPPd<6>(v);
		case 7:
			return mProjectionPPd<7>(v);
		case 8:
			return mProjectionPPd<8>(v);
		case 9:
			return mProjectionPPd<9>(v);
		default:
			cout<<"ADMM not supported"<<endl;
			exit ( 0 );
		}
		return NULL;
	}


	template <int length=6>
	T* mProjectionPPd(float v[])
	{
		T zero_tol_offset = 1e-10;
		bool AllZero = true;
		bool AllOne  = true;

		#pragma unroll
		for(int i = 0; i < length; i++)
		{
#if 0
			AllZero = AllZero && (v[i].value <= 0);
			AllOne  = AllOne  && (v[i].value >  1);
#else
			if(v[i] >  0) AllZero = false;
			if(v[i] <= 1) AllOne  = false;
#endif
		}

		if(AllZero) // exit if the vector has all negative values, which means that the projection is all zero vector
		{
			#pragma unroll
			for(int i = 0;i < length; i++)
				results[i] = 0;
			#ifdef PROFILE_ON
				counter_1 += 1;
			#endif
			return results;
		}

		if (AllOne && (length%2 == 0)) // exit if the vector should be all one vector.
		{
			#pragma unroll
			for(int i = 0;i < length; i++)
				results[i] = 1;
			#ifdef PROFILE_ON
				counter_2 += 1;
			#endif
			return results;
		}

	//	NODE *zSort = new NODE[length]; // vector that save the sorted input vector

		#ifdef PROFILE_ON
               		const auto start_sort1= timer();
		#endif

		#pragma unroll
 		for(int i = 0; i < length; i++) // keep the original indices while sorting,
		{
			zSort[i].index = i;
			zSort[i].value = v[i];
		}

		sort(zSort, zSort + length, the_compare_fx);//sort input vector, decreasing order

                #ifdef PROFILE_ON
                	t_sort1 += timer() - start_sort1;
                	counter_sort1 ++ ;
		#endif

		int clip_idx = 0, zero_idx= 0;
		T constituent = 0;
	//	NODE *zClip = new NODE[length];
		#pragma unroll
		for(int i = 0;i < length; i++)// project on the [0,1]^d cube
		{
			zClip[i].value = min( max( zSort[i].value, 0.0f ), 1.0f);
			zClip[i].index = zSort[i].index;
			constituent   += zClip[i].value;
		}
		int r = (int)floor(constituent);
		if (r & 1)
			r--;
		// calculate constituent parity

		// calculate sum_Clip = $f_r^T z$
		T sum_Clip = 0;
		for(int i = 0; i < r+1; i++)
			sum_Clip += zClip[i].value;

		for(int i = r + 1; i < length; i++)
			sum_Clip -= zClip[i].value;


		if (sum_Clip <= r) // then done, return projection, beta = 0
		{
			#pragma unroll
			for(int i = 0; i < length; i++)
				results[zClip[i].index] = zClip[i].value;
	//		delete [] zSort;
	//		delete [] zClip;

			#ifdef PROFILE_ON
				counter_3 += 1;
			#endif
			return results;
		}

		T beta     = 0;
		T beta_max = 0;
		if (r + 2 <= length)
			beta_max = (zSort[r].value - zSort[r+1].value)/2; // assign beta_max
		else
			beta_max = zSort[r].value;

	//	NODE *zBetaRep = new NODE[length];

		// merge beta, save sort
		int left_idx  = r, right_idx = r + 1;
		int count_idx = 0;

		#ifdef PROFILE_ON
			const auto start_sort2= timer();
 		#endif

		while(count_idx < length)
		{
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

			T temp_a =  zSort[ left_idx].value - 1;
			T temp_b = -zSort[right_idx].value;
			if( (temp_a > temp_b) || (left_idx < 0) )
			{
				zBetaRep[count_idx].index = right_idx;
				zBetaRep[count_idx].value = temp_b;
				right_idx++;
				count_idx++;
			}
			else
			{
				zBetaRep[count_idx].index = left_idx;
				zBetaRep[count_idx].value = temp_a;
				left_idx --;
				count_idx++;
			}
		}

		#ifdef PROFILE_ON
	               t_sort2 += timer() - start_sort2;
        	       counter_sort2 ++ ;
		#endif

		#pragma unroll
		for(int i = 0; i < length; i++)
		{
			if (zSort[i].value > 1)
				clip_idx++;
			if (zSort[i].value >= 0 - zero_tol_offset)
				zero_idx++;
		}
		clip_idx--;

		bool first_positive = true;
		bool below_beta_max = true;
		int idx_start = 0;
		int idx_end   = 0;

		#pragma unroll
		for(int i = 0; i < length; i++)
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
		T active_sum = 0;

		#pragma unroll
		for(int i = 0;i < length; i++)
		{
			if(i > clip_idx && i <= r)
				active_sum += zSort[i].value;
			if(i > r && i < zero_idx)
				active_sum -= zSort[i].value;
		}
		T total_sum = 0;
		total_sum = active_sum + clip_idx + 1;

		int previous_clip_idx, previous_zero_idx;
		T previous_active_sum;
		bool change_pre     = true;
		previous_clip_idx   = clip_idx;
		previous_zero_idx   = zero_idx;
		previous_active_sum = active_sum;

		for(int i = idx_start; i <= idx_end; i++)
		{
			if(change_pre)
			{
				// save previous things
				previous_clip_idx   = clip_idx;
				previous_zero_idx   = zero_idx;
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


		for(int i = 0; i < length; i++)
		{
			if (i <= r)
			{
				results[zSort[i].index] = min(max(zSort[i].value - beta, 0.0f),1.0f);
			}
			else
			{
				results[zSort[i].index] = min(max(zSort[i].value + beta, 0.0f),1.0f);
			}

		}
	//	delete [] zSort; delete [] zClip; delete [] zBetaRep;
		#ifdef PROFILE_ON
			counter_4 += 1;
		#endif
		return results;
	}
};

#endif /* DECODER_HPP_ */
