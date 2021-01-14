#ifndef _FRIENDS_H_
#define _FRIENDS_H_

#include <iostream>
#include <algorithm>
#include "bicsb.h"
#include "utility.h"
#include "timer.gettimeofday.c"

using namespace std;	

template <class NU, class IU>	
class BiCsb;


double prescantime;


#if (__GNUC__ == 4 && (__GNUC_MINOR__ < 7) )
#define emplace_back push_back
#endif

/**
  * Operation y = A*x+y on a semiring SR
  * A: a general CSB matrix (no specialization on booleans is necessary as this loop is independent of numerical values) 
  * x: a column vector or a set of column vectors (i.e. array of structs, array of std:arrays, etc))
  * SR::multiply() handles the multiple rhs and type promotions, etc. 
 **/
template <typename SR, typename NT, typename IT, typename RHS, typename LHS>
void bicsb_gespmv (const BiCsb<NT, IT> & A, const RHS * __restrict x, LHS * __restrict y)
{
	IT ysize = A.lowrowmask + 1;			// size of the output subarray (per block row - except the last)

#omp parallel for
    for (IT i = 0 ; i < A.nbr ; ++i)    // for all block rows of A
    {
        IT * btop = A.top [i];                       // get the pointer to this block row
        IT rhi = ((i << A.rowlowbits) & A.highrowmask);
        LHS * suby = &y[rhi];
        A.template SubSpMV<SR>(btop, 0, A.nbc, x, suby);
    }
}



#endif

