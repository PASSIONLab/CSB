#include <cassert>
#include "bicsb.h"
#include "utility.h"

// Choose block size as big as possible given the following constraints
// 1) The bot array is addressible by IT
// 2) The parts of x & y vectors that a block touches fits into L2 cache [assuming a saxpy() operation]
// 3) There's enough parallel slackness for block rows (at least SLACKNESS * CILK_NPROC)
template <class NT, class IT>
void BiCsb<NT, IT>::Init(int workers, IT forcelogbeta)
{
	ispar = (workers > 1);
	IT roundrowup = nextpoweroftwo(m);
	IT roundcolup = nextpoweroftwo(n);

	// if indices are negative, highestbitset returns -1, 
	// but that will be caught by the sizereq below
	IT rowbits = highestbitset(roundrowup);
	IT colbits = highestbitset(roundcolup);
	bool sizereq;
	if (ispar)
	{
		sizereq = ((IntPower<2>(rowbits) > SLACKNESS * workers) 
			&& (IntPower<2>(colbits) > SLACKNESS * workers));
	}
	else
	{
		sizereq = ((rowbits > 1) && (colbits > 1));
	}

	if(!sizereq)
	{
		cerr << "Matrix too small for this library" << endl;
		return;
	}

	rowlowbits = rowbits-1;	
	collowbits = colbits-1;	
	IT inf = numeric_limits<IT>::max();
	IT maxbits = highestbitset(inf);

	rowhighbits = rowbits-rowlowbits;	// # higher order bits for rows (has at least one bit)
	colhighbits = colbits-collowbits;	// # higher order bits for cols (has at least one bit)
	if(ispar)
	{
		while(IntPower<2>(rowhighbits) < SLACKNESS * workers)
		{
			rowhighbits++;
			rowlowbits--;
		}
	}

	// calculate the space that suby occupies in L2 cache
	IT yL2 = IntPower<2>(rowlowbits) * sizeof(NT);
	while(yL2 > L2SIZE)
	{
		yL2 /= 2;
		rowhighbits++;
		rowlowbits--;
	}

	// calculate the space that subx occupies in L2 cache
	IT xL2 = IntPower<2>(collowbits) * sizeof(NT);
	while(xL2 > L2SIZE)
	{
		xL2 /= 2;
		colhighbits++;
		collowbits--;
	}
	
	// blocks need to be square for correctness (maybe generalize this later?) 
	while(rowlowbits+collowbits > maxbits)
	{
		if(rowlowbits > collowbits)
		{
			rowhighbits++;
			rowlowbits--;
		}
		else
		{
			colhighbits++;
			collowbits--;
		}
	}
	while(rowlowbits > collowbits)
	{
		rowhighbits++;
		rowlowbits--;
	}
	while(rowlowbits < collowbits)
	{
		colhighbits++;
		collowbits--;
	}
	assert (collowbits == rowlowbits);

	lowrowmask = IntPower<2>(rowlowbits) - 1;
	lowcolmask = IntPower<2>(collowbits) - 1;
	if(forcelogbeta != 0)
	{
		IT candlowmask  = IntPower<2>(forcelogbeta) -1;
		cout << "Forcing beta to "<< (candlowmask+1) << " instead of the chosen " << (lowrowmask+1) << endl;
		cout << "Warning : No checks are performed on the beta you have forced, anything can happen !" << endl;
		lowrowmask = lowcolmask = candlowmask;
		rowlowbits = collowbits = forcelogbeta;
		rowhighbits = rowbits-rowlowbits; 
		colhighbits = colbits-collowbits; 
	}
	else 
	{
		double sqrtn = sqrt(sqrt(static_cast<double>(m) * static_cast<double>(n)));
		IT logbeta = static_cast<IT>(ceil(log2(sqrtn))) + 2;
		if(rowlowbits > logbeta)
		{
			rowlowbits = collowbits = logbeta;
			lowrowmask = lowcolmask = IntPower<2>(logbeta) -1;
			rowhighbits = rowbits-rowlowbits;
	                colhighbits = colbits-collowbits;
		}
		cout << "Beta chosen to be "<< (lowrowmask+1) << endl;
	}
	highrowmask = ((roundrowup - 1) ^ lowrowmask);
	highcolmask = ((roundcolup - 1) ^ lowcolmask);
	
	// nbc = #{block columns} = #{blocks in any block row},  nbr = #{block rows)
	IT blcdimrow = lowrowmask + 1;
    	IT blcdimcol = lowcolmask + 1;
    	nbr = static_cast<IT>(ceil(static_cast<double>(m) / static_cast<double>(blcdimrow)));
    	nbc = static_cast<IT>(ceil(static_cast<double>(n) / static_cast<double>(blcdimcol)));
	
	blcrange = (lowrowmask+1) * (lowcolmask+1);	// range indexed by one block
	mortoncmp = MortonCompare<IT>(rowlowbits, collowbits, lowrowmask, lowcolmask);
}

// Partial template specialization for booleans
// Does not check cache considerations as this is mostly likely 
// to be used for gaxpy() with multiple rhs vectors (we don't know how many and what type at this point) 
template <class IT>
void BiCsb<bool,IT>::Init(int workers, IT forcelogbeta)
{
	ispar = (workers > 1);
	IT roundrowup = nextpoweroftwo(m);
	IT roundcolup = nextpoweroftwo(n);

	// if indices are negative, highestbitset returns -1, 
	// but that will be caught by the sizereq below
	IT rowbits = highestbitset(roundrowup);
	IT colbits = highestbitset(roundcolup);
	bool sizereq;
	if (ispar)
	{
		sizereq = ((IntPower<2>(rowbits) > SLACKNESS * workers) 
			&& (IntPower<2>(colbits) > SLACKNESS * workers));
	}
	else
	{
		sizereq = ((rowbits > 1) && (colbits > 1));
	}

	if(!sizereq)
	{
		cerr << "Matrix too small for this library" << endl;
		return;
	}

	rowlowbits = rowbits-1;	
	collowbits = colbits-1;	
	IT inf = numeric_limits<IT>::max();
	IT maxbits = highestbitset(inf);

	rowhighbits = rowbits-rowlowbits;	// # higher order bits for rows (has at least one bit)
	colhighbits = colbits-collowbits;	// # higher order bits for cols (has at least one bit)
	if(ispar)
	{
		while(IntPower<2>(rowhighbits) < SLACKNESS * workers)
		{
			rowhighbits++;
			rowlowbits--;
		}
	}

	// blocks need to be square for correctness (maybe generalize this later?) 
	while(rowlowbits+collowbits > maxbits)
	{
		if(rowlowbits > collowbits)
		{
			rowhighbits++;
			rowlowbits--;
		}
		else
		{
			colhighbits++;
			collowbits--;
		}
	}
	while(rowlowbits > collowbits)
	{
		rowhighbits++;
		rowlowbits--;
	}
	while(rowlowbits < collowbits)
	{
		colhighbits++;
		collowbits--;
	}
	assert (collowbits == rowlowbits);

	lowrowmask = IntPower<2>(rowlowbits) - 1;
	lowcolmask = IntPower<2>(collowbits) - 1;
	if(forcelogbeta != 0)
	{
		IT candlowmask  = IntPower<2>(forcelogbeta) -1;
		cout << "Forcing beta to "<< (candlowmask+1) << " instead of the chosen " << (lowrowmask+1) << endl;
		cout << "Warning : No checks are performed on the beta you have forced, anything can happen !" << endl;
		lowrowmask = lowcolmask = candlowmask;
		rowlowbits = collowbits = forcelogbeta;
		rowhighbits = rowbits-rowlowbits; 
		colhighbits = colbits-collowbits; 
	}
	else 
	{
		double sqrtn = sqrt(sqrt(static_cast<double>(m) * static_cast<double>(n)));
		IT logbeta = static_cast<IT>(ceil(log2(sqrtn))) + 2;
		if(rowlowbits > logbeta)
		{
			rowlowbits = collowbits = logbeta;
			lowrowmask = lowcolmask = IntPower<2>(logbeta) -1;
			rowhighbits = rowbits-rowlowbits;
	                colhighbits = colbits-collowbits;
		}
		cout << "Beta chosen to be "<< (lowrowmask+1) << endl;
	}
	highrowmask = ((roundrowup - 1) ^ lowrowmask);
	highcolmask = ((roundcolup - 1) ^ lowcolmask);
	
	// nbc = #{block columns} = #{blocks in any block row},  nbr = #{block rows)
	IT blcdimrow = lowrowmask + 1;
        IT blcdimcol = lowcolmask + 1;
        nbr = static_cast<IT>(ceil(static_cast<double>(m) / static_cast<double>(blcdimrow)));
        nbc = static_cast<IT>(ceil(static_cast<double>(n) / static_cast<double>(blcdimcol)));
	
	blcrange = (lowrowmask+1) * (lowcolmask+1);	// range indexed by one block
	mortoncmp = MortonCompare<IT>(rowlowbits, collowbits, lowrowmask, lowcolmask);
}


// Constructing empty BiCsb objects (size = 0) are not allowed.
template <class NT, class IT>
BiCsb<NT, IT>::BiCsb (IT size, IT rows, IT cols, int workers): nz(size),m(rows),n(cols)
{
	assert(nz != 0 && n != 0 && m != 0);
	Init(workers);

	num = (NT*) aligned_malloc( nz * sizeof(NT));
	bot = (IT*) aligned_malloc( nz * sizeof(IT));
	top = allocate2D<IT>(nbr, nbc+1);
}

// Partial template specialization for booleans
template <class IT>
BiCsb<bool, IT>::BiCsb (IT size, IT rows, IT cols, int workers): nz(size),m(rows),n(cols)
{
	assert(nz != 0 && n != 0 && m != 0);
	Init(workers);
	bot = (IT*) aligned_malloc( nz * sizeof(IT));
	top = allocate2D<IT>(nbr, nbc+1);
}

// copy constructor
template <class NT, class IT>
BiCsb<NT, IT>::BiCsb (const BiCsb<NT,IT> & rhs)
: nz(rhs.nz), m(rhs.m), n(rhs.n), blcrange(rhs.blcrange), nbr(rhs.nbr), nbc(rhs.nbc), 
rowhighbits(rhs.rowhighbits), rowlowbits(rhs.rowlowbits), highrowmask(rhs.highrowmask), lowrowmask(rhs.lowrowmask), 
colhighbits(rhs.colhighbits), collowbits(rhs.collowbits), highcolmask(rhs.highcolmask), lowcolmask(rhs.lowcolmask),
mortoncmp(rhs.mortoncmp), ispar(rhs.ispar)
{
	if(nz > 0)
	{
		num = (NT*) aligned_malloc( nz * sizeof(NT));
		bot = (IT*) aligned_malloc( nz * sizeof(IT));

		copy (rhs.num, rhs.num + nz, num);	
		copy (rhs.bot, rhs.bot + nz, bot);	
	}
	if ( nbr > 0)
	{
		top = allocate2D<IT>(nbr, nbc+1);
		for(IT i=0; i<nbr; ++i)
			copy (rhs.top[i], rhs.top[i] + nbc + 1, top[i]);
	}
}

// copy constructor for partial NT=boolean specialization
template <class IT>
BiCsb<bool, IT>::BiCsb (const BiCsb<bool,IT> & rhs)
: nz(rhs.nz), m(rhs.m), n(rhs.n), blcrange(rhs.blcrange), nbr(rhs.nbr), nbc(rhs.nbc), 
rowhighbits(rhs.rowhighbits), rowlowbits(rhs.rowlowbits), highrowmask(rhs.highrowmask), lowrowmask(rhs.lowrowmask), 
colhighbits(rhs.colhighbits), collowbits(rhs.collowbits), highcolmask(rhs.highcolmask), lowcolmask(rhs.lowcolmask),
mortoncmp(rhs.mortoncmp), ispar(rhs.ispar)
{
	if(nz > 0)
	{
		bot = (IT*) aligned_malloc( nz * sizeof(IT));
		copy (rhs.bot, rhs.bot + nz, bot);	
	}
	if ( nbr > 0)
	{
		top = allocate2D<IT>(nbr, nbc+1);
		for(IT i=0; i<nbr; ++i)
			copy (rhs.top[i], rhs.top[i] + nbc + 1, top[i]);
	}
}

template <class NT, class IT>
BiCsb<NT, IT> & BiCsb<NT, IT>::operator= (const BiCsb<NT, IT> & rhs)
{
	if(this != &rhs)		
	{
		if(nz > 0)	// if the existing object is not empty, make it empty
		{
			aligned_free(bot);
			aligned_free(num);
		}
		if(nbr > 0)
		{
			deallocate2D(top, nbr);
		}
		ispar 	= rhs.ispar;
		nz	= rhs.nz;
		n	= rhs.n;
		m   	= rhs.m;
		nbr 	= rhs.nbr;
		nbc 	= rhs.nbc;
		blcrange = rhs.blcrange;
		rowhighbits = rhs.rowhighbits;
		rowlowbits = rhs.rowlowbits;
		highrowmask = rhs.highrowmask;
		lowrowmask = rhs.lowrowmask;
		colhighbits = rhs.colhighbits;
		collowbits = rhs.collowbits;
		highcolmask = rhs.highcolmask;
		lowcolmask= rhs.lowcolmask;
		mortoncmp = rhs.mortoncmp;
		if(nz > 0)	// if the copied object is not empty
		{
			num = (NT*) aligned_malloc( nz * sizeof(NT));
			bot = (IT*) aligned_malloc( nz * sizeof(IT));
			copy (rhs.num, rhs.num + nz, num);	
			copy (rhs.bot, rhs.bot + nz, bot);	
		}
		if ( nbr > 0)
		{
			top = allocate2D<IT>(nbr, nbc+1);
			for(IT i=0; i<nbr; ++i)
				copy (rhs.top[i], rhs.top[i] + nbc + 1, top[i]);
		}
	}
	return *this;
}

template <class IT>
BiCsb<bool, IT> & BiCsb<bool, IT>::operator= (const BiCsb<bool, IT> & rhs)
{
	if(this != &rhs)		
	{
		if(nz > 0)	// if the existing object is not empty, make it empty
		{
			aligned_free(bot);
		}
		if(nbr > 0)
		{
			deallocate2D(top, nbr);
		}
		ispar 	= rhs.ispar;
		nz	= rhs.nz;
		n	= rhs.n;
		m   	= rhs.m;
		nbr 	= rhs.nbr;
		nbc 	= rhs.nbc;
		blcrange = rhs.blcrange;
		rowhighbits = rhs.rowhighbits;
		rowlowbits = rhs.rowlowbits;
		highrowmask = rhs.highrowmask;
		lowrowmask = rhs.lowrowmask;
		colhighbits = rhs.colhighbits;
		collowbits = rhs.collowbits;
		highcolmask = rhs.highcolmask;
		lowcolmask= rhs.lowcolmask;
		mortoncmp = rhs.mortoncmp;
		if(nz > 0)	// if the copied object is not empty
		{
			bot = (IT*) aligned_malloc( nz * sizeof(IT));
			copy (rhs.bot, rhs.bot + nz, bot);	
		}
		if ( nbr > 0)
		{
			top = allocate2D<IT>(nbr, nbc+1);
			for(IT i=0; i<nbr; ++i)
				copy (rhs.top[i], rhs.top[i] + nbc + 1, top[i]);
		}
	}
	return *this;
}

template <class NT, class IT>
BiCsb<NT, IT>::~BiCsb()
{
	if( nz > 0)
	{
		aligned_free((unsigned char*) num);
		aligned_free((unsigned char*) bot);
	}
	if ( nbr > 0)
	{	
		deallocate2D(top, nbr);	
	}
}

template <class IT>
BiCsb<bool, IT>::~BiCsb()
{
	if( nz > 0)
	{
		aligned_free((unsigned char*) bot);
	}
	if ( nbr > 0)
	{
		deallocate2D(top, nbr);		
	}
}

template <class NT, class IT>
BiCsb<NT, IT>::BiCsb (Csc<NT, IT> & csc, int workers, IT forcelogbeta):nz(csc.nz),m(csc.m),n(csc.n)
{
	typedef std::pair<IT, IT> ipair;
	typedef std::pair<IT, ipair> mypair;
	assert(nz != 0 && n != 0 && m != 0);
	if(forcelogbeta == 0)
		Init(workers);
	else
		Init(workers, forcelogbeta);	

	num = (NT*) aligned_malloc( nz * sizeof(NT));
	bot = (IT*) aligned_malloc( nz * sizeof(IT));
	top = allocate2D<IT>(nbr, nbc+1);
	mypair * pairarray = new mypair[nz];
	IT k = 0;
	for(IT j = 0; j < n; ++j)
	{
		for (IT i = csc.jc [j] ; i < csc.jc[j+1] ; ++i)	// scan the jth column
		{
			// concatenate the higher/lower order half of both row (first) index and col (second) index bits 
			IT hindex = (((highrowmask &  csc.ir[i] ) >> rowlowbits) << colhighbits)
										| ((highcolmask & j) >> collowbits);
			IT lindex = ((lowrowmask &  csc.ir[i]) << collowbits) | (lowcolmask & j) ;

			// i => location of that nonzero in csc.ir and csc.num arrays
			pairarray[k++] = mypair(hindex, ipair(lindex,i));
		}
	}
	sort(pairarray, pairarray+nz);	// sort according to hindex
	SortBlocks(pairarray, csc.num);
	delete [] pairarray;
}

template <class IT>
template <typename NT>	// to provide conversion from arbitrary Csc<> to specialized BiCsb<bool>
BiCsb<bool, IT>::BiCsb (Csc<NT, IT> & csc, int workers):nz(csc.nz),m(csc.m),n(csc.n)
{
	typedef std::pair<IT, IT> ipair;
	typedef std::pair<IT, ipair> mypair;
	assert(nz != 0 && n != 0 && m != 0);
	Init(workers);
	
	bot = (IT*) aligned_malloc( nz * sizeof(IT));
	top = allocate2D<IT>(nbr, nbc+1);
	mypair * pairarray = new mypair[nz];
	IT k = 0;
	for(IT j = 0; j < n; ++j)
	{
		for (IT i = csc.jc [j] ; i < csc.jc[j+1] ; ++i)	// scan the jth column
		{
			// concatenate the higher/lower order half of both row (first) index and col (second) index bits 
			IT hindex = (((highrowmask &  csc.ir[i] ) >> rowlowbits) << colhighbits)
										| ((highcolmask & j) >> collowbits);
			IT lindex = ((lowrowmask &  csc.ir[i]) << collowbits) | (lowcolmask & j) ;

			// i => location of that nonzero in csc.ir and csc.num arrays
			pairarray[k++] = mypair(hindex, ipair(lindex,i));
		}
	}
	sort(pairarray, pairarray+nz);	// sort according to hindex
	SortBlocks(pairarray);
	delete [] pairarray;
}

// Assumption: rowindices (ri) and colindices(ci) are "parallel arrays" sorted w.r.t. column index values
template <class NT, class IT>
BiCsb<NT, IT>::BiCsb (IT size, IT rows, IT cols, IT * ri, IT * ci, NT * val, int workers, IT forcelogbeta)
:nz(size),m(rows),n(cols)
{
	typedef std::pair<IT, IT> ipair;
	typedef std::pair<IT, ipair> mypair;
	assert(nz != 0 && n != 0 && m != 0);
	Init(workers, forcelogbeta);

	num = (NT*) aligned_malloc( nz * sizeof(NT));
	bot = (IT*) aligned_malloc( nz * sizeof(IT));
	top = allocate2D<IT>(nbr, nbc+1);
	mypair * pairarray = new mypair[nz];
	for(IT k = 0; k < nz; ++k)
	{
		// concatenate the higher/lower order half of both row (first) index and col (second) index bits 
		IT hindex = (((highrowmask &  ri[k] ) >> rowlowbits) << colhighbits)	| ((highcolmask & ci[k]) >> collowbits);	
		IT lindex = ((lowrowmask &  ri[k]) << collowbits) | (lowcolmask & ci[k]) ;

		// k is stored in order to retrieve the location of this nonzero in val array
		pairarray[k] = mypair(hindex, ipair(lindex, k));
	}
	sort(pairarray, pairarray+nz);	// sort according to hindex
	SortBlocks(pairarray, val);
	delete [] pairarray;
}

template <class IT>
BiCsb<bool, IT>::BiCsb (IT size, IT rows, IT cols, IT * ri, IT * ci, int workers, IT forcelogbeta)
:nz(size),m(rows),n(cols)
{
	typedef std::pair<IT, IT> ipair;
	typedef std::pair<IT, ipair> mypair;
	assert(nz != 0 && n != 0 && m != 0);
	Init(workers, forcelogbeta);

	bot = (IT*) aligned_malloc( nz * sizeof(IT));
	top = allocate2D<IT>(nbr, nbc+1);
	mypair * pairarray = new mypair[nz];
	for(IT k = 0; k < nz; ++k)
	{
		// concatenate the higher/lower order half of both row (first) index and col (second) index bits 
		IT hindex = (((highrowmask &  ri[k] ) >> rowlowbits) << colhighbits)	| ((highcolmask & ci[k]) >> collowbits);	
		IT lindex = ((lowrowmask &  ri[k]) << collowbits) | (lowcolmask & ci[k]) ;

		// k is stored in order to retrieve the location of this nonzero in val array
		pairarray[k] = mypair(hindex, ipair(lindex, k));
	}
	sort(pairarray, pairarray+nz);	// sort according to hindex
	SortBlocks(pairarray);
	delete [] pairarray;
}

template <class NT, class IT>
void BiCsb<NT, IT>::SortBlocks(pair<IT, pair<IT,IT> > * pairarray, NT * val)
{
 	typedef typename std::pair<IT, std::pair<IT, IT> > mypair;	
	IT cnz = 0;
	IT ldim = IntPower<2>(colhighbits);	// leading dimension (not always equal to nbc)
	for(IT i = 0; i < nbr; ++i)
	{
		for(IT j = 0; j < nbc; ++j)
		{
			top[i][j] = cnz;
			IT prevcnz = cnz; 
			vector< mypair > blocknz;
			while(cnz < nz && pairarray[cnz].first == ((i*ldim)+j) )	// as long as we're in this block
			{
				IT lowbits = pairarray[cnz].second.first;
				IT rlowbits = ((lowbits >> collowbits) & lowrowmask);
				IT clowbits = (lowbits & lowcolmask);
				IT bikey = BitInterleaveLow(rlowbits, clowbits);
				
				blocknz.push_back(mypair(bikey, pairarray[cnz++].second));
			}
			// sort the block into bitinterleaved order
			sort(blocknz.begin(), blocknz.end());

			for(IT k=prevcnz; k<cnz ; ++k)
			{
				bot[k] = blocknz[k-prevcnz].second.first;
				num[k] = val[blocknz[k-prevcnz].second.second];
			}
		}
		top[i][nbc] = cnz;  // hence equal to top[i+1][0] if i+1 < nbr
	}
	assert(cnz == nz);
}

template <class IT>
void BiCsb<bool, IT>::SortBlocks(pair<IT, pair<IT,IT> > * pairarray)
{
 	typedef pair<IT, pair<IT, IT> > mypair;	
	IT cnz = 0;
	IT ldim = IntPower<2>(colhighbits);	// leading dimension (not always equal to nbc)
	for(IT i = 0; i < nbr; ++i)
	{
		for(IT j = 0; j < nbc; ++j)
		{
			top[i][j] = cnz;
			IT prevcnz = cnz; 
			std::vector<mypair> blocknz;
			while(cnz < nz && pairarray[cnz].first == ((i*ldim)+j) )	// as long as we're in this block
			{
				IT lowbits = pairarray[cnz].second.first;
				IT rlowbits = ((lowbits >> collowbits) & lowrowmask);
				IT clowbits = (lowbits & lowcolmask);
				IT bikey = BitInterleaveLow(rlowbits, clowbits);
				
				blocknz.push_back(mypair(bikey, pairarray[cnz++].second));
			}
			// sort the block into bitinterleaved order
			sort(blocknz.begin(), blocknz.end());

			for(IT k=prevcnz; k<cnz ; ++k)
				bot[k] = blocknz[k-prevcnz].second.first;
		}
		top[i][nbc] = cnz;
	}
	assert(cnz == nz);
}


// Print stats to an ofstream object
template <class NT, class IT>
ofstream & BiCsb<NT, IT>::PrintStats(ofstream & outfile) const 
{
	if(nz == 0)
	{
		outfile << "## Matrix Doesn't have any nonzeros" <<endl;
		return outfile;
	}
	const IT ntop = nbr * nbc; 	

	outfile << "## Average block is of dimensions "<< lowrowmask+1 << "-by-" << lowcolmask+1 << endl;
	outfile << "## Number of real blocks is "<< ntop << endl;
	outfile << "## Row imbalance is " << RowImbalance(*this) << endl;
	outfile << "## Col imbalance is " << ColImbalance(*this) << endl;
	outfile << "## Block parallel calls is " << blockparcalls.get_value() << endl;
	
	std::vector<int> blocksizes(ntop);
	for(IT i=0; i<nbr; ++i)
	{
		for(IT j=0; j < nbc; ++j) 
		{
			blocksizes[i*nbc+j] = static_cast<int> (top[i][j+1]-top[i][j]);
		}
	}	
	sort(blocksizes.begin(), blocksizes.end());
	outfile<< "## Total nonzeros: "<< accumulate(blocksizes.begin(), blocksizes.end(), 0) << endl;

	outfile << "## Nonzero distribution (sorted) of blocks follows: \n" ;
	for(IT i=0; i< ntop; ++i)
	{	
		outfile << blocksizes[i] << "\n";
	}
	outfile << endl;
	return outfile;
}



// double* restrict a; --> No aliases for a[0], a[1], ...
// bstart/bend: block start/end index (to the top array)
template <class NT, class IT>
template <typename SR, typename RHS, typename LHS>
void BiCsb<NT, IT>::SubSpMV(IT * __restrict btop, IT bstart, IT bend, const RHS * __restrict x, LHS * __restrict suby) const
{
    IT * __restrict r_bot = bot;
    NT * __restrict r_num = num;

    for (IT j = bstart ; j < bend ; ++j)        // for all blocks inside that block row
    {
        // get higher order bits for column indices
        IT chi = (j << collowbits);
        const RHS * __restrict subx = &x[chi];
        for(IT k=btop[j]; k<btop[j+1]; ++k)
        {
            IT rli = ((r_bot[k] >> collowbits) & lowrowmask);
            IT cli = (r_bot[k] & lowcolmask);
            SR::axpy(r_num[k], subx[cli], suby[rli]);
        }
    }
}

// Partial boolean specialization on NT=bool
template <class IT>
template <typename SR, typename RHS, typename LHS>
void BiCsb<bool, IT>::SubSpMV(IT * __restrict btop, IT bstart, IT bend, const RHS * __restrict x, LHS * __restrict suby) const
{
    IT * __restrict r_bot = bot;
    for (IT j = bstart ; j < bend ; ++j)        // for all blocks inside that block row or chunk
    {
        // get higher order bits for column indices
        IT chi = (j << collowbits);
        const RHS * __restrict subx = &x[chi];
        for (IT k = btop[j] ; k < btop[j+1] ; ++k)    // for all nonzeros within ith block (expected =~ nnz/n = c)
        {
            IT rli = ((r_bot[k] >> collowbits) & lowrowmask);
            IT cli = (r_bot[k] & lowcolmask);
            SR::axpy(subx[cli], suby[rli]);        // suby [rli] += subx [cli]  where subx and suby are vectors.
        }
    }
}
