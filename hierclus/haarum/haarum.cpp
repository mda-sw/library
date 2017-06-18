//////////////////////////////////////////////////////////////////////////////////////////////
//
// haarum.cpp 23/03/2005
//
//////////////////////////////////////////////////////////////////////////////////////////////

// INCLUDES //////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <limits.h>
#include <stdlib.h> 
#include <memory.h>
#include <math.h>
#include <assert.h>

// DEFINES ///////////////////////////////////////////////////////////////////////////////////

#define EPS (double)1.e-10

// CLASSES ///////////////////////////////////////////////////////////////////////////////////

template <class T>
class vector_t
{
private:
	T*		buffer;
	int		size;

public:

	// constructors/ destructors
	vector_t(): buffer(NULL), size(0) {}
	vector_t(const vector_t<T>& v);
	vector_t(int _size);
	vector_t(const T *v, int _size);	

	~vector_t() { clear(); }

	
	// members
	void clear(void);
	void set_size(int new_size);

	int length(void) const { return(size); }
	void push_back(T elem) { set_size(size+1); buffer[size - 1] = elem; }


	T& operator[] (int i) 
	{
		i--;
		assert((i >= 0) && (i < size)); 
		return buffer[i]; 
	}

	T& operator[] (int i) const
	{
		i--;
		assert((i >= 0) && (i < size)); 
		return buffer[i]; 
	} 

	vector_t<T>& operator= (const vector_t<T>& v);

};

template <class T>
vector_t<T>::vector_t(const vector_t<T>& v)
{
	size = v.length();
	buffer = new T[size];	
	memcpy(buffer, &v[1], sizeof(T) * size);
}

template <class T>
vector_t<T>::vector_t(int _size)
{
	size = _size;
	buffer = new T[size];
	memset(buffer, 0, sizeof(T) * size);
}

template <class T>
vector_t<T>::vector_t(const T *v, int _size)
{
	size = _size;
	buffer = new T[size];	
	memcpy(buffer, v, sizeof(T) * size);
}


template <class T> 
void vector_t<T>::clear(void)
{
	if (buffer) 
	{
		delete [] buffer;
		buffer = NULL;
		size = 0;
	}		
}

template <class T> 
void vector_t<T>::set_size(int new_size)
{
	T *temp;

	assert(new_size >= 1);

	if (size == 0)
	{
		buffer = new T[new_size];
		size = new_size;
	}
	else
	{
		temp = new T[size];
		memcpy(temp, buffer, sizeof(T) * size);
		delete [] buffer;
		buffer = new T[new_size];

		if (new_size < size)
			memcpy(buffer, temp, sizeof(T) * new_size);
		else
			memcpy(buffer, temp, sizeof(T) * size);

		size = new_size;

		delete [] temp;
	}
	
} // end set range


template <class T>
vector_t<T>& vector_t<T>::operator= (const vector_t<T>& v) 
{ 
	set_size(v.length()); 
	memcpy(buffer, (void*)&v[1], sizeof(T) * v.length());
	return *this; 
} 

// MATRIX CLASS ///////////////


template<class T> 
class matrix_t
{
private:
	T		*buffer;
	int		num_rows;
	int		num_cols;

	void clear(void);

	matrix_t():buffer(NULL), num_rows(0), num_cols(0){};

public:
	
	matrix_t(const matrix_t<T>& m);
	matrix_t(int n, int m);
	matrix_t(const T *v, int n, int m);

	~matrix_t() { clear(); }

	
	int nrow(void) const  { return(num_rows); }
	int ncol(void) const { return(num_cols); }

	T *get_vector(void) const { return(buffer); }

	T& operator[](int i) 
	{ 
		return buffer[i]; 
	} 

	T& operator[](int i) const
	{ 
		return buffer[i]; 
	} 

	T& operator () ( int row, int col  )
	{
		row--;
		col--;
		assert((row >= 0) && (row < num_rows));
		assert((col >= 0) && (col < num_cols));
		
		return buffer[row + col*num_rows];
	}

	T& operator () ( int row, int col ) const 
	{
		row--;
		col--;
		assert((row >= 0) && (row < num_rows));
		assert((col >= 0) && (col < num_cols));
		
		return buffer[row + col*num_rows];
	}

	matrix_t<T>& operator= (const matrix_t<T>& m);
	

};


template<class T> 
matrix_t<T>::matrix_t(const matrix_t<T>& m)
{
	
	num_rows = m.nrow();
	num_cols = m.ncol(); 
	
	buffer = new T[num_rows * num_cols];

	memcpy(buffer, &m[0], sizeof(T) * m.nrow() *  m.ncol());
}

template<class T> 
matrix_t<T>::matrix_t(int n, int m)
{
	num_rows = n;
	num_cols = m;
	
	buffer = new T[num_rows * num_cols];

	memset(buffer, 0, sizeof(T) * num_rows * num_cols);
}


template<class T> 
matrix_t<T>::matrix_t(const T *v, int n, int m)
{
	num_rows = n;
	num_cols = m;
	
	buffer = new T[num_rows * num_cols];

	memcpy(buffer, v, sizeof(T) * n *  m);
	
}

template<class T>
void matrix_t<T>::clear(void)
{
	if (buffer)
	{
		delete [] buffer;
		buffer = NULL;
		num_rows = 0;
		num_cols = 0;
		
	}	
}

template<class T>
matrix_t<T>& matrix_t<T>::operator= (const matrix_t<T>& m) 
{ 
	clear();

	num_rows = m.nrow();
	num_cols = m.ncol(); 
	
	buffer = new T[num_rows * num_cols];

	memcpy(buffer, &m[0], sizeof(T) * m.nrow() *  m.ncol());

	return *this; 
} 

// TYPES /////////////////////////////////////////////////////////////////////////////////////

typedef vector_t<double> fvector;
typedef vector_t<int>   ivector;

typedef matrix_t<double> fmatrix;
typedef matrix_t<int>   imatrix;

// FUNCTIONS /////////////////////////////////////////////////////////////////////////////////

int compare(const void *p1, const void *p2)
{
	int *num1, *num2;

	num1 = (int*)p1;
	num2 = (int*)p2;

	if (*num1 <  *num2) 
		return -1;
  	else
  	if (*num1 > *num2) 
		return  1;
	else
  		return 0;
} // end compare

//////////////////////////////////////////////////////////////////////////////////////////////

extern "C" {

void haarum(double *pxout,   // output
			double *psmooth, // output
			double *pdetail, // output
			const double *px, 
			const int *pnrow, 
			const int *pncol, 
			const int *pmemberships,  
			const double* pthreshold)
{
	int n, m, i, j;
	int clusval;
	int aine, benjamin;
	int waszero, nowzero;
	int obs;
	int wasclusnum, clusnum, k;
	

	n = (*pnrow);
	m = (*pncol);

	// Reformat membership so that they have the following format.
	// In the following, column 8 contains singletons
	// Column 7 indicates cluster 9 agglomerating 1 and 8
	// Column 6 entails agglomeration of 1, 8 and 7 into new cluster 10
	// And so on, until column 1 shows one cluster with all observesions
	// this data was obtained from median hierarchic clusterring of
	// the first 8 observesins in Fishers's iris data (oroginaly 4-dimensional)
	// 15 14 10 10 10 10 9 1
	// 15 14 13  3  3  3 3 3 
	// 15 14 13 12 11  4 4 4
	// 15 14 13 12 11  6 6 6
	// 15 14 10 10 10 10 9 8
	// 15  2  2  2 	2  2 2 2
	// 15 14 13 12  5  5 5 5
	// 15 14 10 10 10 10 7 7

	imatrix memberships(pmemberships, n, n);
	imatrix a(memberships);

	//ivector changedvals(n);

	for (j=(n-1); j >= 1; j--)
	{
		clusval = INT_MAX;
		for (i=1; i<=n; i++)
		{
			if (memberships(i,j) != memberships(i,j+1))
			{
				if (memberships(i,j) < clusval)
					clusval = memberships(i,j);
			} // end if			
				
		} // end for i

		for (i=1; i<=n; i++)
		{
			if (memberships(i,j) == clusval) 
				a(i,j) = 2*n-j;
			else
				a(i,j) = a(i,j+1);
		} // end for i	
		
	} /// end for j


	// First, forward transform

	fmatrix ident(px, n, m);
	fmatrix xout(px, n, m);

	fmatrix smooth(m, n-1);
	fmatrix detail(m, n-1);

	
	// In the succession of signal smooth and detail( vectors, we have this scheme:
	// Col n-1 contains in particular details of new cluster, n+1        = q1
	// Col n-2                                                n+2        = q2
	// Col n-3                                                n+3        = q3
	// ...
	// Col n-(n-2)                                            n+(n-2)    = q(n-2)
	// Col n-(n-1)                                            n+(n-1)    = q(n-1)

	// We must also cater for the left and right branching in the dendrogram
	// We call left branch: aine (elder),
	// and right branch: benjamin (younger)

	fvector sub1(m);
	fvector sub2(m);

	for (j=1; j<=(n-1); j++)				// Set of hierarchy levels, 1 through n-1
	{
			
		// Get cluster menbers, proceeding through q1, q2, ... q(n-1)

		// Find subclusters: these are only two; and they will be found
		// in the column that follows. We do this to characterize aine, benjamin
		
		ivector nextcolmmbrs;

		for (i=1; i<=n; i++)
		{
			if (a(i, n-j) == n+j)
			{
				nextcolmmbrs.push_back(a(i,n-j+1));		
			}
		}

		qsort(&nextcolmmbrs[1], nextcolmmbrs.length(), sizeof(int),  compare);	// Always, by convention, aine < benjamin
		
		aine     = nextcolmmbrs[1];
		benjamin = nextcolmmbrs[nextcolmmbrs.length()];
		
		nextcolmmbrs.clear();

		// Aine is a singleton
		if (aine <= n)
		{
			for (i=1; i<=m; i++)
				sub1[i] = ident(aine,i);
		} // end if
		else
		{
			for (i=1; i<=m; i++)
				sub1[i] = smooth(i,2*n-aine);			
		}

		// Benjamin is a singleton
		if (benjamin <= n)
		{
			for (i=1; i<=m; i++)
				sub2[i] = ident(benjamin,i);			
		} // end if
		else
		{
			for (i=1; i<=m; i++)
				sub2[i] = smooth(i,2*n-benjamin);		
		}
		

		for (i=1; i<=m; i++)
		{
			smooth(i,n-j) = 0.5 * (sub1[i] + sub2[i]);  // Smooth s(qj)
			detail(i,n-j) = 0.5 * (sub1[i] - sub2[i]);	// Detail d(qj)
			
		} // end for i

		//mmbrs.clear();

	} // end for j-loop (hierarchy levels)
	
	sub1.clear();
	sub2.clear();


	// Filter
	waszero = 0;
	nowzero	= 0;
	for (i=1; i<=m; i++)
	{
		
		for (j=1; j<=(n-1); j++)
		{
			if (fabs(detail(i,j)) < EPS) waszero++;
			if (fabs(detail(i,j)) <= (*pthreshold)) detail(i,j) = 0;
			if (fabs(detail(i,j)) < EPS) nowzero++;
			
			
		} // end for j
	} // end for i 
	

	memcpy(psmooth, smooth.get_vector(), sizeof(double) * smooth.nrow() * smooth.ncol());
	memcpy(pdetail, detail.get_vector(), sizeof(double) * detail.nrow() * detail.ncol());



	// Second, inverse ultrametric Haar wavelet transform
	// We need: cluster member array, a; and transform: smooth(,1]; and detail(.


	fvector reconstruction(m);	// m<-m Reconstructed vectors

	for (obs = 1; obs<=n; obs++)				// For all objervations
	{
		wasclusnum = obs;

		for (i=1; i<=reconstruction.length(); i++)
		{
			reconstruction[i] = 0; 
		} // end for i		

		for (j=1; j<=(n-1); j++)				// For all levels in the hierarchy
		{
			// Get cluster members, proceeding through q1, q2, ... q(n-1)
			ivector mmbrs;

			for (i=1; i<=n; i++)
			{
				if (a(i,n-j) == n+j)
				{
					mmbrs.push_back(i);
				} // end if
			} // end for i

			// Find subclusters: these are only two; and they will be found
			// in the column that follows. We do this to characterize aine, benjamin
			
			ivector nextcolmmbrs;

			nextcolmmbrs.set_size(mmbrs.length());
		
			for (i=1; i<=mmbrs.length(); i++)
			{
				nextcolmmbrs[i] = a(mmbrs[i],n-j+1);		
			}

			qsort(&nextcolmmbrs[1], nextcolmmbrs.length(), sizeof(int),  compare);	// Always, by convention, aine < benjamin
			
			aine     = nextcolmmbrs[1];
			benjamin = nextcolmmbrs[nextcolmmbrs.length()];

			clusnum = 0;
			for (k=1; k<=mmbrs.length(); k++)
			{
				if (obs == a(mmbrs[k],n))
				{
					// Have found that obs is in cluster at level j (0 <= j < n-1)
					// there can be only one occurence of this, in set mmbrs
					clusnum = j + n; // Change cluster numbering convention
					// Aine branch - In the following add or subtract detail( signal
					if (wasclusnum == aine)
					{
						for (i=1; i<=reconstruction.length(); i++)
						{
							reconstruction[i] += detail(i,2*n-clusnum);
						} // end for i
					}
					else
					if (wasclusnum == benjamin)
					{
						for (i=1; i<=reconstruction.length(); i++)
						{
							reconstruction[i] -= detail(i,2*n-clusnum);
						} // end for i
					}

					wasclusnum = clusnum;
				} // end if
			} // End of K-loop, checking through cluster memmbers
				
			mmbrs.clear();
			nextcolmmbrs.clear();

		} // end of j-loop (levels of tree)


		for (i=1; i<=reconstruction.length(); i++)
		{
			reconstruction[i] += smooth(i,1); // Add continuum (DC component)
		} // end for i		

		for (i=1; i<=reconstruction.length(); i++)
		{
			xout(obs,i) = reconstruction[i]; 
		} // end for i	

	} // end of obs-loop (set of all observetions)

	memcpy(pxout, xout.get_vector(), sizeof(double) * xout.nrow() * xout.ncol());

} // end haarum

}  // extern "C"
