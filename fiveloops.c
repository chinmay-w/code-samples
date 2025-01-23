/*
    This code was written for my class "Programming for Correctness and Performance," which
    was focused around writing pseudocode algorithms for microkernel linear-algebra operations a la
    BLAS (https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms), writing logical proofs of their correctness, and then implementing them in code. (This was for the third part).

    The specific algorithm being implemented is double-precision general matrix-matrix multiplication, aka DGEMM
    (https://www.intel.com/content/www/us/en/docs/onemkl/tutorial-c/2021-4/multiplying-matrices-using-dgemm.html).
    Intel intrinsics are used to specify assembly-level code without actually having to write assembly.

    The point was not just to achieve a "correct" algorithm, but to gradually implement various optimizations while maintaining correctness, to see how closely our code could approach the reference implementation, measured in average GFLOPs achieved by the CPU. I have added inline comments explaining the optimizations I implemented, and with this algorithm I was able to achieve 90-95% of the BLIS implementation's performance (https://github.com/flame/blis), which we were testing against.

    The code makes reference to five constants: MR, NR, MC, NC, and KC: these are the dimensions by which we divide each matrix into submatrices: one division per loop.
    Their exact values are not very relevant, since we calculated them based on our personal CPU L1/2/3 cache sizes.
*/

/*
 * Matrix multiplication requires all three matrices to have certain dimensions.
 * For the operation C = AB + C, where A B and C are all matrices:
 * Let A be an m x n matrix, B must be an n x k matrix, C must be an m x k matrix
 *
    @param m              the column-length of matrices A and C
    @param n              the row-length of matrix A, the column-length of matrix B
    @param k              the row-length of matrices B and C

    @params A, B, C       pointers to the matrices A, B and C
    @params rsA, rsB, rsC the "row stride" of matrices A, B and C
    @params csA, csB, csC the "column stride" of matrices A, B and C

    Row stride is the number of bytes in memory between separate columns of a matrix
    Column stride is the number of bytes in memory between separate rows of a matrix
*/
void fiveloops( int m, int n, int k, double *A, int rsA, int csA, 
	     double *B, int rsB, int csB,  double *C, int rsC, int csC )

{
  // fifth loop - A is passed in completely, B and C are split up into NC column-wide chunks
  for (int j=0; j<n; j+=NC) {

    // This line allows for arbitrary-size matrices to be passed through, since NC is only an upper bound
    int jb = min(NC, n-j);

    // We were provided macros for matrix element access, eg
    // #define alpha( i,j ) A[ (i)*rsA + (j)*csA ]
    // So as to not clutter our code with references to matrix striding / pointer arithmetic
    fourloops( m, jb, k, A, rsA, csA, &beta(0,j), rsB, csB, &gamma(0,j), rsC, csC );
  }
}

void fourloops( int m, int n, int k, double *A, int rsA, int csA,
       double *B, int rsB, int csB,  double *C, int rsC, int csC )

{
  // fourth loop - C is passed in completely, A is split up into KC column-wide chunks, B is buffered into Bt (temporary)

  // Aligned malloc, coupled with aligned free later, allows quicker CPU access in the kernel
  // There is of course an accompanying slight loss in space efficiency
  // '_mm_malloc' aligns with an 64-byte boundary in memory
  double *Bt = (double *) _mm_malloc(KC*NC*sizeof(double), (size_t) 64);

  for (int p=0; p<k; p+=KC) {
    int pb = min(KC, k-p);

    packB_KCxNC( pb, n, &beta(p,0), rsB, csB, Bt);

    threeloops( m, n, pb, &alpha(0,p), rsA, csA, Bt, C, rsC, csC );
  }

  _mm_free(Bt);
}

// The following two functions "pack" the submatrix Bt so that it can be accessed contiguously in memory, for increased access speed later
// A diagram of this process can be seen at https://www.cs.utexas.edu/~flame/laff/pfhp/images/Week3/BLISPicturePack.png
void packB_KCxNC( int k, int n, double *B, int rsB, int csB, double *Bt )

{
  for ( int j=0; j<n; j+= NR ){
    int jb = min( NR, n-j );
    packB_KCxNR( k, jb, &beta( 0, j ), rsB, csB, Bt );
    Bt += k * jb;
  }
}

void packB_KCxNR( int k, int n, double *B, int rsB, int csB, double *Bt )
{
  if (n == NR)
  { // in the case of a "full-size" NR panel
    for (int p=0;p<k;p++)
      for (int j=0;j<NR;j++)
        *Bt++ = beta(p,j); // put elements from B to Bt, one by one
  }
  else
  {
    for ( int p=0; p<k; p++ ) {
      for ( int j=0; j<n; j++ )
	      *Bt++ = beta( p, j );
      for ( int j=n; j<NR; j++ ) 
	      *Bt++ = 0.0; // after packing all elements, pad the rest of the panel with zeros
    }
  }
}


void threeloops( int m, int n, int k, double *A, int rsA, int csA, double *Bt,  double *C, int rsC, int csC )
{
  // third loop - Bt is passed in completely, C is split up into MC row-length chunks, A is buffered into At

  // aligned malloc, coupled with aligned free later
  double *At = (double *) _mm_malloc(MC*KC*sizeof(double), (size_t) 64);

  for (int i=0; i<m; i+=MC) {
    int ib = min(MC, m-i);

    packA_MCxKC( ib, k, &alpha(i,0), rsA, csA, At);

    twoloops( ib, n, k, At, Bt, &gamma(i,0), rsC, csC );
  }

  _mm_free(At);
}

// The same packing process is done with matrix A
void packA_MCxKC( int m, int k, double *A, int rsA, int csA, double *At )
{
  for ( int i=0; i<m; i+=MR ) {
    int ib = min(MR, m-i);
    packA_MRxKC( ib, k, &alpha(i,0), rsA, csA, At);
    At += ib * k;
  }
}


void packA_MRxKC( int m, int k, double *A, int rsA, int csA, double *At )
{
  if (m = MR) {
    for (int p=0; p<k;p++) 
      for (int i=0; i<MR; i++)
        *At++ = alpha(i,p);

  } else {
    for (int p=0; p<k;p++) 
      for (int i=0; i<m; i++)
        *At++ = alpha(i,p);
      for (int i=m; i<MR; i++)
        *At++ = 0.0;
  }
}


void twoloops( int m, int n, int k, double *At, double *Bt,  double *C, int rsC, int csC )
{
  // second loop - At is passed in completely, C and Bt are split up into NR column-wide chunks
  for (int j=0; j<n; j+=NR) {
    int jb = min(NR, n-j);

    oneloop( m, jb, k, At, &Bt[j*k], &gamma(0,j), rsC, csC );
  }
}


void oneloop( int m, int n, int k, double *At, double *Bt, double *C, int rsC, int csC )
{
  // first loop - Bt is passed in completely, At and C are split up into MR row-length chunks
  for (int i=0; i<m; i+=MR) {
    int ib = min(MR, m-i);

    dgemm_ukernel_packed(k, &At[i*k], Bt, &gamma(i,0), rsC, csC );
  }
}


 void dgemm_ukernel_packed(int k, double *mpA, double *mpB, double *C, int rsC, int csC) {

    // Finally, in the kernel, we have a guarantee that all matrices can be stored in memory

    // We allocate 4 registers for matrix C, and one for mpA and mpB (our slices of A and B) respectively
    __m256d gamma_0123_0, gamma_0123_1, gamma_0123_2, gamma_0123_3;
    __m256d alpha_0123_p, beta_p_j;

    // We load the current contents of C into our registers
    // This is what makes the math work out even when subdividing the matrices like this--C can contain a "partial result" that is used later
    // Also, even for the tiniest matrices, C = AB + C still cares about the previous contents of C
    gamma_0123_0 = _mm256_loadu_pd( &gamma(0, 0) ) ;
    gamma_0123_1 = _mm256_loadu_pd( &gamma(0, 1) ) ;
    gamma_0123_2 = _mm256_loadu_pd( &gamma(0, 2) ) ;
    gamma_0123_3 = _mm256_loadu_pd( &gamma(0, 3) ) ;


    for ( int p=0; p < k; p++){

      // aligned load is used ('load' instead of 'loadu'), since A and B are packed, for faster access
      alpha_0123_p = _mm256_load_pd( mpA ) ;
      beta_p_j     = _mm256_broadcast_sd( mpB );
      // loads values into the registers holding values for A and B
      gamma_0123_0 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_0 );
      // 'fmadd' is a batch multiply-add operation for packed vectors: gamma = alpha * beta + gamma

      beta_p_j     = _mm256_broadcast_sd( mpB+1 );
      gamma_0123_1 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_1 );

      beta_p_j     = _mm256_broadcast_sd( mpB+2 );
      gamma_0123_2 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_2 );

      beta_p_j     = _mm256_broadcast_sd( mpB+3 );
      gamma_0123_3 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_3 );
      mpA += MR;
      mpB += NR;
  }


  // stores the result to memory
  _mm256_storeu_pd( &gamma(0,0), gamma_0123_0 );
  _mm256_storeu_pd( &gamma(0,1), gamma_0123_1 );
  _mm256_storeu_pd( &gamma(0,2), gamma_0123_2 );
  _mm256_storeu_pd( &gamma(0,3), gamma_0123_3 );


  // By performing the computation in many pieces, with all matrices reduced to sizes where they can be stored directly in registers,
  // with all matrices aligned and packed, and using Intel's intrinsics to specify assembly in C,
  // I was able to achieve massive speedups while maintaining the correctness of my algorithm.
 }

