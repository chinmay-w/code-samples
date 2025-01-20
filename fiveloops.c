/*
    This code was written for my class "Programming for Correctness and Performance," which 
    was focused around writing pseudocode algorithms for microkernel linear-algebra operations a la 
    BLAS (https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms), writing logical proofs of their correctness, and then implementing them in code. (This was for the third part).

    The specific algorithm being implemented is double-precision general matrix-matrix multiplication, aka DGEMM
    (https://www.intel.com/content/www/us/en/docs/onemkl/tutorial-c/2021-4/multiplying-matrices-using-dgemm.html).

    The point was not just to achieve a "correct" algorithm, but to gradually implement various optimizations while maintaining correctness, to see how closely our code could approach the reference implementation, measured in average GFLOPs achieved by the CPU. I have added inline comments explaining each optimization I implemented, and with this algorithm I was able to achieve ~90% performance of the BLIS implementation (https://github.com/flame/blis) that we were testing against. 
*/

void fiveloops( int m, int n, int k, double *A, int rsA, int csA, 
	     double *B, int rsB, int csB,  double *C, int rsC, int csC )

{
  // fifth loop - A is passed in completely, B and C are split up into NC column-wide chunks

  for (int j=0; j<n; j+=NC) {
    int jb = min(NC, n-j);
    fourloops( m, jb, k, A, rsA, csA, &beta(0,j), rsB, csB, &gamma(0,j), rsC, csC );

  }

}

//debug only method
void print4x4matrix(double *A, int rsA, int csA, double *At) {
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      printf("%f ", alpha(i,j));
      //printf("%f ", At[j+i*4]);
    }
    printf("\n");

  }

}


void fourloops( int m, int n, int k, double *A, int rsA, int csA,
       double *B, int rsB, int csB,  double *C, int rsC, int csC )

{
  // fourth loop - c is passed in completely, A is split up into KC column-wide chunks, B is buffered into Bt
  // aligned malloc, coupled with aligned free later
  double *Bt = (double *) _mm_malloc(KC*NC*sizeof(double), (size_t) 64);
  for (int p=0; p<k; p+=KC) {
    int pb = min(KC, k-p);
    packB_KCxNC( pb, n, &beta(p,0), rsB, csB, Bt);
    #ifdef DEBUG
    printf("\n\nB:\n");
    print4x4matrix(B, rsB, csB, Bt);
    #endif

    threeloops( m, n, pb, &alpha(0,p), rsA, csA, Bt, C, rsC, csC );
  }

  _mm_free(Bt);
}


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

  if (n == NR) {
    for (int p=0;p<k;p++)
      for (int j=0;j<NR;j++)
        *Bt++ = beta(p,j);

  } else {
    for ( int p=0; p<k; p++ ) {
      for ( int j=0; j<n; j++ )
	      *Bt++ = beta( p, j );
      for ( int j=n; j<NR; j++ ) 
	      *Bt++ = 0.0;
    }
  }
}


void threeloops( int m, int n, int k, double *A, int rsA, int csA, double *Bt,  double *C, int rsC, int csC )
{

  // third loop - B is passed in completely, C is split up into MC row-length chunks, A is buffered into At
  // aligned malloc, coupled with aligned free later
  double *At = (double *) _mm_malloc(MC*KC*sizeof(double), (size_t) 64);
  for (int i=0; i<m; i+=MC) {
    int ib = min(MC, m-i);

    packA_MCxKC( ib, k, &alpha(i,0), rsA, csA, At);
    #ifdef DEBUG
    printf("\n\nA:\n");
    print4x4matrix(A, rsA, csA, At);
    #endif

    twoloops( ib, n, k, At, Bt, &gamma(i,0), rsC, csC );
  }

  _mm_free(At);
}


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
  // second loop - A is passed in completely, C and B are split up into NR column-wide chunks
  for (int j=0; j<n; j+=NR) {
    int jb = min(NR, n-j);
    oneloop( m, jb, k, At, &Bt[j*k], &gamma(0,j), rsC, csC );
  }
}


void oneloop( int m, int n, int k, double *At, double *mpB, double *C, int rsC, int csC )

{
  // first loop - B is passed in completely, A and C are split up into MR row-length chunks
  for (int i=0; i<m; i+=MR) {
    int ib = min(MR, m-i);
    dgemm_ukernel_packed(k, &At[i*k], mpB, &gamma(i,0), rsC, csC );
  }
}

// debug only method
void print16elementBuffer(double *A) {
  for (int i=0; i<16; i++) {
    printf("%f ", A[i]);
  }

  printf("\n")
}


 void dgemm_ukernel_packed(int k, double *mpA, double *mpB, double *C, int rsC, int csC) {

    #ifdef DEBUG
    printf("\n\nA:\n");
    print16elementBuffer(mpA);
    printf("\n\nB:\n");
    print16elementBuffer(mpB);
    #endif

    __m256d gamma_0123_0, gamma_0123_1, gamma_0123_2, gamma_0123_3;
    __m256d alpha_0123_p, beta_p_j;


    #ifdef DEBUG
    //printf("\nentering the ukernel_packed\n\n");
    #endif


    gamma_0123_0 = _mm256_loadu_pd( &gamma(0, 0) ) ;
    gamma_0123_1 = _mm256_loadu_pd( &gamma(0, 1) ) ;
    gamma_0123_2 = _mm256_loadu_pd( &gamma(0, 2) ) ;
    gamma_0123_3 = _mm256_loadu_pd( &gamma(0, 3) ) ;


    for ( int p=0; p < k; p++){

      //aligned load
      alpha_0123_p = _mm256_load_pd( mpA ) ;
      beta_p_j     = _mm256_broadcast_sd( mpB );
      gamma_0123_0 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_0 );

      beta_p_j     = _mm256_broadcast_sd( mpB+1 );
      gamma_0123_1 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_1 );

      beta_p_j     = _mm256_broadcast_sd( mpB+2 );
      gamma_0123_2 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_2 );

      beta_p_j     = _mm256_broadcast_sd( mpB+3 );
      gamma_0123_3 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_3 );

      mpA += MR;
      mpB += NR;
  }


  _mm256_storeu_pd( &gamma(0,0), gamma_0123_0 );
  _mm256_storeu_pd( &gamma(0,1), gamma_0123_1 );
  _mm256_storeu_pd( &gamma(0,2), gamma_0123_2 );
  _mm256_storeu_pd( &gamma(0,3), gamma_0123_3 );

 }


void dgemm_ukernel( int m, int n, int k, double *A, int rsA, int csA,
             double *B, int rsB, int csB, double *C, int rsC, int csC )

{

  #ifdef DEBUG/*
  printf("\nentering the ukernel\n\n");
  printf("m = %d, n = %d, k = %d\n", m, n, k);
  printf("rsA = %d, csA = %d, rsB = %d, csB = %d, rsC = %d, csC = %d\n", rsA, csA, rsB, csB, rsC, csC);
  printf("A = %p, B = %p, C = %p\n", A, B, C);*/

  #endif

  __m256d gamma_0123_0, gamma_0123_1, gamma_0123_2, gamma_0123_3;
  __m256d alpha_0123_p, beta_p_j;

  gamma_0123_0 = _mm256_loadu_pd( &gamma(0, 0) ) ;
  gamma_0123_1 = _mm256_loadu_pd( &gamma(0, 1) ) ;
  gamma_0123_2 = _mm256_loadu_pd( &gamma(0, 2) ) ;
  gamma_0123_3 = _mm256_loadu_pd( &gamma(0, 3) ) ;

  for ( int p=0; p < k; p++){
    alpha_0123_p = _mm256_loadu_pd( &alpha(0, p) ) ;
    
    beta_p_j     = _mm256_broadcast_sd( &beta( p, 0) );
    gamma_0123_0 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_0 );

    beta_p_j     = _mm256_broadcast_sd( &beta( p, 1) );
    gamma_0123_1 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_1 );

    beta_p_j     = _mm256_broadcast_sd( &beta( p, 2) );
    gamma_0123_2 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_2 );

    beta_p_j     = _mm256_broadcast_sd( &beta( p, 3) );
    gamma_0123_3 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_3 );
  }

  _mm256_storeu_pd( &gamma(0,0), gamma_0123_0 );
  _mm256_storeu_pd( &gamma(0,1), gamma_0123_1 );
  _mm256_storeu_pd( &gamma(0,2), gamma_0123_2 );
  _mm256_storeu_pd( &gamma(0,3), gamma_0123_3 );
}


