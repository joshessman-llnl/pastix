/**
 * @file compress.c
 *
 * @brief A compression example that factorizes the matrix with the Just-In-Time
 * strategy and Rank-Revealing kernels.
 *
 * @copyright 2015-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Gregoire Pichon
 * @date 2017-06-31
 *
 * @ingroup pastix_examples
 * @code
 *
 */
#include <pastix.h>
#include <spm.h>

int main (int argc, char **argv)
{
    pastix_data_t  *pastix_data = NULL; /*< Pointer to the storage structure required by pastix */
    pastix_int_t    iparm[IPARM_SIZE];  /*< Integer in/out parameters for pastix                */
    double          dparm[DPARM_SIZE];  /*< Floating in/out parameters for pastix               */
    pastix_driver_t driver;
    char           *filename;
    pastix_spm_t   *spm, *spm2;
    void           *x, *b, *x0 = NULL;
    size_t          size;
    int             check = 1;
    int             nrhs = 1;
    double          normA;

    /**
     * Initialize parameters to default values
     */
    pastixInitParam( iparm, dparm );

    /**
     * Get options from command line
     */
    pastixGetOptions( argc, argv,
                      iparm, dparm,
                      &check, &driver, &filename );

    /**
     * Read the sparse matrix with the driver
     */
    spm = malloc( sizeof( pastix_spm_t ) );
    spmReadDriver( driver, filename, spm, MPI_COMM_WORLD );
    free(filename);

    spmPrintInfo( spm, stdout );

    spm2 = spmCheckAndCorrect( spm );
    if ( spm2 != spm ) {
        spmExit( spm );
        free(spm);
        spm = spm2;
    }

    /**
     * Generate a Fake values array if needed for the numerical part
     */
    if ( spm->flttype == PastixPattern ) {
        spmGenFakeValues( spm );
    }

    /**
     * Set low-rank parameters
     */
    iparm[IPARM_COMPRESS_WHEN]       = PastixCompressWhenEnd;
    iparm[IPARM_COMPRESS_METHOD]     = PastixCompressMethodRRQR;
    iparm[IPARM_COMPRESS_MIN_WIDTH]  = 128;
    iparm[IPARM_COMPRESS_MIN_HEIGHT] = 25;
    dparm[DPARM_COMPRESS_TOLERANCE]  = 1e-8;

    /**
     * Scal the matrix to avoid unexpected rouding errors
     */
    normA = spmNorm( PastixFrobeniusNorm, spm );
    spmScalMatrix( 1./normA, spm );

    /**
     * Startup PaStiX
     */
    pastixInit( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    /**
     * Perform ordering, symbolic factorization, and analyze steps
     */
    pastix_task_analyze( pastix_data, spm );

    /**
     * Perform the numerical factorization
     */
    pastix_task_numfact( pastix_data, spm );

    /**
     * Generates the b and x vector such that A * x = b
     * Compute the norms of the initial vectors if checking purpose.
     */
    size = pastix_size_of( spm->flttype ) * spm->n;
    x = malloc( size );
    b = malloc( size );

    if ( check )
    {
        if ( check > 1 ) {
            x0 = malloc( size );
        }
        spmGenRHS( PastixRhsRndX, nrhs, spm, x0, spm->n, b, spm->n );
        memcpy( x, b, size );
    }
    else {
        spmGenRHS( PastixRhsRndB, nrhs, spm, NULL, spm->n, x, spm->n );

        /* Save b for refinement: TODO: make 2 examples w/ or w/o refinement */
        memcpy( b, x, size );
    }

    /**
     * Solve the linear system
     */
    pastix_task_solve( pastix_data, nrhs, x, spm->n );
    pastix_task_refine( pastix_data, x, nrhs, b );

    if ( check )
    {
        spmCheckAxb( nrhs, spm, x0, spm->n, b, spm->n, x, spm->n );

        if (x0) free(x0);
    }

    /**
     * Unscal the solution
     */
    spmScalVector( normA, spm, x );

    spmExit( spm );
    free( spm );
    free(b); free(x);
    pastixFinalize( &pastix_data );

    return EXIT_SUCCESS;
}

/**
 * @endcode
 */
