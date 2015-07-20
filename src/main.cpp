
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Project:				ADMM Decoder 
//                                      modif 1: Bertrand LE GAL 
//                                      modif 2: Imen DEBBABI
// Date:				20/05/2015
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//compile:
// g++ ./src/main.cpp -std=c++0x -o ./release/IEEE_ADMM_decoder
//run
// ./release/IEEE_ADMM_decoder -N .. -NmK .. -fe_limit .. -min .. -max .. -pas .. -ADMMtype .. -rho .. -mu .. -alpha1 .. -alpha2 .. -thr ..
//
//

#include "LDPC_Class.hpp"
#include <string.h>


using namespace std;

int main(int argc, char* argv[])
{
    double snrMin   = 2.00;
    double snrMax   = 2.10;
    double snrStep  = 0.5;

    unsigned int fe_limit        = 100;
    unsigned int codeword_limit  = 100000000;
    unsigned int time_limit      = 0;
    
    unsigned int num_threads     = 1;
    unsigned int maxIts          = 200;

    double mu_in  = 5.5;
    double rho_in = 1.9;
    double alpha1 = 0.6;
    double alpha2 = 0.1;
    double thr    = 1;// threshold for entropy penalized ADMM
    int N        = 576;
    int NmK      = 288;
    int ADMMtype = 2;     // LP        = 0
                          // LP1       = 1
                          // LP2       = 2
                          // LPentropy = 3

    for (int p = 1; p < argc; p++)
    {
        //
        //  SIMULATION PARAMETERS
        //
        if (strcmp(argv[p], "-N") == 0) {
            N = atoi(argv[p + 1]);
            p += 1;
        
        }else if (strcmp(argv[p], "-NmK") == 0) {
            NmK = atoi(argv[p + 1]);
            p += 1;

        }else if (strcmp(argv[p], "-fe_limit") == 0) {
            fe_limit = atoi(argv[p + 1]);
            p += 1;

        }else if (strcmp(argv[p], "-min") == 0) {
            snrMin = atof(argv[p + 1]);
            p += 1;
            
        }else if (strcmp(argv[p], "-threads") == 0) {
            num_threads = atoi(argv[p + 1]);
            p += 1;
            
        }else if (strcmp(argv[p], "-max") == 0) {
            snrMax = atof(argv[p + 1]);
            p += 1;
        
        }else if (strcmp(argv[p], "-pas") == 0) {
            snrStep = atof(argv[p + 1]);
            p += 1;
        
        }else if (strcmp(argv[p], "-ADMMtype") == 0) {
            ADMMtype = atoi(argv[p + 1]);
            p += 1;
        
        }else if (strcmp(argv[p], "-rho") == 0) {
            rho_in = atof(argv[p + 1]);
            p += 1;
        
        }else if (strcmp(argv[p], "-mu") == 0) {
            mu_in = atof(argv[p + 1]);
            p += 1;
        
        }else if (strcmp(argv[p], "-alpha1") == 0) {
            alpha1 = atof(argv[p + 1]);
            p += 1;
        
        }else if (strcmp(argv[p], "-alpha2") == 0) {
            alpha2 = atof(argv[p + 1]);
            p += 1;
        
        }else if (strcmp(argv[p], "-thr") == 0) {
            thr = atof(argv[p + 1]);
            p += 1;
        
        }else{
                exit( 0 );
        }
    }

    
    char ldpcFile[1024], codwFile[1024];
    sprintf(ldpcFile, "./mat/H_%dx%d.txt", N, NmK);
    sprintf(codwFile, "./cwd/codeword_%dx%d.txt", N, NmK);
    double rendement = (float) (N-NmK) / (float) (N);
 
    printf("(II) VersioN: (%s, %s)\n", __DATE__, __TIME__);
    printf("(II) Code LDPC (N, K, N-K): (%d, %d, %d)\n", N, N-NmK, NmK);
    printf("(II) Rendement du code    : %.3f\n", rendement);
    printf("(II) # ITERATIONs du CODE : %d\n", maxIts);
    printf("(II) FER LIMIT FOR SIMU   : %d\n", fe_limit);
    printf("(II) SIMULATION  RANGE    : [%.2f, %.2f], STEP = %.2f\n", snrMin, snrMax, snrStep);
    printf("(II) LOADED H FILENAME    : %s\n", ldpcFile);
    printf("(II) codeword file        : %s\n", codwFile);

    if (ADMMtype == 0)    
    {    
 
        printf("\n**********************************************************\n");
        printf("(II) ADMM LP decoding, rho: %.3f, mu: %.2f \n", rho_in, mu_in);
        printf("\n**********************************************************\n");
	for(double snr = snrMin; snr < snrMax; snr += snrStep)
	{
		Simulator ldpcsim(N, NmK, ldpcFile, "AWGN", snr, num_threads);
        ldpcsim.SetCodeword(codwFile);
		ldpcsim.SetCommandLineOutput(100000,   1);// create a command line output every 10000 simulations.
		int seed = rand();
		ldpcsim.SetTargets(fe_limit, codeword_limit, time_limit);// set target number of frame errors = 100
		ldpcsim.SetDecoder("ADMMLP");
		ldpcsim.SetDecoderParameters(maxIts, 1e-5, mu_in, rho_in);
		ldpcsim.SetChannelSeed(seed);
		ldpcsim.RunSim();
	}exit( 0 );
      }

    if (ADMMtype == 1)
    {
        printf("\n***************************************************************\n");
        printf("(II) ADMM L1 penalized decoding, rho: %.3f, mu: %.3f , alpha: %.3f\n", rho_in, mu_in,alpha1);
        printf("\n***************************************************************\n");
	for(double snr = snrMin; snr < snrMax; snr += snrStep)
	{
		Simulator ldpcsim(N, NmK, ldpcFile, "AWGN", snr, num_threads);
        ldpcsim.SetCodeword(codwFile);
		ldpcsim.SetCommandLineOutput(10000,   1);// create a command line output every 10000 simulations. 
		int seed = rand();
		ldpcsim.SetTargets(fe_limit, codeword_limit, time_limit);// set target number of frame errors = 100
		ldpcsim.SetDecoder("ADMML1");
		ldpcsim.SetDecoderParameters(maxIts, 1e-5, mu_in, rho_in,alpha1);
		ldpcsim.SetChannelSeed(seed);
		ldpcsim.RunSim();
	}exit( 0 );
      }
     if (ADMMtype == 2)    
    {
        printf("\n***************************************************************\n");
        printf("(II) ADMM L2 penalized decoding, rho: %.3f, mu: %.3f , alpha: %.3f\n", rho_in, mu_in,alpha2);
        printf("\n***************************************************************\n");
	for(double snr = snrMin; snr < snrMax; snr += snrStep)
	{
		Simulator ldpcsim(N, NmK, ldpcFile, "AWGN", snr, num_threads);
        ldpcsim.SetCodeword(codwFile);
		ldpcsim.SetCommandLineOutput(10000,   1);// create a command line output every 10000 simulations. 
		int seed = rand();
		ldpcsim.SetTargets(fe_limit, codeword_limit, time_limit);// set target number of frame errors = 100
		ldpcsim.SetDecoder("ADMML2");
		ldpcsim.SetDecoderParameters(maxIts, 1e-5, mu_in, rho_in,alpha2);
		ldpcsim.SetChannelSeed(seed);
		ldpcsim.RunSim();
	}exit( 0 );
      }

    if (ADMMtype == 3)
    {
        printf("\n*********************************************************************************************************\n");
        printf("(II) ADMM entropy penalized decoding, rho: %.3f, mu: %.3f , alpha: %.3f , threshold: %.3f\n", rho_in, mu_in,alpha1,thr);
        printf("\n*********************************************************************************************************\n");
        for(double snr = snrMin; snr < snrMax; snr += snrStep)
        {
            Simulator ldpcsim(N, NmK, ldpcFile, "AWGN", snr, num_threads);
            ldpcsim.SetCodeword(codwFile);
            ldpcsim.SetCommandLineOutput(10000,   1);// create a command line output every 10000 simulations.
            int seed = rand();
            ldpcsim.SetTargets(fe_limit, codeword_limit, time_limit);// set target number of frame errors = 100
            ldpcsim.SetDecoder("ADMMEntropy");
            ldpcsim.SetDecoderParameters(maxIts, 1e-5, mu_in, rho_in, alpha1,thr);
            ldpcsim.SetChannelSeed(seed);
            ldpcsim.RunSim();
        }exit( 0 );
    }

	return 0;
}
