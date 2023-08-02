#include <stdio.h>
#include <iostream>
#include <ctime>
#include <climits>
//********* choose just one of these pairs **********
//#define MR_PAIRING_CP      // AES-80 security   
//#define AES_SECURITY 80

//#define MR_PAIRING_MNT	// AES-80 security
//#define AES_SECURITY 80

#define MR_PAIRING_BN    // AES-128 or AES-192 security
#define AES_SECURITY 128
//#define AES_SECURITY 192

//#define MR_PAIRING_KSS    // AES-192 security
//#define AES_SECURITY 192

//#define MR_PAIRING_BLS    // AES-256 security
//#define AES_SECURITY 256
//*********************************************

#include "pairing_3.h"

#define MIN_TIME 10.0
#define MIN_ITERS 20

using namespace std;

int main()
{   
	
    PFC pfc(AES_SECURITY);  // initialise pairing-friendly curve
    
    Big q;
    q = pfc.order();

    time_t seed;
    time(&seed);
    irand((long)seed);

    clock_t start;                                        // variable for measuring time
    double elapsed;

    int iterations;
    cout << "Individual operation testbed on BN curve." << endl;
    cout << endl;
    
    
    
    
    // Bilinear pairing
    G1 g2;
    G2 g1;
    GT Z;
    
    pfc.random(g1);
	pfc.random(g2);
    
    iterations = 0;
    start=clock();
    do {
       Z = pfc.pairing(g1, g2);
       iterations++;
       elapsed=(clock()-start)/(double)CLOCKS_PER_SEC;
    } while (elapsed<MIN_TIME || iterations<MIN_ITERS);
    elapsed=1000.0*elapsed/iterations;
    printf("After %d times iteration, each Bilinear Pairing costs %5.7lf (ms) in average.\n", iterations, elapsed);
    cout << endl;
    
    
    
    
    
    
    // Operation in G1, G2, GT
	Big alpha, beta, gamma;
	G2 g1p, temp1;
	G1 g2p, temp2;
	GT Zp;
	
	pfc.random(alpha);
	pfc.random(beta);
	pfc.random(gamma);
    
    iterations = 0;
    start=clock();
    do {
       g1p = pfc.mult(g1, alpha);
       iterations++;
       elapsed=(clock()-start)/(double)CLOCKS_PER_SEC;
    } while (elapsed<MIN_TIME || iterations<MIN_ITERS);
    elapsed=1000.0*elapsed/iterations;
    printf("After %d times iteration, each scale multiplication on G1 costs %5.7lf (ms) in average.\n", iterations, elapsed);
    cout << endl;
    
    iterations = 0;
    start=clock();
    do {
       g2p = pfc.mult(g2, beta);
       iterations++;
       elapsed=(clock()-start)/(double)CLOCKS_PER_SEC;
    } while (elapsed<MIN_TIME || iterations<MIN_ITERS);
    elapsed=1000.0*elapsed/iterations;
    printf("After %d times iteration, each scale multiplication on G2 costs %5.7lf (ms) in average.\n", iterations, elapsed);
    cout << endl;
    
    iterations = 0;
    start=clock();
    do {
       temp1 = g1+g1p;
       iterations++;
       elapsed=(clock()-start)/(double)CLOCKS_PER_SEC;
    } while (elapsed<MIN_TIME || iterations<MIN_ITERS);
    elapsed=1000.0*elapsed/iterations;
    printf("After %d times iteration, each point addition on G1 costs %5.7lf (ms) in average.\n", iterations, elapsed);
    cout << endl;
    
    iterations = 0;
    start=clock();
    do {
       temp2 = g2+g2p;
       iterations++;
       elapsed=(clock()-start)/(double)CLOCKS_PER_SEC;
    } while (elapsed<MIN_TIME || iterations<MIN_ITERS);
    elapsed=1000.0*elapsed/iterations;
    printf("After %d times iteration, each point addition on G2 costs %5.7lf (ms) in average.\n", iterations, elapsed);
    cout << endl;
    
    iterations = 0;
    start=clock();
    do {
       Zp = pfc.power(Z, gamma);
       iterations++;
       elapsed=(clock()-start)/(double)CLOCKS_PER_SEC;
    } while (elapsed<MIN_TIME || iterations<MIN_ITERS);
    elapsed=1000.0*elapsed/iterations;
    printf("After %d times iteration, each scale multiplication on GT costs %5.7lf (ms) in average.\n", iterations, elapsed);
    cout << endl;
    
    
	
	
	
	// Modular inverse
	Big a;
	G1 g3, g4;
	
	pfc.random(a);
    pfc.random(g3);
	
	iterations = 0;
    start=clock();
    do {
       g4 = pfc.mult(g3, inverse(a,q));
       iterations++;
       elapsed=(clock()-start)/(double)CLOCKS_PER_SEC;
    } while (elapsed<MIN_TIME || iterations<MIN_ITERS);
    elapsed=1000.0*elapsed/iterations;
    printf("After %d times iteration, each Modular Inverse costs %5.7lf (ms) in average.\n", iterations, elapsed);
    cout << endl;
	
	
	
	
	
	// Hash and xor 
	Big hashOuput;
	
	iterations = 0;
    start=clock();
    do {
	   hashOuput = pfc.hash_to_aes_key(Zp);
       iterations++;
       elapsed=(clock()-start)/(double)CLOCKS_PER_SEC;
    } while (elapsed<MIN_TIME || iterations<MIN_ITERS);
    elapsed=1000.0*elapsed/iterations;
	printf("After %d times iteration, each Hash (to AES key) costs %5.7lf (ms) in average.\n", iterations, elapsed);
    cout << endl;
    
    
    
    Big A, B, C;
    
    pfc.random(A);
	pfc.random(B);
    
	do {
       C = lxor(A, B);
       iterations++;
       elapsed=(clock()-start)/(double)CLOCKS_PER_SEC;
    } while (elapsed<MIN_TIME || iterations<MIN_ITERS);
    elapsed=1000.0*elapsed/iterations;
    printf("After %d times iteration, each Xor costs %5.7lf (ms) in average.\n", iterations, elapsed);
	cout << endl;
	
	
	/* 
	cout << "The char size of elements in G1: " << sizeof(G2) << endl;
	cout << "The char size of elements in G2: " << sizeof(G1) << endl;
	cout << "The char size of elements in GT: " << sizeof(GT) << endl;
	cout << "The char size of type big: " << sizeof(big) << endl;
	cout << "The char size of hash (to AES key) output: " << sizeof(hashOuput) << endl;
	cout << "CHAR_BIT: " << CHAR_BIT << endl;
	*/
	
	cout << "Order q: " << q << endl;
	cout << "Big: " << A << endl;
	cout << "Hash output: " << hashOuput << endl;
	cout << "G1: " << g1.g << endl;
	cout << "G2: " << g2.g << endl;
	cout << "GT: " << Z.g << endl;
	
	// measure bytes
	pfc.precomp_for_mult(g1p);
	pfc.precomp_for_mult(g2p);
	pfc.precomp_for_power(Zp);
	
	char *bytes1;
	char *bytes2;
	char *bytesT;
	
	int leng1 = g1p.spill(bytes1);
	int leng2 = g2p.spill(bytes2);
	int lengT = Zp.spill(bytesT);
	
	cout << "Bytes of elements in G1: " << leng1 << endl;
	cout << "Bytes of elements in G2: " << leng2 << endl;
	cout << "Bytes of elements in GT: " << lengT << endl;
	
	g1p.restore(bytes1);
	g2p.restore(bytes2);
	Zp.restore(bytesT);







    return 0;
}
