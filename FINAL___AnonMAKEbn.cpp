#include <stdio.h>
#include <iostream>
#include <ctime>
#include <string>
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

using namespace std;

int main()
{   
	
    PFC pfc(AES_SECURITY);                                // initialise pairing-friendly curve

    time_t seed;
    time(&seed);
    irand((long)seed);

    clock_t startTA, startRSU, startV;                    // variable for measuring time
    double elapsedTA, elapsedRSU, elapsedV;


	
	//--------------------------------------------------------------------------------------
    int numberOfVehicle;
    cout << "please input number of vehicle" << endl;
    cin >> numberOfVehicle;
    //--------------------------------------------------------------------------------------
	
	
	
	// AnonMAKE starts
	// setup and initialization (by TA)
	//--------------------------------------------------------------------------------------
    Big alpha, beta, temp, q;
    G2 g1, g1p;
    G1 g2;
    GT Z;
    
    q = pfc.order();
    
    pfc.random(alpha);
    pfc.random(beta);
    pfc.random(g1);
    pfc.random(g2);
    
    startTA = clock();//---------------------------------------------time-------------------
    g1p = pfc.mult(g1, alpha);
    g1p = pfc.mult(g1p, beta);
    Z = pfc.pairing(g1, g2);
    Z = pfc.power(Z, alpha);
    elapsedTA = (clock()-startTA)/(double)CLOCKS_PER_SEC;//----------time-------------------
    
    elapsedTA = 1000.0*elapsedTA;                                   // millisecond
    printf("In AnonMAKE, TA costs %5.7lf (ms) during Setup.\n", elapsedTA);
    //--------------------------------------------------------------------------------------
    
    
    
    // vehicle registration (by TA)
    //--------------------------------------------------------------------------------------
    
	Big v[numberOfVehicle];                                     // v_i
    G1 PK_V[numberOfVehicle], BAK[numberOfVehicle];				// PK_Vi, BAK_i
    
    startTA = clock();//---------------------------------------------time-------------------
	for(int i = 0; i<numberOfVehicle; i++){
		
		// single authentication parameters
        char buffer[20] = "vehicle";
        char tempv[5];
        sprintf(tempv, "%d", i);
        strncat(buffer, tempv, 3);
        char *p = buffer;
        temp = pfc.hash_to_group(p);
        pfc.start_hash();
        pfc.add_to_hash(alpha);
        pfc.add_to_hash(temp);
        v[i] = pfc.finish_hash_to_group();
        temp = (v[i] + modmult(alpha, beta, q))%q;
		PK_V[i] = pfc.mult(g2, moddiv(alpha, temp, q));
		
		// batch authentication parameters
		BAK[i] = pfc.mult(g2, inverse(modmult(v[i], beta, q), q));
    }
    elapsedTA = (clock()-startTA)/(double)CLOCKS_PER_SEC;//----------time--------------------
    //--------------------------------------------------------------------------------------
    
    
    
    // RSU registration (by TA)
    //--------------------------------------------------------------------------------------
    Big rj;                                                    // r_j
    G1 PK_Rj; 
	G2 BVKj;                                            // PK_Rj, BVK_j
    
    startTA = clock();//---------------------------------------------time-------------------
    temp = pfc.hash_to_group((char*)"RSU");
	pfc.start_hash();
    pfc.add_to_hash(alpha);
    pfc.add_to_hash(temp);
    rj = pfc.finish_hash_to_group();
    temp = (rj + modmult(alpha, beta, q))%q;
    PK_Rj = pfc.mult(g2, moddiv(alpha, temp, q));
    temp = (alpha + moddiv(rj, beta, q))%q;
    BVKj = pfc.mult(g1, temp);
    elapsedTA += (clock()-startTA)/(double)CLOCKS_PER_SEC;//---------time--------------------
    
    elapsedTA = 1000.0*elapsedTA;                                   // millisecond
    printf("In AnonMAKE, TA costs %5.7lf (ms) during Registration.\n", elapsedTA);
    //--------------------------------------------------------------------------------------
    
    
    
    // do single authentication
    if(numberOfVehicle==1){
    	cout << "Only 1 vehicle. Performing single authentication" << endl;
        cout << endl;
        // AKE
		// partial AKE (by RSU)
        //--------------------------------------------------------------------------------------
        Big b;                                                      // b from RSU
        G1 Bj;
        
        Bj = pfc.mult(PK_Rj, b);
        
        startRSU = clock();//------------------------------------------------time--------------
        pfc.random(b);
        elapsedRSU = (clock()-startRSU)/(double)CLOCKS_PER_SEC;//------------time--------------
        //--------------------------------------------------------------------------------------
        
        
        
        // partial AKE to compute L_s1 (by vehicle)
        //--------------------------------------------------------------------------------------
        Big a;                                                     // a from vehicle
        G2 A, Ap;                                                  // A, A'
        G1 Btilde;                                                 // B tilde
        
        startV = clock();//------------------------------------------------time---------------
        pfc.random(a);
        Btilde = pfc.mult(Bj, a);
        A = pfc.mult(g1, inverse(a, q));
        Ap = pfc.mult(g1p, inverse(a, q));
        elapsedV = (clock()-startV)/(double)CLOCKS_PER_SEC;//--------------time---------------
        //--------------------------------------------------------------------------------------
        
        
        
        // partial AKE to verify L_s1 and compute L_s2 (by RSU)
        //--------------------------------------------------------------------------------------
        Big xj, xp;                                                // x_j, x'
        G2 Atilde, g1b;                                            // A tilde, g_1^b
        GT Zb;													   // Z^b
        
        Zb = pfc.power(Z, b);
        g1b = pfc.mult(g1, b);
        
        startRSU = clock();//------------------------------------------------time--------------
        if (pfc.pairing(Ap + pfc.mult(A, rj), Btilde) != Zb){
            cout << "R_j rejects single v's request" << endl;
            cout << endl;
        }

        pfc.random(xj);
        xp = moddiv(xj, b, q);
        Atilde = pfc.mult((pfc.mult(A, rj) + Ap), xp);
        elapsedRSU += (clock()-startRSU)/(double)CLOCKS_PER_SEC;//-----------time--------------
        //--------------------------------------------------------------------------------------
        
        
        
        // partial AKE to verify L_s2 and compute L_s3 (by vehicle)
        //--------------------------------------------------------------------------------------
        Big yi, mask, C1;                                          // y_i, MASK, C1
        GT Zxy;                                                    // Z^{x_{j}y_{i}}
        G2 C2;                                                     // C2
        G1 C3;                                                     // C3
        
        pfc.random(yi);
        
        startV = clock();//------------------------------------------------time---------------
        if (pfc.pairing(g1b, PK_Rj) != pfc.pairing(g1, Bj)){
        	cout << "Vehicle accpets incorrect reply of R_j. Single authentication fails." << endl;
            cout << endl;
        }
        Zxy = pfc.power(pfc.pairing(Atilde, Btilde), yi);
        mask = pfc.hash_to_aes_key(Zxy);
        C1 = lxor(mask, a);
        C2 = pfc.mult(g1, moddiv(v[0], a, q));
        C3 = pfc.mult(PK_V[0], modmult(a, yi, q));
        elapsedV += (clock()-startV)/(double)CLOCKS_PER_SEC;//-------------time---------------
        
        elapsedV = 1000.0*elapsedV;                                   // millisecond
		printf("In AnonMAKE, Vehicle costs %5.7lf (ms) during Single Authentication.\n", elapsedV);
        //--------------------------------------------------------------------------------------
        
        
        
        // partial AKE to retrieve MASK from L_s3 (by RSU)
        //--------------------------------------------------------------------------------------
        Big maskp, ap;                                             // MASK', a'
        G2 tt;													   // g1^(1/a')
        GT Zyx;                                                    // Z^{yixj}
        
        startRSU = clock();//------------------------------------------------time--------------
        Zyx = pfc.power(pfc.pairing((Ap + C2), C3), xj);
        maskp = pfc.hash_to_aes_key(Zyx);
        ap = lxor(maskp, C1);
        tt = pfc.mult(g1, inverse(ap, q));
        if(A == tt){
            cout << "AKE between RSU and vehicle is success." << endl;
            cout << "MASK: " << mask <<  endl;
            cout << "will be used for a key in this session." << endl;
        }
        cout << endl;
        elapsedRSU += (clock()-startRSU)/(double)CLOCKS_PER_SEC;//-----------time--------------
        
        elapsedRSU = 1000.0*elapsedRSU;                                   // millisecond
		printf("In AnonMAKE, RSU costs %5.7lf (ms) during Single Authentication.\n", elapsedRSU);
        //--------------------------------------------------------------------------------------
    }
    
    
    
    
    
    // do batch authentication
    else if(numberOfVehicle > 1){
    	cout << "Performing batch authentication" << endl;
        cout << endl;
    	// AKE
        // partial AKE (by RSU)
        //--------------------------------------------------------------------------------------
        Big b;														// b
		G1 Bj;                                                      // B_j
        
        startRSU = clock();//------------------------------------------------time--------------
        pfc.random(b);
        Bj = pfc.mult(PK_Rj, b);
        pfc.precomp_for_mult(Bj);									// pre-comp
        elapsedRSU = (clock()-startRSU)/(double)CLOCKS_PER_SEC;//------------time--------------
        //--------------------------------------------------------------------------------------
        
        
        
        // partial AKE (by vehicle)
        //--------------------------------------------------------------------------------------
        Big a[numberOfVehicle];                                   // a_i
        G1 Btilde[numberOfVehicle], A[numberOfVehicle];           // Btilde_i, A_i
        
        startV = clock();//------------------------------------------------time---------------
        for(int i = 0; i<numberOfVehicle; i++){
            pfc.random(a[i]);
            Btilde[i] = Bj + pfc.mult(BAK[i], a[i]);
            A[i] = pfc.mult(g2, moddiv(a[i], v[i], q));
        }
        elapsedV = (clock()-startV)/(double)CLOCKS_PER_SEC;//-------------time---------------
        
        elapsedV = 1000.0*elapsedV;                                   // millisecond
		printf("In AnonMAKE, Vehicle costs %5.7lf (ms) during Batch Authentication.\n", elapsedV);
        //--------------------------------------------------------------------------------------



        // batch authentication (by RSU)
        //--------------------------------------------------------------------------------------
        G1 righttemp, lefttemp;
        GT left, right;
        startRSU = clock();//------------------------------------------------time--------------
        for(int i = 0; i<numberOfVehicle; i++){
            lefttemp = lefttemp + Btilde[i];
            righttemp = righttemp + A[i];
        }
        left = pfc.pairing((g1p + pfc.mult(g1, rj)), lefttemp);
        right = pfc.power(Z, (b*numberOfVehicle)%q) * pfc.pairing(BVKj, righttemp);
        if (left==right){
            cout << "R_j accepts all vehicles' requests" << endl;
            cout << "Batch authentication successes." << endl;
        }else{
            cout << "R_j rejects all vehicles' requests" << endl;
            cout << "Batch authentication fails." << endl;
        }
        cout << endl;
        elapsedRSU += (clock()-startRSU)/(double)CLOCKS_PER_SEC;//-----------time--------------
        
        elapsedRSU = 1000.0*elapsedRSU;                                   // millisecond
		printf("In AnonMAKE, RSU costs %5.7lf (ms) during Batch Authentication.\n", elapsedRSU);
    }





    return 0;
}
