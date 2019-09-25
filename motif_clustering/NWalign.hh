#ifndef NWALIGN_H_
#define NWALIGN_H_

#include <vector>

#define getmax(a,b) a>b?a:b
#define getmin(a,b) a<b?a:b

class NWalign {

public:

    inline float dist(float x[], float y[]);


    /* Dynamic Programming */
    void initialize( );
    void traceback( );
    void align( );
    float hack_TMscore(  );
    int get_first_encoded_alignment_pair( );

    void setup( float *x, float *y, int len1, int len2);
    void allocate( );
    void deallocate( );

    
    NWalign();
    ~NWalign();

private:
    const float d02;
    const float gap_open;
    int xlen, ylen;
    int Lnorm;
    float *x_ca, *y_ca;
    // matrix involved in DP
    float **score; /* Scoring table for dynamic programming (DP) */
    bool **path; /* A table to track path for DP */
    float **val; /* DP matrix */
    int * j2i;
 
};

#endif /* NWALIGN_H_ */
