
#include "defs.hh"
#include "NWalign.hh"


inline float NWalign::dist(float x[], float y[]) {

    float dx = x[0] - y[0];
    float dy = x[1] - y[1];
    float dz = x[2] - y[2];

    return (dx * dx + dy * dy + dz * dz);

}


void NWalign::initialize( ) {

    for (int ii = 0; ii <= xlen; ++ii) {
        val[ii][0] = 0;
        path[ii][0] = false;
    }
    for (int jj = 0; jj <= ylen; ++jj) {
        val[0][jj] = 0;
        path[0][jj] = false;
    }
    for (int jj = 0; jj < ylen; ++jj) {
        j2i[jj] = -1;
    }
}

void NWalign::traceback(  ) {

    int i = xlen;
    int j = ylen;
    float h, v;
    while (i > 0 && j > 0) {
        if (path[i][j]) //from diagonal
        {
            j2i[j - 1] = i - 1;
            i--;
            j--;
        } else {
            h = val[i - 1][j];
            if (path[i - 1][j])
                h += gap_open;

            v = val[i][j - 1];
            if (path[i][j - 1])
                v += gap_open;

            if (v >= h)
                j--;
            else
                i--;
        }
    }
}

void NWalign::align(  ) {

    int i, j;
    float h, v, d;

    //initialization
    initialize( );

    float dij;

    //decide matrix and path
    for (i = 1; i <= xlen; i++) {
        for (j = 1; j <= ylen; j++) {
            //d=val[i-1][j-1]+score[i][j]; //diagonal

            dij = dist( x_ca + (i - 1) * 3, y_ca + (j - 1) * 3 );
            d = val[i - 1][j - 1] + d02 / (d02 + dij);

            //symbol insertion in horizontal (= a gap in vertical)
            h = val[i - 1][j];
            if (path[i - 1][j]) //aligned in last position
                h += gap_open;

            //symbol insertion in vertical
            v = val[i][j - 1];
            if (path[i][j - 1]) //aligned in last position
                v += gap_open;

            if (d >= h && d >= v) {
                path[i][j] = true; //from diagonal
                val[i][j] = d;
            } else {
                path[i][j] = false; //from horizontal
                if (v >= h)
                    val[i][j] = v;
                else
                    val[i][j] = h;
            }
        } //for i
    } //for j

    //trace back to extract the alignment
    traceback( );

}

float NWalign::hack_TMscore(  ) {

    float TMscore = 0.0;
    float d;

    for (int ii = 0; ii < ylen; ++ii) {
        if ( -1 != j2i[ii] ) {
            d = dist(x_ca + j2i[ii] * 3, y_ca + ii * 3 );
            TMscore += d02 / (d02 + d);
        } 
    }
    return TMscore / Lnorm;
}

// this function return the first alignment pair, indexed from 0 and encoded as 1000*i+j
// it is very easy to decode
// i = x / 1000
// j = x % 1000
int NWalign::get_first_encoded_alignment_pair( ) {
    for ( int ii = 0; ii < ylen; ++ii ) {
        if ( -1 != j2i[ii] ) {
            return j2i[ii] * 1000 + ii;
        }
    }
    return -1;
}

void NWalign::setup( float *x, float *y, int len1, int len2 ) {

    Lnorm = getmin( len1, len2 );
    x_ca = x, y_ca = y;
    xlen = len1, ylen = len2;
    
}


void NWalign::allocate( ) {

    score = new float*[ MAX_LENGTH + 1 ];
    path  = new bool*  [ MAX_LENGTH + 1 ];
    val   = new float*[ MAX_LENGTH + 1 ];
    for ( int ii = 0; ii <= MAX_LENGTH; ++ii ) {
        score[ ii ] = new float[ MAX_LENGTH + 1 ];
        path [ ii ] = new bool  [ MAX_LENGTH + 1 ];
        val  [ ii ] = new float[ MAX_LENGTH + 1 ];
    }

    j2i = new int[MAX_LENGTH];

}

void NWalign::deallocate( ) {

    if ( nullptr != score ) {
        for (int ii = 0; ii <= MAX_LENGTH; ++ii) {
            delete [] score[ ii ];
        }
        delete [] score;
        score = nullptr;
    }

    if ( nullptr != path ) {
        for (int ii = 0; ii <= MAX_LENGTH; ++ii) {
        delete [] path [ ii ];
        }
        delete [] path ;
        path = nullptr;
    }

    if ( nullptr != val ) {
        for (int ii = 0; ii <= MAX_LENGTH; ++ii) {
        delete [] val  [ ii ];
        }
        delete [] val  ;
        val = nullptr;
    }

    if ( nullptr !=j2i ) { delete [] j2i; j2i = nullptr; }

}

NWalign::NWalign() : d02(4.0 /* from TMalign */), gap_open(-0.6 /* from TMalign */ ), xlen(0), ylen(0), Lnorm(0),
            x_ca(nullptr), y_ca(nullptr), score(nullptr), path(nullptr), val(nullptr), j2i(nullptr) {

    allocate();

}

NWalign::~NWalign() {

    /* empty */
    deallocate();

}
