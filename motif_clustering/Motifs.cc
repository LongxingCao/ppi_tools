#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>


#include "defs.hh"
#include "util.hh"
#include "Motifs.hh"


Motifs::Motifs() : motif_nums(0), cursor_p(nullptr), coords_p(nullptr), file_buffer_p(nullptr) {

  // fix the logic here.
  exit(0);

}

Motifs::Motifs( const std::vector<std::string> & motif_fnames) : cursor_p(nullptr), coords_p(nullptr), file_buffer_p(nullptr) {

    motif_nums = motif_fnames.size();
    if ( 0 == motif_nums ) {
        std::cout << "No motifs loaded, why?" << std::endl;
        deallocate( );
        exit(0);
    }
    motif_names = motif_fnames;
    motif_lengths.resize( motif_nums );
    motif_pointers.resize( motif_nums );

    allocate();

    load_motifs( );

}

void Motifs::allocate( ) {
    try {
        // buffer for file reading
        file_buffer_p = new char[MAX_FILE_SIZE]; 
        // mem for coordinates
        coords_p = new float[ MAX_LENGTH * 3 * 4 * motif_nums ];
        cursor_p  = coords_p;
    }
    catch ( const std::bad_alloc & ex ) {
        std::cout << ex.what() << std::endl;
        deallocate( );
        exit(0);
    }
}

void Motifs::deallocate( ) {
    if ( nullptr != coords_p ) {
        delete [] coords_p;
    }

    if ( nullptr != file_buffer_p ) {
        delete [] file_buffer_p;
    }
}


void Motifs::load_motifs(  ) {

    for ( int imotif = 0; imotif < motif_nums; ++imotif ) {
        // now the cursor points to the start point of the coordinates
        motif_pointers[imotif] = cursor_p;

        const std::string & fname = motif_names[imotif];

        FILE * f = fopen(fname.c_str(), "r");
        if ( nullptr == f ) { 
            std::cout << "Error: failed to open file " << fname << " " << std::endl; 
            deallocate(); 
            exit (0); 
        }
        int total_bytes = fread( file_buffer_p, 1, MAX_FILE_SIZE, f );
        fclose(f);

        char * p = file_buffer_p;
        int count = 0;
        int len = 0;
        while ( count <= total_bytes ) {
            if ( *p == 'A' && *(p+1) =='T' && *(p+2) == 'O' && *(p+3) == 'M') {
                if ( *(p+13) == 'C' && *(p+14) == 'A' ) {
                    cursor_p[0] = str_to_float(p+30);
                    cursor_p[1] = str_to_float(p+38);
                    cursor_p[2] = str_to_float(p+46);
                    len += 1;
                    cursor_p += 3; // move forward by 3.
                }
                p += 81; // each ATOM line has 81 characters
                count += 81; // each ATOM line has 81 characters
            } else {
                ++p;
                ++count;
            }
        }

        motif_lengths[imotif] = len;
    }
}

Motifs::~Motifs() {

    deallocate();

}
