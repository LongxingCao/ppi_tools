#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <sstream>
#include<iomanip>


#include "defs.hh"
#include "util.hh"
#include "Motifs.hh"
#include "gzstream.h"


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
    jsons.resize( motif_nums );
    sorted_indices.resize( motif_nums );

    std::cout << "Allocating..." << std::endl;
    allocate();

    std::cout << "Loading jsons..." << std::endl;
    load_jsons();

    std::cout << "Sorting jsons..." << std::endl;
    sort_motifs();

    std::cout << "Loading motifs..." << std::endl;
    load_motifs( );

}

void Motifs::allocate( ) {
    try {
        // buffer for file reading
        file_buffer_p = new char[MAX_FILE_SIZE]; 
        // mem for coordinates
        coords_p = new float[ MAX_LENGTH * 3 * motif_nums ];
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

struct IsLessThan {
    IsLessThan( Motifs const & motifs ) : motifs_( motifs ){}
    bool operator() ( int const & il, int const & ir ) {
        nlohmann::json const & jl = motifs_.direct_access_json( il );
        nlohmann::json const & jr = motifs_.direct_access_json( ir );

        return jl["ddg"] < jr["ddg"];
    }

    Motifs const & motifs_;
};


void Motifs::sort_motifs( ) {
    for ( int i = 0; i < sorted_indices.size(); i++ ) sorted_indices[i] = i;
    std::sort( sorted_indices.begin(), sorted_indices.end(), IsLessThan( *this ) );
}

void Motifs::load_jsons(  ) {

    int pdb_per_dot = motif_nums / 109;
    for ( int imotif = 0; imotif < motif_nums; ++imotif ) {
        if ( imotif % pdb_per_dot == 0 ) std::cout << "*" << std::flush;

        const std::string & fname = motif_names[imotif];

// Now read the json file that must be in the same folder

        std::string json_fname = fname;
        if ( endswith(json_fname, ".gz") ) {
            json_fname = json_fname.substr(0, json_fname.length()-3);
        }
        if ( endswith(json_fname, ".pdb") ) {
            json_fname = json_fname.substr(0, json_fname.length()-4);
        }

        json_fname += ".json";

        std::ifstream jfile( json_fname );
        if ( ! jfile ) { 
            std::cout << "Error: failed to open file " << json_fname << " " << std::endl; 
            deallocate(); 
            exit (0); 
        }
        std::stringstream buffer;
        buffer << jfile.rdbuf();
        jsons[imotif] = nlohmann::json::parse( buffer.str() );

        jfile.close();

    }
    std::cout << std::endl;
}


void Motifs::load_motifs(  ) {

    int pdb_per_dot = motif_nums / 109;
    for ( int imotif = 0; imotif < motif_nums; ++imotif ) {
        if ( imotif % pdb_per_dot == 0 ) std::cout << "*" << std::flush;

        // now the cursor points to the start point of the coordinates
        motif_pointers[imotif] = cursor_p;

        int name_index = sorted_indices[imotif];
        const std::string & fname = motif_names[name_index];

        igzstream f;
        f.open( fname.c_str() );

        if ( ! f ) { 
            std::cout << "Error: failed to open file " << fname << " " << std::endl; 
            deallocate(); 
            exit (0); 
        }

        f.read( file_buffer_p, MAX_FILE_SIZE );
        int total_bytes = f.gcount();

        // int total_bytes = fread( file_buffer_p, 1, MAX_FILE_SIZE, f );
        // fclose(f);
        f.close();

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

    std::cout << std::endl;
}

Motifs::~Motifs() {

    deallocate();

}



void
Motifs::write_cluster_info( std::ostream & out, std::vector<int> const & indices  ) {


    float average_length = 0;
    for ( int ind : indices ) {
        average_length += get_json( ind )["end"].get<int>() - get_json( ind )["start"].get<int>() + 1;
    }
    average_length /= indices.size();
    out << "Length: " << std::setw(12) << std::left << average_length;



    std::vector<std::string> to_write {
        "interface_sc_median_dist",
        "ddg_hydrophobic",
        "hydrophobic_residue_contacts",
        "interface_buried_sasa",
        "ddg_norepack",
        "interface_sc",
        "score_per_res",
        "vbuns5.5_heavy_ball_1.1"
    };


    for ( std::string const & write : to_write ) {
        float average = 0;

        for ( int ind : indices ) {
            average += get_json( ind )[write].get<float>();
        }

        average /= indices.size();

        out << write << ": " << std::setw(12) << std::left << average;
    }


}












