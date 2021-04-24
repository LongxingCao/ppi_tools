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
    if ( pdb_per_dot == 0 ) {
        pdb_per_dot = 1;
    }
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

        nlohmann::json & js = jsons[imotif];
        js["length_short"] = 0;
        js["length_long"] = 0;
        js["length"] = 0;
        if ( js.find("start") != js.end() ) {
            js["length"] = js["end"].get<int>() - js["start"].get<int>() + 1;
        }
        if ( js.find("start1") != js.end() ) {
            float val1 = js["end1"].get<int>() - js["start1"].get<int>() + 1;
            float val2 = js["end2"].get<int>() - js["start2"].get<int>() + 1;
            js["length_short"] = std::min(val1, val2);
            js["length_long"] = std::max(val1, val2);
        }


        js["length_hash"] = js["length_short"].get<int>() + js["length_long"].get<int>()*100 + js["length"].get<int>()*10000;

    }
    std::cout << std::endl;
}


void Motifs::load_motifs(  ) {

    int pdb_per_dot = motif_nums / 109;
    if ( pdb_per_dot == 0 ) {
        pdb_per_dot = 1;
    }
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


int
Motifs::find_the_best( std::vector<int> const & indices  ) {

    // Calculate this term
    // for ( int ind : indices ) {
    //     nlohmann::json & js = get_json( ind );
    //     js["ddg_per_1000_sasa"] = js["ddg_norepack"].get<float>() / js["interface_buried_sasa"].get<float>() * 1000;
    // }

    // These are the filters we care about for sorting
    std::vector<std::pair<std::string,bool>> sort_term_high_good {
        // { "interface_buried_sasa", true },
        { "ddg", false },
        //{ "hbond_to_buried_crit", true },
        //{ "hbond_to_crit", true },
        //{ "buried_crit_unsat", false },
        // { "interface_sc", true },
        // { "interface_sc_median_dist", false}
        // { "ddg_per_1000_sasa", false }
    };

    // Initialize the bests array
    std::vector<float> bests;
    for ( int i = 0; i < sort_term_high_good.size(); i++ ) {
        bool high_good = sort_term_high_good[i].second;
        if ( high_good ) {
            bests.push_back( -1000 );
        } else {
            bests.push_back( 1000 );
        }
    }

    // Find the actual bests

    for ( int ind : indices ) {
        nlohmann::json & js = get_json( ind );

        for ( int i = 0; i < sort_term_high_good.size(); i++ ) {
            std::string const & term = sort_term_high_good[i].first;
            bool high_good = sort_term_high_good[i].second;
            if ( high_good ) {
                bests[i] = std::max<float>( bests[i], js[term].get<float>() );
            } else {
                bests[i] = std::min<float>( bests[i], js[term].get<float>() );
            }
        }
    }

    std::vector<double> random_by_index( indices.size(), 0 );

    // std::vector<int> debug_sizes;
    // std::vector<float> debug_fracs;


    // Do a simple little algorithm
    // Filter at frac_of best for each of the different terms
    // If multiple pdbs exists, lower frac_of by frac_step
    // If no pdbs exists, increase frac_of by frac_step
    // If the direction changes, frac_step /= 2

    float frac_of = 1.0;
    float frac_step = 0.01;
    bool frac_down = true;  // what direction we are currently going

    int the_best = -1;
    std::vector<int> last_nonzero_contenders;

    int iter = 0;
    while ( the_best == -1 ) {

        std::vector<int> best_contenders;


        for ( int iind = 0; iind < indices.size(); iind++ ) {
            int ind = indices[iind];

            nlohmann::json & js = get_json( ind );

            bool passed_filters = true;
            for ( int i = 0; i < sort_term_high_good.size(); i++ ) {
                bool high_good = sort_term_high_good[i].second;

                float abs_frac = std::abs(bests[i]*(1-frac_of));
                float cut_value = 0;
                if ( high_good ) {
                    cut_value = bests[i] - abs_frac;
                } else {
                    cut_value = bests[i] + abs_frac;
                }
                
                std::string const & term = sort_term_high_good[i].first;

                float use_value = js[term];
                use_value += random_by_index[iind];

                if ( high_good ) {
                    passed_filters &= use_value >= cut_value;
                } else {
                    passed_filters &= use_value <= cut_value;
                }
            }

            if ( passed_filters ) {
                best_contenders.push_back( ind );
            }
        }
        // debug_sizes.push_back(best_contenders.size());
        // debug_fracs.push_back(frac_of);

        if ( best_contenders.size() == 1 ) {
            the_best = best_contenders[0];
        } else if ( best_contenders.size() == 0 ) {
            // frac_of is too big!!

            if ( ! frac_down ) frac_step /= 2;
            frac_down = true;

            frac_of -= frac_step;

        } else {
            // frac_of is too small!!

            if ( frac_down ) frac_step /= 2;
            frac_down = false;

            frac_of += frac_step;
        }

        if ( best_contenders.size() > 0 ) {
            last_nonzero_contenders = best_contenders;
        }

        iter++;
        // std::cout << frac_of << " " << best_contenders.size() << std::endl;
        // if ( iter == 1000 ) {
        //     std::cout << "Identical motifs detected, adding noise" << std::endl;

        //     for ( int i = 0; i < random_by_index.size(); i++ ) {
        //         random_by_index[i] = 0.001 * std::rand();
        //     }
        //     frac_of = 1;
        //     frac_step = 0.01;
        //     frac_down = true;
        // }

        if (iter > 10000 && last_nonzero_contenders.size() > 0) {
            std::cout << "Identical motifs?" << std::endl;
            the_best = last_nonzero_contenders[0];
            // break
            // std::cout << "bests" << std::endl;
            // for ( float best : bests ) {
            //     std::cout << "  " << best << std::endl;
            // }
            // exit(0);
        }

        if (iter > 100000 ) {
            std::cout << "Unresolvable" << frac_of<< std::endl;
            // the_best = best_contenders[0];
            // break
            std::cout << "bests" << std::endl;
            for ( float best : bests ) {
                std::cout << "  " << best << std::endl;
            }

            // for ( int i = 0; i < debug_sizes.size(); i++ ) {
            //     std::cout << "  " << debug_sizes[i] << " " << debug_fracs[i] << std::endl;
            // }
            // std::cout << "bests" << std::endl;
            // for ( float best : bests ) {
            //     std::cout << "  " << best << std::endl;
            // }
            // for ( int ind : indices ) {
            //     std::cout << "  " << get_pdb_name( ind ) << std::endl;
            // }
            exit(1);
        }
    }

    for ( int i = 0; i < indices.size(); i++ ) {
        if ( the_best == indices[i] ) return i;
    }

    std::cout << "Error" << std::endl;
    exit(0);
}



int
Motifs::write_best_info( std::ostream & out, std::vector<int> const & indices  ) {

    int local_best = find_the_best( indices );

    float average_length = 0;
    // for ( int ind : indices ) {
    //     average_length += get_json( ind )["end"].get<int>() - get_json( ind )["start"].get<int>() + 1;
    // }
    average_length /= indices.size();
    // out << "Length: " << std::setw(12) << std::left << average_length;
    out << " Sec_type: " << get_json( indices.front() )["sec_type"].get<std::string>() << " ";



    std::vector<std::string> to_write {
        // "interface_sc_median_dist",
        // "ddg_hydrophobic",
        // "hydrophobic_residue_contacts",
        // "interface_buried_sasa",
        // "ddg_norepack",
        // "interface_sc",
        // "score_per_res",
        // "vbuns5.5_heavy_ball_1.1",
        // "ddg_per_1000_sasa",
        // "fa_atr_pocket",
        // "contact_molecular_surface",
        "length_short",
        "length_long",
        "length",
        //"buried_crit_unsat",
        //"hbond_to_buried_crit",
        //"hbond_to_crit",
        "ddg"
    };


    for ( std::string const & write : to_write ) {
        float average = 0;

        for ( int ind : indices ) {
            average += get_json( ind )[write].get<float>();
        }

        average /= indices.size();

        out << write << ": " << std::setw(12) << std::left << average << " ";
    }


    nlohmann::json const & best_js = get_json( indices[local_best] );

    // out << "Best_Length: " << std::setw(12) << std::left << best_js["end"].get<int>() - best_js["start"].get<int>() + 1 << " ";


    for ( std::string const & write : to_write ) {

        out << "Best_" << write << ": " << std::setw(12) << std::left << best_js[write].get<float>() << " ";
    }


    return local_best;
}












