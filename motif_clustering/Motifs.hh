#ifndef MOTIFS_H_
#define MOTIFS_H_

#include <string>
#include <vector>

#include "json.hpp"

class Motifs {
public:
    Motifs();
    Motifs( const std::vector<std::string> & motif_fnames );
    ~Motifs();

    void allocate( );
    void deallocate( );

    void load_motifs( );
    int total_num_motifs() { return motif_nums; };
    int length( int index ) { return motif_lengths[index]; }
    float * get_coords( int index ) { return motif_pointers[index]; }

    void load_jsons();
    void sort_motifs();

    nlohmann::json & get_json( int index ) { 
        return jsons[ sorted_indices[index] ];
    }

    nlohmann::json const & get_json( int index ) const { 
        return jsons[ sorted_indices[index] ];
    }
    int length_hash( int index ) { return get_json(index)["length_hash"]; }

    nlohmann::json const & direct_access_json( int pdb_index ) const { return jsons[pdb_index]; }

    std::string get_pdb_name( int index ) {
        return motif_names[ sorted_indices[index] ];
    }


    int
    find_the_best( std::vector<int> const & indices  );

    int
    write_best_info( std::ostream & out, std::vector<int> const & indices );


private:
    std::vector< std::string > motif_names;
    std::vector< int >         motif_lengths;
    std::vector< float * >     motif_pointers;
    int motif_nums;

    float * cursor_p;
    float * coords_p;
    char * file_buffer_p;

    // The jsons are loaded in the order that they are present in
    // the original file. Then we sort the indices.
    // From then on, the indices that are reported are those 
    // of sorted_indices[ index ].
    std::vector< int > sorted_indices;
    std::vector< nlohmann::json > jsons;
};

#endif /* MOTIFS_H_ */
