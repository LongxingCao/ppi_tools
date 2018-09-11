#ifndef MOTIFS_H_
#define MOTIFS_H_

#include <string>
#include <vector>

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

private:
    std::vector< std::string > motif_names;
    std::vector< int >         motif_lengths;
    std::vector< float * >     motif_pointers;
    int motif_nums;

    float * cursor_p;
    float * coords_p;
    char * file_buffer_p;
};

#endif /* MOTIFS_H_ */
