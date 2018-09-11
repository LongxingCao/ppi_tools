#include <iostream>
#include <fstream>
#include <omp.h>

#include "NWalign.hh"
#include "Motif.hh"


int main( int argc, char * argv[] )
{
    std::vector<std::string> pdbs;
    if ( argc != 2 ) {
        std::cout << "Usage: " << argv[0] << " pdblist.txt" << std::endl;
    } else {
        std::fstream f( argv[1], std::ios::in );
        if (!f.is_open()) {
            std::cout << "Error opening file " << argv[1] << std::endl;
            exit(0);
        }
        std::string line;
        while( !f.eof() ) {
            std::getline(f, line);
            if ( line.length() != 0 ) {
                pdbs.push_back(line);    
            }
        }
        f.close();
    }

    const int num_pdbs = pdbs.size();
    const int total_num_scores = num_pdbs * (num_pdbs - 1) / 2;

    int num_threads = omp_get_max_threads();
    std::vector<Motif> motifs(num_threads * 2);
    std::vector<NWalign> tmaligns(num_threads);

    std::vector<uint8_t> scores(total_num_scores);

    for (int ii = 0; ii < num_pdbs; ++ii) {
        #pragma omp parallel for schedule (static)
        for ( int jj = ii + 1; jj < num_pdbs; ++jj ) {

            int thread = omp_get_thread_num();

            Motif & motif1 = motifs[thread * 2 ];
            Motif & motif2 = motifs[thread * 2 + 1];

            NWalign &tm = tmaligns[thread];

            motif1.load_motif(pdbs[ii].c_str());
            motif2.load_motif(pdbs[jj].c_str());

            tm.setup( motif1.get_coords(), motif2.get_coords(), motif1.length(), motif2.length() );
            tm.align( );
            
            int index = total_num_scores - (num_pdbs - ii) * (num_pdbs - ii - 1) / 2 + (jj - ii) - 1;

            scores[ index ] = uint8_t( (tm.hack_TMscore()) * 256 );
        }
    }

    std::fstream f( "cluster.dist.bin", std::ios::out);
    if (!f.is_open()) {
        std::cout << "Error opening file cluster.dist.bin" << std::endl;
            exit(0);
    }
    for( auto sc : scores ) f.write((char*)(&sc), 1);
    f.close();


    return 0;
}
