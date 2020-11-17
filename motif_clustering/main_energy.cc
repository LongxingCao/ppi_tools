#include <iostream>
#include<iomanip>
#include <fstream>
#include <algorithm>
#include <omp.h>
#include <sstream>

#include "NWalign.hh"
#include "Motifs.hh"

int main( int argc, char * argv[] )
{
    std::vector<std::string> pdbs;
    double cluster_cutoff;
    if ( argc < 3 ) {
        std::cout << "Usage: " << argv[0] << " pdblist.txt cluster_cutoff" << std::endl;
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

        cluster_cutoff = atof( argv[2] );
    }

    bool same_length = argc >= 4;
    if ( same_length ) {
        std::cout << "Only clustering pdbs of same length!!!" << std::endl;
    }

    // This sorts the pdbs!!!!!
    Motifs motifs(pdbs);

    // Don't use pdbs again in this file!!!!
    // Only use motifs.get_pdb_name()

    std::cout << "Clustering..." << std::endl;

    const int num_pdbs = pdbs.size();
    int num_threads = omp_get_max_threads();
    
    std::vector<NWalign> tmaligns(num_threads);

    int pdb_per_dot = num_pdbs / 109;
    if ( pdb_per_dot == 0 ) {
        pdb_per_dot = 1;
    }

    // assign the first pdb to cluster 1.
    // intialize all the clusters
    // the cluster index starts from 0.
    std::vector<int> clusters(num_pdbs, -1);
    std::vector<int> alignments(num_pdbs, -1);
    int cluster_nums = 0;
    for ( int ii = 0; ii < num_pdbs; ++ii ) {
        if ( ii % pdb_per_dot == 0 ) std::cout << "*" << std::flush;
        if ( -1 != clusters[ii] ) continue;

        clusters[ii] = cluster_nums;
        alignments[ii] = 0;
        #pragma omp parallel for schedule (static)
        for ( int jj = ii+1; jj < num_pdbs; ++jj ) {
            if (-1 != clusters[jj]) continue;

            if ( same_length && motifs.length_hash(ii) != motifs.length_hash(jj) ) continue;

            int thread = omp_get_thread_num();

            NWalign &tm = tmaligns[thread];

            tm.setup( motifs.get_coords(ii), motifs.get_coords(jj), motifs.length(ii), motifs.length(jj) );
            tm.align( );


            if ( tm.hack_TMscore() >= cluster_cutoff ) {
                clusters[jj] = cluster_nums;
                alignments[jj] = tm.get_first_encoded_alignment_pair( );
            }
        }
        ++cluster_nums;
    }

    std::cout << std::endl;

    if ( 0 == cluster_nums ) {
        std::cout << "No clusters found!" << std::endl;
        exit(0);
    }

    // collect all the cluster infos
    std::vector<std::vector< std::string > > clustered_pdbs( cluster_nums );
    std::vector<std::vector< int > > clustered_indices( cluster_nums );
    std::vector<std::vector< int> > clustered_alignments( cluster_nums );
    for ( int ii = 0; ii < num_pdbs; ++ii ) {
        clustered_pdbs[ clusters[ii] ].push_back( motifs.get_pdb_name(ii) );
        clustered_indices[ clusters[ii] ].push_back( ii );
        clustered_alignments[ clusters[ii] ].push_back( alignments[ii] );
    }

    std::vector<std::string> all_strings(clustered_pdbs.size());

    int iters = clustered_pdbs.size();
    #pragma omp parallel for schedule (dynamic, 1)
    for ( int ii = 0; ii < iters; ++ii ) {

        std::stringstream f;

        int    size     = clustered_pdbs[ii].size();
        f << "Cluster: " << std::setw(12) << std::left << ii + 1 << "Size: " << std::setw(12) << std::left << size;

        int local_best = motifs.write_best_info( f, clustered_indices[ii] );
        f << std::endl;

        // dump the alignment info
        for( int local = 0; local < clustered_pdbs[ii].size(); local++ ) {
            f << clustered_alignments[ii][local] << " ";
        }
        f << std::endl;

        //f << clustered_pdbs[ii][local_best] << " ";
        for ( int local = 0; local < clustered_pdbs[ii].size(); local++ ) {
            //if (local == local_best) continue;

            const std::string & fn = clustered_pdbs[ii][local];

            f << fn << " ";
        }

        f << std::endl;

        all_strings[ii] = f.str();

    }

    std::fstream f( "cluster_results.list", std::ios::out);
    if (!f.is_open()) {
        std::cout << "Error opening file cluster_results.liist" << std::endl;
            exit(0);
    }

    for ( std::string const & str : all_strings ) {
        f << str;
    }

    f.close();


    return 0;
}
