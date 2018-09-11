#include <iostream>
#include<iomanip>
#include <fstream>
#include <algorithm>
#include <omp.h>

#include "NWalign.hh"
#include "Motifs.hh"

bool comp(std::string fn1, std::string fn2){
    int x, y;
    double ddg1, ddg2;

    x = fn1.find("ddg", 0);
    y = fn1.find("_", x);
    ddg1 = atof(fn1.substr(x+3, y-x-3).c_str());

    x = fn2.find("ddg", 0);
    y = fn2.find("_", x);
    ddg2 = atof(fn2.substr(x+3, y-x-3).c_str());

    return ddg1 < ddg2;

}

double calc_average_ddg( const std::vector<std::string> & pdbs ) {
    int size = pdbs.size();
    int x, y;
    double total_ddg = 0.0;
    for ( const std::string & fn : pdbs ) {
        x = fn.find("ddg", 0);
        y = fn.find("_",   0);
        total_ddg += atof( fn.substr(x+3, y-x-3).c_str() );
    }

    return total_ddg / size;
}
double calc_average_sasa( const std::vector<std::string> & pdbs ) {
    int size = pdbs.size();
    int x, y;
    double total_sasa = 0.0;
    for ( const std::string & fn : pdbs ) {
        x = fn.find("sasa", 0);
        y = fn.find("_",   0);
        total_sasa += atof( fn.substr(x+4, y-x-4).c_str() );
    }

    return total_sasa / size;
}

int main( int argc, char * argv[] )
{
    std::vector<std::string> pdbs;
    double cluster_cutoff;
    if ( argc != 3 ) {
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

    std::sort(pdbs.begin(), pdbs.end(), comp);


    const int num_pdbs = pdbs.size();
    int num_threads = omp_get_max_threads();
    Motifs motifs(pdbs);
    std::vector<NWalign> tmaligns(num_threads);

    // assign the first pdb to cluster 1.
    // intialize all the clusters
    // the cluster index starts from 0.
    std::vector<int> clusters(num_pdbs, -1);
    int cluster_nums = 0;
    for ( int ii = 0; ii < num_pdbs; ++ii ) {
        if ( -1 != clusters[ii] ) continue;

        clusters[ii] = cluster_nums;
        #pragma omp parallel for schedule (static)
        for ( int jj = ii+1; jj < num_pdbs; ++jj ) {
            if (-1 != clusters[jj]) continue;

            int thread = omp_get_thread_num();

            NWalign &tm = tmaligns[thread];

            tm.setup( motifs.get_coords(ii), motifs.get_coords(jj), motifs.length(ii), motifs.length(jj) );
            tm.align( );

            if ( tm.hack_TMscore() >= cluster_cutoff ) {
                clusters[jj] = cluster_nums;
            }
        }
        ++cluster_nums;
    }

    if ( 0 == cluster_nums ) {
        std::cout << "No clusters found!" << std::endl;
        exit(0);
    }

    // collect all the cluster infos
    std::vector<std::vector< std::string > > clustered_pdbs( cluster_nums );
    for ( int ii = 0; ii < num_pdbs; ++ii ) {
        clustered_pdbs[ clusters[ii] ].push_back( pdbs[ii] );
    }

    std::fstream f( "cluster_results.list", std::ios::out);
    if (!f.is_open()) {
        std::cout << "Error opening file cluster_results.liist" << std::endl;
            exit(0);
    }
    for ( int ii = 0; ii < clustered_pdbs.size(); ++ii ) {
        int    size     = clustered_pdbs[ii].size();
        double avg_ddg  = calc_average_ddg ( clustered_pdbs[ii] );
        double avg_sasa = calc_average_sasa( clustered_pdbs[ii] );
        f << "Cluster: " << std::setw(12) << std::left << ii + 1 << "Size: " << std::setw(12) << std::left << size << "Ddg: " << std::setw(14) << std::left << avg_ddg << "Sasa: " << std::setw(14) << std::left << avg_sasa << std::endl;
        for ( const std::string & fn : clustered_pdbs[ii] ) f << fn << " ";
        f << std::endl;
    }
    f.close();


    return 0;
}
