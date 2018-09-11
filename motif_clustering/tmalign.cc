#include <pybind11/pybind11.h>

namespace py = pybind11;


#include <iostream>
#include <cmath>

#include "NWalign.hh"
#include "Motif.hh"

double TMscore( std::string fname1, std::string fname2 ) {

    Motif motif1 ( fname1 );
    Motif motif2 ( fname2 );

    NWalign rmsd;
    rmsd.setup( motif1.get_coords(), motif2.get_coords(), motif1.length(), motif2.length() );
    rmsd.align( );
    double v = rmsd.hack_TMscore( );
    return v;
}

PYBIND11_MODULE(tmscore, m) {
    m.doc() = "a quick method to calculate the TMscore of two motifs without alignment";

    m.def("TMscore", &TMscore, "main funciton to calculate the TMscore given two file names");
}
