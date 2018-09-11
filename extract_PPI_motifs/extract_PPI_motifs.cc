// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// headers


#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/motif/reference_frames.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <devel/init.hh>
#include <numeric/conversions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <numeric/xyzVector.hh>
#include <protocols/idealize/IdealizeMover.hh>
// #include <protocols/sicdock/Assay.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/string_util.hh>

#include <boost/foreach.hpp>


// protocols
#include <protocols/simple_moves/AlignChainMover.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.fwd.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <protocols/simple_filters/InterfaceSasaFilter.hh>

#include <protocols/simple_filters/ShapeComplementarityFilter.cc>

#include <algorithm>

#include <core/kinematics/FoldTree.hh>
#include <core/pose/subpose_manipulation_util.hh>

static basic::Tracer TR( "pilot.longxing.extract_PPI_motifs" );


struct Segment {
    char sec_type;
    core::Size start;
    core::Size end;
    core::Real ddg = 0.0;
    core::Real sc = 0.0;
    core::Real sasa = 0.0;

    std::string tag() {
        using namespace ObjexxFCL::format;
        std::string tag = "ss" + std::string(1, this->sec_type) + "_" + "ddg" + F(7,2,this->ddg) + "_" + "sasa" + F(7,2,this->sasa) + "_" + "sc" + F(5,2,this->sc);
        size_t index = 0;
        while( (index = tag.find(' ',index)) != std::string::npos)
        {
           tag.erase(index,1);
        }
        return tag;
    }
};
//bool comp (Segment seg1, Segment seg2) { return ( seg1.ddg < seg2.ddg ); }

OPT_1GRP_KEY( String, longxing, ref_pdb )
OPT_1GRP_KEY( Real, longxing, ddg_threshold )
OPT_1GRP_KEY( Boolean, longxing, multi_segs)
void register_options() {
    using namespace basic::options;
    using namespace basic::options::OptionKeys;
    NEW_OPT( longxing::ref_pdb                 , "the reference pdb file"                                      , ""    );
    NEW_OPT( longxing::ddg_threshold           , "the ddg cutoff value of the motif"                           , -15.0 );
    NEW_OPT( longxing::multi_segs              , "dump multi segments or not if they are connected by a loop"  , false );
}


int main(int argc, char *argv[]) {
    try{

        using core::Real;
        using core::Size;
        
        register_options();
        devel::init(argc,argv);

        std::string fname, ref_fname, ftag;
        utility::vector1< std::string > filenames = basic::options::option[ basic::options::OptionKeys::in::file::s ]();
        if ( filenames.size() != 1 || !basic::options::option[ basic::options::OptionKeys::longxing::ref_pdb ].user() ) {
            TR.Error << "You should only give me one file each time, and also the reference pdb file!" << std::endl;
            std::exit( 1 );
        } else {
            fname = filenames[1];
						utility::vector1<std::string> temp_split = utility::string_split(fname, '/');
						std::string basename = temp_split.back();
            if   ( fname.substr(fname.size()-4,4)==".pdb" )         ftag = basename.substr(0,basename.size()-4);
            else if ( fname.substr(fname.size()-7,7)==".pdb.gz" )   ftag = basename.substr(0,basename.size()-7);
            else { TR.Error << "Unknown pdb format!" << std::endl; std::exit( 1 ); }
            ref_fname = basic::options::option[ basic::options::OptionKeys::longxing::ref_pdb ]();
            TR << "Processing pdb file " << fname << std::endl;
        }

        Real ddg_threshold = basic::options::option[ basic::options::OptionKeys::longxing::ddg_threshold ]();
				bool multi_segs    = basic::options::option[ basic::options::OptionKeys::longxing::multi_segs    ]();

        
        // load the pose
        core::pose::PoseOP pose_op     = core::import_pose::pose_from_file(fname);
        core::pose::PoseOP ref_pose_op = core::import_pose::pose_from_file(ref_fname);

        protocols::simple_moves::AlignChainMoverOP align( new protocols::simple_moves::AlignChainMover() );
        align->pose( ref_pose_op );
        align->source_chain( 2 );
        align->target_chain( 1 );
        align->apply( *pose_op );

        //pose_op->dump_pdb("test.pdb");


        // register the sasa calculator
        core::pose::metrics::CalculatorFactory & calculator_factory = core::pose::metrics::CalculatorFactory::Instance();
        if ( ! calculator_factory.check_calculator_exists("sasa") ) {
            core::pose::metrics::PoseMetricCalculatorOP sasa_calculator( new core::pose::metrics::simple_calculators::SasaCalculatorLegacy() );
            calculator_factory.register_calculator("sasa", sasa_calculator);
        }

        protocols::simple_filters::InterfaceSasaFilterOP interface_sasa( new protocols::simple_filters::InterfaceSasaFilter() );
        protocols::simple_filters::ShapeComplementarityFilterOP interface_sc( new protocols::simple_filters::ShapeComplementarityFilter() );
        //TR << "sasa: " << interface_sasa->compute(*pose_op) << std::endl;
        //TR << "interface_sc: " << interface_sc->compute(*pose_op).sc << std::endl;


        // split the chain and get the binder and target
        utility::vector1<core::pose::PoseOP> binder_target = pose_op->split_by_chain();
        core::pose::PoseOP binder_op = binder_target[1];
        core::pose::PoseOP target_op = binder_target[2];
        Size chainA_len = binder_op->size();
        TR << "length of the binder " << chainA_len << std::endl;


        // scoring function
        core::scoring::ScoreFunctionOP sf = core::scoring::get_score_function();
        core::scoring::methods::EnergyMethodOptions myopt = sf->energy_method_options();
        myopt.hbond_options().decompose_bb_hb_into_pair_energies(true);
        sf->set_energy_method_options(myopt);

        // score the pose
        Real score = sf->score(*pose_op);
        TR << "total score of the complex: " << score << std::endl;
        core::scoring::Energies    const & energies       ( pose_op->energies() );
        core::scoring::EnergyGraph const & energy_graph   ( energies.energy_graph() );
        core::scoring::EMapVector  const & energy_weights ( energies.weights() );
        utility::vector1<Real> res_ddg( chainA_len );

        // calculate the binding energy of each residue, this is just a more efficent way to calculate ddg. Why should I do this?????
        // using the code below, I can easily control what kind of interaction I want to calculate.
        for ( Size ires = 1; ires <=chainA_len; ++ires ){
            Real sum = 0;
            for ( utility::graph::Graph::EdgeListConstIter
                    iru  = energy_graph.get_node(ires)->const_upper_edge_list_begin(),
                    irue = energy_graph.get_node(ires)->const_upper_edge_list_end();
                    iru != irue; ++iru
                    ) {
                core::scoring::EnergyEdge const & edge( static_cast< core::scoring::EnergyEdge const & > (**iru) );
                Size const jres( edge.get_second_node_ind() );
                if ( jres <= chainA_len ) continue;

                sum += edge.dot( energy_weights );
            }
            
            res_ddg[ires] = sum;
        }

        // dssp
        core::scoring::dssp::Dssp dssp(*binder_op);
        std::string dssp_str = dssp.get_dssp_secstruct();
        // fix the N-ter and C-ter problems
        assert( dssp_str.length() >=2 );
        if( dssp_str[0] == 'L' && dssp_str[0] != dssp_str[1] ) dssp_str[0] = dssp_str[1];
        if( dssp_str[ dssp_str.length() -1 ] == 'L' && dssp_str[ dssp_str.length() -1 ] != dssp_str[ dssp_str.length() -2 ] ) dssp_str[ dssp_str.length() -1 ] = dssp_str[ dssp_str.length() -2 ];

        TR << "Secondary structure of the binder: " << dssp_str << std::endl;
        
        utility::vector1<Segment> segs;
        Segment seg_temp;
        for ( unsigned int ii = 0; ii < dssp_str.length(); ++ii ) {

            if ( ii == 0 ) { seg_temp.start = 1; seg_temp.sec_type = dssp_str[ii]; continue; }

            if ( dssp_str[ii] != seg_temp.sec_type ) {
                seg_temp.end = ii;
                segs.push_back( seg_temp );

                seg_temp.start = ii + 1;
                seg_temp.sec_type = dssp_str[ii];
            } else {
                //
            }

            if ( ii == dssp_str.length()-1 ) {
                seg_temp.end = ii + 1;
                segs.push_back( seg_temp );
            }
        }
        for ( utility::vector1<Segment>::iterator s=segs.begin(); s!=segs.end(); ++s ) {
            s->ddg = 0.0;
            for ( Size ires = s->start; ires <= s->end; ++ires ) s->ddg += res_ddg[ires];
        }

        utility::vector1<Size> good_motifs;
        for ( Size ii = 1; ii <= segs.size(); ++ii ) {
            if ( segs[ii].ddg < ddg_threshold ) {
                good_motifs.push_back( ii );

                core::pose::Pose motif_alone;

                core::kinematics::FoldTree ft( segs[ii].end - segs[ii].start + 1 );
                utility::vector1<Size> motif_res;
                for ( Size ires = segs[ii].start; ires <= segs[ii].end; ++ires ) { motif_res.push_back(ires); };
                core::pose::create_subpose( *binder_op, motif_res, ft, motif_alone );
                core::pose::Pose motif_complex( motif_alone );
                core::pose::append_pose_to_pose( motif_complex, *target_op, true );
                segs[ii].sc   = interface_sc->compute( motif_complex ).sc;
                segs[ii].sasa = interface_sasa->compute( motif_complex );

                std::string fout =  ftag + "_" + segs[ii].tag() + ".pdb";
                motif_alone.dump_pdb( fout );
                //motif_alone.dump_pdb("test.pdb");
            }
        }

        // dump the combination of good motifs
        if ( multi_segs &&  good_motifs.size() >= 2 ) {
            Segment comb_seg;
            comb_seg.sec_type = 'X'; // dummy sec type
            for ( uint ii = 1; ii < good_motifs.size(); ++ii ) {
                if ( good_motifs[ii] + 2 == good_motifs[ii+1] && segs[good_motifs[ii]+1].sec_type == 'L' ) {
                    good_motifs.insert( good_motifs.begin()+ii, good_motifs[ii] + 1 );
                }
            }
            utility::vector1<Size> motif_res;
            for ( uint ii : good_motifs ) {
                for ( Size ires = segs[ii].start; ires <= segs[ii].end; ++ires ) { motif_res.push_back(ires); };
            }

            comb_seg.ddg = 0.0;
            for ( Size ires : motif_res ) { comb_seg.ddg += res_ddg[ires]; }

            core::pose::Pose motif_alone;

            core::kinematics::FoldTree ft( motif_res.size() );
            core::pose::create_subpose( *binder_op, motif_res, ft, motif_alone );
            core::pose::Pose motif_complex( motif_alone );
            core::pose::append_pose_to_pose( motif_complex, *target_op, true );
            comb_seg.sc   = interface_sc->compute( motif_complex ).sc;
            comb_seg.sasa = interface_sasa->compute( motif_complex );

            std::string comb_tag = comb_seg.tag();
            std::string ss_tag;
            for ( uint ii : good_motifs ) ss_tag += std::string( 1, segs[ii].sec_type );
            comb_tag = "ss" + ss_tag + comb_tag.substr(3,comb_tag.length()-3);

            std::string fout =  ftag + "_" + comb_tag + ".pdb";
            motif_alone.dump_pdb( fout );
        }


        TR << "Job done!" << std::endl;

        //for ( utility::vector1<Segment>::const_iterator s=segs.begin(); s!=segs.end(); ++s ) {
        //    TR << s->sec_type << "  " << s->start << ", " << s->end << ", ddg: " << s->ddg << std::endl;
        //}


    } catch ( utility::excn::Exception const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        std::exit( 1 );
    }
    return 0;
}
