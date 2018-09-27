#!/usr/bin/env python

import os
import sys
import json
import argparse

from pyrosetta import *
from pyrosetta.rosetta import *

init("-corrections::beta_nov16 -holes:dalphaball /work/tlinsky/Rosetta/main/source/external/DAlpahBall/DAlphaBall.macgcc -use_truncated_termini")


script_dir = os.path.dirname(os.path.realpath(__file__))
xml = script_dir + "/extract_PPI_motifs_helper.xml"

objs = protocols.rosetta_scripts.XmlObjects.create_from_file(xml)
scorefxn = objs.get_score_function("sfxn")
remove_polars = objs.get_mover("delete_polar")

def compute_filter(pose, filter, filtername):
    print("protocols.rosetta_scripts.ParsedProtocol.REPORT: ============Begin report for " + filtername + "=======================")
    value = filter.compute(pose)
    print("============End report for " + filtername + "=======================")
    extra_value = None

    # interface_sc returns a weird type
    if ( isinstance(value, pyrosetta.rosetta.core.scoring.sc._RESULTS)):
        extra_value = value.distance
        value = value.sc

    return value, extra_value


filters_to_apply = [
"interface_buried_sasa",
"ddg_norepack",
# "ddg_hydrophobic",    // we have to do this manually because it doesn't have compute
"interface_sc",
"score_per_res",
#"buns_heavy_ball_1.1",
"vbuns5.5_heavy_ball_1.1",
"hydrophobic_residue_contacts",
"interface_vdw",
"net_interface_vdw",
"net_ddg_norepack",
]

def move_chainA_far_away(pose):
    pose = pose.clone()
    sel = core.select.residue_selector.ChainSelector("A")
    subset = sel.apply(pose)

    x_unit = numeric.xyzVector_double_t(1, 0, 0)
    far_away = numeric.xyzVector_double_t(10000, 0, 0)

    protocols.toolbox.pose_manipulation.rigid_body_move(x_unit, 0, far_away, pose, subset)

    return pose





parser = argparse.ArgumentParser()
parser.add_argument("-ref_pdb", type=str, default="", help="the reference pdb file")
parser.add_argument("-ddg_threshold", type=float, default=-20, help="the ddg cutoff value of the motif")
parser.add_argument("-sc_threshold", type=float, default=0.60, help="the interface shape complementary value of the motif")
parser.add_argument("-sasa_threshold", type=float, default=750, help="the interface buried sasa cutoff value of the motif")
parser.add_argument("-multi_segs", type=bool, default=False, help="dump multi segments or not if they are connected by a loop")
parser.add_argument("-out_prefix", type=str, default="", help="prefix on out files")
parser.add_argument("-s", type=str, nargs="*", help="name of the input pdb file")
parser.add_argument("-l", type=str, help="a list file of all pdb files")
args = parser.parse_args(sys.argv[1:])

if args.s == None:
    assert (args.l != None)
    with open(args.l) as f:
        fnames = [line.strip() for line in f]
else:
    fnames = args.s

# global settings
ref_fname = args.ref_pdb
ref_pose  = pose_from_file(ref_fname)
align = protocols.simple_moves.AlignChainMover()
align.pose( ref_pose )
align.source_chain( 2)
align.target_chain( 1)

sasa_threshold = args.sasa_threshold
sc_threshold = args.sc_threshold
ddg_threshold = args.ddg_threshold
multi_segs = args.multi_segs

for fname in fnames:
    try:
        basename = os.path.basename(fname)

        ftag = ""
        if (fname.endswith(".pdb")):
            ftag = basename[:-4]
        if (fname.endswith(".pdb.gz")):
            ftag = basename[:-7]
        assert(ftag)
        ftag = args.out_prefix + ftag

        print("Processing pdb file %s"%fname )

        pose = pose_from_file(fname)
        align.apply(pose)

        binder_target = pose.split_by_chain()
        binder = binder_target[1]
        target = binder_target[2]
        chainA_len = binder.size()
        print("length of the binder %i"%chainA_len)


        score = scorefxn(pose)
        print("total score of the complex: %.3f"%score)


        separate_pose = move_chainA_far_away(pose)
        scorefxn(separate_pose)


        per_res_ddg = utility.vector1_double()
        # per_res_score = utility.vector1_double()
        # per_res_bound = utility.vector1_double()
        for i in range(1, chainA_len+1):
            ddg = 2* ( pose.energies().residue_total_energy(i) - separate_pose.energies().residue_total_energy(i) )
            per_res_ddg.append(ddg)
            # per_res_score.append(separate_pose.energies().residue_total_energy(i))
            # per_res_bound.append(pose.energies().residue_total_energy(i))

        per_res_ddg0 = list(per_res_ddg)

############################################

        all_positions = utility.vector1_unsigned_long()
        for i in range(1, binder.size()+1):
            all_positions.append(i)

        poly_ala_binder_pose = pose.clone()
        protocols.toolbox.pose_manipulation.construct_poly_ala_pose(poly_ala_binder_pose, all_positions, True, True, True)
        scorefxn(poly_ala_binder_pose)

        poly_ala_binder_separate_pose = move_chainA_far_away(poly_ala_binder_pose)
        scorefxn(poly_ala_binder_separate_pose)

        per_res_ddg_ala = utility.vector1_double()
        # per_res_score = utility.vector1_double()
        # per_res_bound = utility.vector1_double()
        for i in range(1, chainA_len+1):
            ddg = 2* ( poly_ala_binder_pose.energies().residue_total_energy(i) - poly_ala_binder_separate_pose.energies().residue_total_energy(i) )
            per_res_ddg_ala.append(ddg)
            # per_res_score.append(separate_pose.energies().residue_total_energy(i))
            # per_res_bound.append(pose.energies().residue_total_energy(i))

        per_res_ddg0_ala = list(per_res_ddg_ala)


######################################
        dssp = core.scoring.dssp.Dssp( binder )
        dssp_str = dssp.get_dssp_secstruct()

        if ( dssp_str[0] == "L" and dssp_str[1] != "L"):
            tmp = list(dssp_str)
            tmp[0] = tmp[1]
            dssp_str = ''.join(tmp) 
        if (dssp_str[-1] == "L" and dssp_str[-2] != "L"):
            tmp = list(dssp_str)
            tmp[-1] = tmp[-2]
            dssp_str = ''.join(tmp) 

        segs = []
        seg_temp = {}

        for ii in range(len(dssp_str)):
            if (ii == 0):
                seg_temp["start"] = 1
                seg_temp["sec_type"] = dssp_str[ii]
                continue

            if ( dssp_str[ii] != seg_temp["sec_type"]):
                seg_temp["end"] = ii
                segs.append(seg_temp)
                seg_temp = {}
                seg_temp["start"] = ii + 1
                seg_temp["sec_type"] = dssp_str[ii]

            if ( ii == len(dssp_str)-1 ):
                seg_temp["end"] = ii + 1
                segs.append(seg_temp)

        for seg in segs:
            seg["len"] = seg["end"] - seg["start"] + 1
            seg["sequence"] = binder.sequence()[seg["start"]-1:seg["end"]+1-1]


        for seg in segs:
            seg["ddg"] = 0
            for ires in range(seg["start"], seg["end"]+1):    
                seg["ddg"] += per_res_ddg[ires]


        for seg in segs:
            if ( seg["ddg"] > ddg_threshold ):
                continue


            ft = core.kinematics.FoldTree( seg["end"] - seg["start"] + 1 )
            motif_res = utility.vector1_unsigned_long()
            for i in range(seg["start"], seg["end"] + 1 ):
                motif_res.append(i)

            motif_alone = core.pose.Pose()
            # core.pose.create_subpose( binder, motif_res, ft, motif_alone)
            core.pose.pdbslice(motif_alone, binder, motif_res)

            motif_complex = motif_alone.clone()
            motif_complex.append_pose_by_jump(target, 1)
            # core.pose.append_pose_to_pose( motif_complex, target, True)
            # motif_complex.dump_pdb("asdf.pdb")
            scorefxn(motif_complex)


            for filt in filters_to_apply:
                the_filter = objs.get_filter(filt)

                # Get rid of stochastic filter
                if ( isinstance(the_filter, pyrosetta.rosetta.protocols.filters.StochasticFilter) ):
                    the_filter = the_filter.subfilter()

                result, extra = compute_filter(motif_complex, the_filter, filt)
                seg[filt] = result
                if (not extra is None):
                    seg["interface_sc_median_dist"] = extra

            if seg['interface_sc'] < sc_threshold : continue
            if seg['interface_buried_sasa'] < sasa_threshold : continue

            ddg_hydrophobic_pre = objs.get_filter("ddg_hydrophobic_pre").subfilter()
            tmp = pose.clone()
            remove_polars.apply(tmp)
            seg["ddg_hydrophobic"] = ddg_hydrophobic_pre.compute(tmp)

            hotspots_energy = []
            for ii in range(seg["start"], seg["end"]+1):
                hotspots_energy.append(per_res_ddg0[ii-1] - per_res_ddg0_ala[ii-1])
            seg["hotspots"] = hotspots_energy

            fout = ftag + "_%i_%i_%s.pdb.gz"%(seg["start"], seg["end"], seg["sec_type"])
            motif_alone.dump_pdb(fout)

            fname = fout[:-7]
            seg["name"] = fname
            f = open(fname + ".json", "w")
            f.write(json.dumps(seg, indent=4))
            f.close()


        print("Job done!")

    except Exception as e:
        print("Error!!!")
        print(e)



