# @file: make_ideal_peptide.py
# @brief: Make ideal helices given a database of sequences
# @author: Rebecca Alford (ralford3@jhu.edu)

from pyrosetta import *
from string import Template

import sys, os
#import commands
import random

from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

from pyrosetta.rosetta.utility import vector1_bool
from pyrosetta.rosetta.core.chemical import aa_from_oneletter_code
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.pose import PDBInfo
from pyrosetta.rosetta.core.chemical import VariantType
from pyrosetta.rosetta.core.pack.task import TaskFactory
import pyrosetta.rosetta.core.select.residue_selector as residue_selector

def main( args ):

    # Read options from the commadnline
    parser = OptionParser( usage="usage: %prog --helixdb helices.dat" )
    parser.set_description(main.__doc__)

    parser.add_option('--helixdb', '-d',
        action="store",
        help="Path to list of helix sequences",)

    (options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options

    if ( not Options.helixdb ):
        print("Missing database of transmembrane helix sequences! Exiting...")
        sys.exit()
    print(" Reading the fasta file: ")
    print( Options.helixdb )
    with open( Options.helixdb, 'r' ) as f:
        content = f.readlines()
        print(content)
    content = [ x.strip() for x in content ]
    db = [ x.split(' ') for x in content ]

    init()

    for helix in db:
        helixfasta = helix[0]+"/"+helix[0] + ".fasta"
        sequence_objects = rosetta.core.sequence.read_fasta_file( helixfasta )
        sequence = sequence_objects[1].sequence()
        pose = pose_from_sequence( sequence )
        for i in range(1, pose.total_residue()+1):
            pose.set_phi(i, -57.8)
            pose.set_psi(i, -47.0)

        repack_radius = 6
        sfxn=create_score_function('franklin2021')
        fa_water_bilayer_wts = 0.863
        fa_imm_elec_wts = 0.001
        f_elec_lipidlayer_wts = 0.152
        sfxn.set_weight(fa_water_to_bilayer, fa_water_bilayer_wts)
        sfxn.set_weight(fa_imm_elec, fa_imm_elec_wts)
        sfxn.set_weight(f_elec_lipidlayer, f_elec_lipidlayer_wts)
        sfxn.set_weight(fa_elec, 1.00)
        sfxn.set_weight(menv_pH, 0.00)
        
        #TO OD remove the charged end points:
        # ref: vikrams paper
        #c-terminal amidation:
        #ModifyVariantType name="vartype" add_type="CTERM_AMIDATION" residue_selector="select_cterm"
        #Nterminal acetylation or methylation
        # ModifyVariantType name="vartype" add_type="ACETYLATED_NTERMINUS_VARIANT/METHYLATED_NTERM_VARIANT" remove_type="LOWER_TERMINUS_VARIANT" residue_selector="select_nterm"
        
        pose.dump_pdb( helix[0] + ".pdb" )

if __name__ == "__main__" : main(sys.argv)
