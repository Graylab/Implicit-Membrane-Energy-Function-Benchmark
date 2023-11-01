# @file: make_ideal_peptide.py
# @brief: Make ideal helices given a database of sequences
# @author: Rebecca Alford (ralford3@jhu.edu)

from pyrosetta import *
from string import Template

import sys, os
import read_config
import relax_pdbstructures
#import commands
import random

from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

def main( args ):

    # Read options from the commadnline
    parser = OptionParser( usage="usage: %prog --helixdb helices.dat" )
    parser.set_description(main.__doc__)

    parser.add_option('--helixdb', '-d',
        action="store",
        help="Path to list of helix sequences",)

    parser.add_option( '--target_folder', '-t', 
                      action="store",
                      help="Name of target folder", )

    (options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options
    
    if ( not Options.helixdb ):
        print("Missing database of transmembrane helix sequences! Exiting...")
        sys.exit()
    if ( not Options.target_folder ):
        print("Missing database of transmembrane helix sequences! Exiting...")
        sys.exit()
    
    # Read path configuration file
    config = read_config.read_config()
    
    targets_path = config.benchmark_path + "targets/tilt_angle/" + Options.target_folder + "/"
    print("targets path : ", targets_path)
    
    if ( not os.path.isdir( targets_path ) ): 
        print("target folder do not exist")
        sys.exit()
    os.chdir( targets_path )
    
    print(" Reading the fasta file: ")
    print( Options.helixdb )
    with open( Options.helixdb, 'r' ) as f:
        content = f.readlines()
        print(content)
    content = [ x.strip() for x in content ]
    db = [ x.split(' ') for x in content ]

    init()

    for helix in db:
        print(helix[0])
        helixfasta = helix[0] +"/" + helix[0] + ".fasta"
        sequence_objects = rosetta.core.sequence.read_fasta_file( helixfasta )
        sequence = sequence_objects[1].sequence()
        pose = pose_from_sequence( sequence )
        for i in range(1, pose.total_residue()+1):
            pose.set_phi(i, -57.8)
            pose.set_psi(i, -47.0)

        # if ( helix[0] == "pHLIP-v14" ):
            
        # acetylate the n terminus
        # pyrosetta.rosetta.core.conformation.add_variant_type_to_conformation_residue( pose.conformation(), rosetta.core.chemical.AcetylatedNtermProteinFull, 1 )
        # ACETYLATED_NTERMINUS_CONNECTION_VARIANT
        # ACETYLATED_NTERMINUS_VARIANT
        # METHYLATED_CTERMINUS
        # N_ACETYLATION
        # N_METHYLATION
        # aminate the c terminus
        # pyrosetta.rosetta.core.conformation.add_variant_type_to_conformation_residue( pose.conformation(), rosetta.core.chemical.METHYLATED_CTERMINUS, pose.total_residue() )

        pose.dump_pdb( helix[0] + "/" + helix[0] + ".pdb" )
        
    relax_pdbstructures.run_relax_on_crystal(Options.target_folder, targetlist = Options.helixdb)   

if __name__ == "__main__" : main(sys.argv)
