#!/usr/bin/env python
""" Generate abinitio structures given the fasta file

This script runs Rosetta abinitio for a given sequence. 
This script is trying to generate a disordered state. 

Authors: 
  Rituparna Samanta <rsamant2@jhu.edu> 

Example: 
  $ import 

Arguments: 
  - config: Container with path to benchmark and rosetta files
  - energy_fxn: Weights file for energy function of interest
  - targets: Location of design targets
  - test_name: Name of test
  - restore: Use talaris2014 energy function parameters

Requirements: 
  - Rosetta release 246 or greater
"""

import random, os, sys
import hpc_util, read_config
from string import Template
import numpy as np

from pyrosetta import *
from pyrosetta.teaching import *
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.select.movemap import *
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task import operation
from pyrosetta.rosetta.core.fragment import *

from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

def run_abinitio_structure_prediction( sfxn, peptide_name, benchmark_path ): 
   
    target_folder = benchmark_path + "targets/stability/C5_pHLIP_helical_peptides"
    with_mers = False
    #fasta_file = target_folder + "/" + peptide_name + "/" + peptide_name + ".fasta"
    disordered_folder = "/home/rsamant2/Softwares/folding_sequence/disordered_states/"
   
    if(~with_mers):
      fasta_file = target_folder + "/" + peptide_name + "/" + peptide_name + ".fasta"
    else:
      fasta_file = disordered_folder + "/" + peptide_name + "/" + peptide_name + ".fasta"
      file_3mer = disordered_folder + "/" + peptide_name + "/" + peptide_name + ".robetta.3mer"
      file_9mer = disordered_folder + "/" + peptide_name + "/" + peptide_name + ".robetta.9mer"
    
    f = open(fasta_file,'r')
    #skips the header
    lines = f.readlines()[1:]
    print(lines[0])
    #closes the file
    f.close()
    
    span_file = target_folder + "/" + peptide_name + "/" + peptide_name + ".span"
    #just initializes to a stretched peptide
    start_pose = pose_from_sequence(lines[0])
    # Add membrane to pose
    add_memb = rosetta.protocols.membrane.AddMembraneMover( span_file )
    add_memb.apply( start_pose )
    transform_into_memb = rosetta.protocols.membrane.TransformIntoMembraneMover()
    transform_into_memb.apply( start_pose)
    
    memb_jump = start_pose.conformation().membrane_info().membrane_jump()
    memb_length = start_pose.conformation().membrane_info().membrane_core()
    initial_move = rosetta.numeric.xyzVector_double_t(0.0,0.0,memb_length)
    initialize = rosetta.protocols.membrane.TranslationMover(initial_move, memb_jump)

    #rotation mover
    peptide_axis = rosetta.numeric.xyzVector_double_t(0.0,0.0,1.0)
    normal_axis = rosetta.numeric.xyzVector_double_t(0.0,1.0,0.0)
    membrane_center = rosetta.numeric.xyzVector_double_t(0.0,0.0,0.0)
    initial_angle_rotation = 90.0 
    init_rotation = rosetta.protocols.membrane.RotationMover( normal_axis, peptide_axis, membrane_center, memb_jump )


    init_rotation.apply(start_pose)
    initialize.apply(start_pose)
    #start position of the pose is along the membrane layer.
    
    print("score function at beginning: " + str(sfxn(start_pose))) 
    
    mm = MoveMap()
    mm.set_bb(True)
    small = rosetta.protocols.simple_moves.SmallMover()
    small.nmoves(10)
    small.temperature(1.0)

    shear = rosetta.protocols.simple_moves.ShearMover()
    shear.nmoves(10)
    shear.temperature(1.0)
    seq = pyrosetta.rosetta.protocols.moves.SequenceMover()
    seq.add_mover(small)
    seq.add_mover(shear)
    mc = MonteCarlo(start_pose, sfxn, 2.0)
    trial = TrialMover(seq, mc)
    
    loop1_move = 1000
    if(~with_mers):
      loop2_move = 200000
    else:
      loop2_move = 50000
    loop2a_3mer_move = 75000
    loop2b_9mer_move = 75000
    loop3_move = 1000
    
    if(~with_mers):
      score = np.zeros((loop1_move+loop2_move+loop3_move,),dtype=float)
    else:
      score = np.zeros((loop1_move+loop2_move+loop2a_3mer_move+loop2b_9mer_move+loop3_move,),dtype=float)
    
    work_pose = start_pose.clone()
    mc.reset(work_pose)
    
    for i in range(loop1_move):
        trial.apply(work_pose)
        score[i] = sfxn(work_pose)


    work_pose.dump_pdb(peptide_name+"_loop1.pdb")
    
    #adding random torsion angles
    rtor = pyrosetta.rosetta.protocols.simple_moves.RandomTorsionMover(mm, 30, 20) 
    seq.add_mover(rtor)
    
    trial = TrialMover(seq, mc)
    
    for i in range(loop2_move):
        trial.apply(work_pose)
        score[loop1_move + i] = sfxn(work_pose)
        
        
    work_pose.dump_pdb(peptide_name+"_loop2.pdb")   
    ##--------adding fragmets here with glycines at the end-------##
    if(with_mers):
        fragset3 = ConstantLengthFragSet(3)
        fragset3.read_fragment_file(file_3mer)
        fragset9 = ConstantLengthFragSet(9)
        fragset9.read_fragment_file(file_9mer)
        
        mover_3mer = pyrosetta.rosetta.protocols.simple_moves.ClassicFragmentMover(fragset3, mm)
        mover_9mer = pyrosetta.rosetta.protocols.simple_moves.ClassicFragmentMover(fragset9, mm)
        
        repeat_9mer_frags = pyrosetta.rosetta.protocols.moves.RepeatMover(mover_9mer, 1)
        repeat_3mer_frags = pyrosetta.rosetta.protocols.moves.RepeatMover(mover_3mer, 1)
        
        trial9 = TrialMover(repeat_9mer_frags, mc)
        trial3 = TrialMover(repeat_3mer_frags, mc)
        
          
        for i in range(loop2b_9mer_move):
            trial9.apply(work_pose)
            score[loop1_move + loop2_move + i] = sfxn(work_pose)
        
        print("score after 9mer: "+ str(sfxn(work_pose)))
        
        for i in range(loop2a_3mer_move):
            trial3.apply(work_pose)
            score[loop1_move + loop2_move + loop2b_9mer_move + i] = sfxn(work_pose)
        
        print("score after 3mer: "+ str(sfxn(work_pose)))    
         
    ##adding small shears. I have removed the fragmets movers here. 
    seq_ss = pyrosetta.rosetta.protocols.moves.SequenceMover()
    seq_ss.add_mover(small)
    seq_ss.add_mover(shear)
    trial_ss = TrialMover(seq_ss, mc)
    mc.set_temperature(1.0)
    if(~with_mers):
      index_shift = 0
    else:
      index_shift = loop2a_3mer_move + loop2b_9mer_move 
    for i in range(loop3_move):
        trial_ss.apply(work_pose)
        score[loop1_move + loop2_move + index_shift+ i] = sfxn(work_pose)
    
    with open(peptide_name+'_score.dat', 'wb') as f:
        np.save(f, score)

    best_pose = Pose()
    mc.recover_low(best_pose)
    print("best pose score:", str(sfxn(best_pose)))
    
    min_mover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
    min_mover.set_movemap(mm)
    min_mover.score_function(sfxn)
    #min_mover.max_iter(100)
    min_mover.apply(best_pose)
    print("after minimization:", str(sfxn(best_pose)))
    
    tf = TaskFactory()
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepacking())
    
    ex1ex2 = pyrosetta.rosetta.core.pack.task.operation.ExtraRotamersGeneric()
    ex1ex2.ex1( True )
    ex1ex2.ex2( True )
    ex1ex2.ex1aro( True )
    ex1ex2.ex2aro( True )
    tf.push_back(ex1ex2)
    
    packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover()
    packer.task_factory(tf)
    packer.score_function(sfxn)
    packer.nloop(10)
    
    n_iterations = 100
    score_relax = np.zeros((n_iterations,),dtype=float)

    for i in range(n_iterations):
        packer.apply(best_pose)
        score_relax[i] = sfxn(best_pose)
        
    with open(peptide_name+'_relax_score.dat', 'wb') as f:
      np.save(f, score_relax)
    
    tmp_best_pose = best_pose.energies().total_energies().weighted_string_of( sfxn.weights() )
    print(tmp_best_pose)
    best_pose.dump_pdb(peptide_name+'_relax_withmer.pdb')
    
    
def main( args ):

    parser = OptionParser(usage="usage: %prog -f franklin2021 -n pHLIP-v1")
    parser.set_description(main.__doc__)

    parser.add_option('--efxn','-f',
    action = "store",
    help="energy function name", )

    parser.add_option('--name','-n',
    action = "store",
    help="peptide name", )
    
    parser.add_option('--path','-p',
    action = "store",
    help = "path to benchmark",)
    
    
    (options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options
    
    if ( not Options.efxn ):
      print("Missing required Energy function name! Exiting...")
      sys.exit()

    if ( not Options.name ):
      print("Missing required pHLIP peptide name! Exiting...")
      sys.exit()
      
    if ( not Options.path ):
      print("Missing required path to the bechmark folder! Exiting...")
      sys.exit()
  
    #initialize 
    init('-ignore_unrecognized_res -relax:constrain_relax_to_start_coords\
     -ex1 -ex2 -ex1aro -ex2aro -overwrite -relax:default_repeats 10 -pH_mode true -value_pH 8.0')
    
    # Create an energy function
    sfxn = ScoreFunction()
    #sfxn = get_fa_scorefxn() #this is for default score function
    if( Options.efxn == "franklin2019" ):
        
        sfxn = create_score_function( Options.efxn )
        

    elif( Options.efxn == "score12" or Options.efxn == "score07"):

        score_file = Options.path + "tests/python/" + energy_fxn + ".wts"
        
        sfxn = create_score_function( score_file )
    elif( Options.efxn =="franklin2021"):
        
        sfxn = create_score_function( Options.efxn )
        #fa_w_b    fa_imm_elec f_elec_lipidlayer 
        #1.375		-0.143		0.286
        ##0.933		0.733		0.094
        #1.179		-0.015		-0.107
        #1.571		-0.025		-0.179
        #1.629	    -1.561	    0.136
        #1.484	    0.571	    -0.462

        sfxn.set_weight(fa_water_to_bilayer, 1.629)
        sfxn.set_weight(fa_imm_elec, -1.561)
        sfxn.set_weight(f_elec_lipidlayer, 0.136)
        sfxn.set_weight(fa_elec, 1.00)
    else:
        #uses the default efxn
        sfxn = create_score_function()

    run_abinitio_structure_prediction( sfxn, Options.name, Options.path )
    
if __name__ == "__main__": main( sys.argv )    
    
    
    
    
    
        
        