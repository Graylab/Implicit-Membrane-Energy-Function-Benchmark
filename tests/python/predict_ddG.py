# @file: predict_ddG.py
""" Predict the ddG of mutation

This script computes the energetic cost of single point
mutations in membrane proteins. To account for the conformational
change, we repack side chains within 8 Angstroms of the host site.
This script generates data for Test 7.

Authors: 
    Rebecca Alford <ralford3@jhu.edu> 

Example: 
    $ import predict_ddG
    $ predict_ddG.run_ddG_of_mutation_calc( 
    config, energy_fxn, targets_dir, native_pdb, 
    native_span, reference_data )
        
Arguments: 
    - config: Container with path to benchmark and rosetta files
    - energy_fxn: Weights file for energy function of interest
    - targets_dir: Name of output data directory
    - native_pdb: Conformation of the native state
    - native_span: Transmembrane topology of the native state
    - reference_data: Text file containing experimental ref data. 

Requirements: 
    - PyRosetta4 and Python 3.6 or 3.7
"""

from pyrosetta import *
from pyrosetta.teaching import *
from string import Template
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

import sys, os
import random

import pyrosetta.rosetta.protocols.membrane
from pyrosetta.rosetta.utility import vector1_bool
from pyrosetta.rosetta.core.chemical import aa_from_oneletter_code
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.pose import PDBInfo
from pyrosetta.rosetta.core.chemical import VariantType
from pyrosetta.rosetta.core.pack.task import TaskFactory
import pyrosetta.rosetta.core.select.residue_selector as residue_selector

# Make a dictionary for classifying residue types (yes I know, global vars are bad practice)
classification = {'A': "nonpolar", 'C' : "special", 'D' : "n-charged", 'E' : "n-charged", \
                         'F' : "aromatic", 'G' : "special", 'H' : "p-charged", 'I' : "nonpolar", \
                         'K' : "p-charged", 'L' : "nonpolar", 'M' : "nonpolar", 'N' : "polar", \
                         'P' : "special", 'Q' : "polar", 'R' : "p-charged", 'S' : "polar", \
                         'T' : "polar", 'V' : "nonpolar", 'W' : "aromatic", 'Y' : "aromatic" }

# @brief Replace the residue at <resid> in <pose> with <new_res> and allows
# repacking within a given <pack_radius>
def mutate_residue( pose, mutant_position, mutant_aa, pack_radius, pack_scorefxn ):

    if pose.is_fullatom() == False:
       IOError( 'mutate_residue only works with fullatom poses' )

    test_pose = Pose()
    test_pose.assign( pose )

    # Create a packer task (standard)
    task = TaskFactory.create_packer_task( test_pose )

    # the Vector1 of booleans (a specific object) is needed for specifying the
    #    mutation, this demonstrates another more direct method of setting
    #    PackerTask options for design
    aa_bool = vector1_bool()

    # PyRosetta uses several ways of tracking amino acids (ResidueTypes)
    # the numbers 1-20 correspond individually to the 20 proteogenic amino acids
    # aa_from_oneletter returns the integer representation of an amino acid
    #    from its one letter code
    # convert mutant_aa to its integer representation
    mutant_aa = aa_from_oneletter_code( mutant_aa )

    # mutation is performed by using a PackerTask with only the mutant
    #    amino acid available during design
    # to do this, construct a Vector1 of booleans indicating which amino acid
    #    (by its numerical designation, see above) to allow
    for i in range( 1 , 21 ):
      # in Python, logical expression are evaluated with priority, thus the
      #    line below appends to aa_bool the truth (True or False) of the
      #    statement i == mutant_aa
       aa_bool.append( i == mutant_aa )

    # modify the mutating residue's assignment in the PackerTask using the
    #    Vector1 of booleans across the proteogenic amino acids
    task.nonconst_residue_task( mutant_position).restrict_absent_canonical_aas( aa_bool )

    # prevent residues from packing by setting the per-residue "options" of
    #    the PackerTask
    center = pose.residue( mutant_position ).nbr_atom_xyz()
    for i in range( 1, pose.total_residue() + 1 ):
        dist = center.distance_squared( test_pose.residue( i ).nbr_atom_xyz() )
       # only pack the mutating residue and any within the pack_radius
        if i != mutant_position and dist > pow( float( pack_radius ), 2 ) :
           task.nonconst_residue_task( i ).prevent_repacking()
        elif i != mutant_position and dist <= pow( float( pack_radius ), 2 ) :
            task.nonconst_residue_task( i ).restrict_to_repacking()
        # apply the mutation and pack nearby residues
    packer = PackRotamersMover( pack_scorefxn , task )
    packer.apply( test_pose )

    return test_pose

def mutate_residue_withgreenpacker( pose, mutant_position, mutant_aa, pack_radius, pack_scorefxn ):

    sel_res = residue_selector.ResidueIndexSelector(str(mutant_position))

    nbr_selector = residue_selector.NeighborhoodResidueSelector()
    nbr_selector.set_distance(pack_radius)
    nbr_selector.set_focus_selector(sel_res)
    nbr_selector.set_include_focus_in_subset(False)

    sel_res_nbr = residue_selector.OrResidueSelector(selector1=sel_res, selector2=nbr_selector)
    
    neighbor_res = pyrosetta.rosetta.core.select.get_residues_from_subset( sel_res_nbr.apply(pose) )
    #print(neighbor_res)

    notsel_res_nbr = residue_selector.NotResidueSelector(sel_res_nbr)
    #print(pyrosetta.rosetta.core.select.get_residues_from_subset(not_chain_H_sel_res_nbr.apply(pose))) 

    target_AA = mutant_aa
    print(target_AA)

    test_pose = Pose()
    test_pose.assign( pose )

    
    init_cmd = pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline()  

    #generate a taskfactory
    mut_res = TaskFactory()
    mut_res.push_back(init_cmd)
    mut_res.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
    # Disable design (i.e. repack only) on not_chain_A_cys_res
    #mut_res.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(\
    #    pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), not_chain_H_sel_res))
    mut_res.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(\
        pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT(), notsel_res_nbr))
    mut_res.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(\
        pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), nbr_selector))

    aa_to_design = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
    aa_to_design.aas_to_keep(target_AA)
    mut_res.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(aa_to_design, sel_res))
    
    ex1ex2 = pyrosetta.rosetta.core.pack.task.operation.ExtraRotamersGeneric()
    ex1ex2.ex1( True )
    ex1ex2.ex2( True )
    ex1ex2.ex1aro( True )
    ex1ex2.ex2aro( True )
    mut_res.push_back(ex1ex2)

    group_ids = pyrosetta.rosetta.utility.vector1_unsigned_long()
    
    packer_task = mut_res.create_task_and_apply_taskoperations(pose)
    #print(packer_task)
    #can be shortened, but keeping it general
    for i in range(1,pose.size()+1):
        if( packer_task.residue_task(i).being_designed() & pose.residue(i).is_protein() ):
            group_ids.append(0)
        else:
            group_ids.append(1)
    
    #print(group_ids)

    user_defined_group_discriminator = pyrosetta.rosetta.protocols.minimization_packing.UserDefinedGroupDiscriminator()
    user_defined_group_discriminator.set_group_ids( group_ids )
 
 
    greenpacker = pyrosetta.rosetta.protocols.minimization_packing.GreenPacker()
    greenpacker.set_reference_round_task_factory( mut_res )
    greenpacker.set_group_discriminator( user_defined_group_discriminator )
    greenpacker.set_scorefunction( pack_scorefxn )
    greenpacker.set_task_factory( mut_res )

    greenpacker.apply( test_pose )

    print("score after GP:", pack_scorefxn( test_pose ))

    for i in range(1, pose.size()+1):
        if( pose.residue(i).name1() != test_pose.residue(i).name1() ):
            print(str(i) + " before: " + pose.residue(i).name1() + " after: " + test_pose.residue(i).name1())
   

    print("Green Packer applied")

    mm_fa = pyrosetta.rosetta.core.kinematics.MoveMap()
    
    mm_fa.set_jump(False)
    for i in range(1, test_pose.size()+1):
        if( packer_task.residue_task(i).being_packed() & test_pose.residue(i).is_protein() ):
            mm_fa.set_chi(i,True)
            mm_fa.set_bb(i,False)
            #only relax the side chains not the bb
        else:
            mm_fa.set_chi(i,False)
            mm_fa.set_bb(False)

    min_mover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
    min_mover.set_movemap(mm_fa)
    min_mover.max_iter(1000)
    #max iter default is 200/not changing
    min_mover.score_function(pack_scorefxn)
    min_mover.apply(test_pose)

    print("score after minimization:", pack_scorefxn( test_pose ))
    return( test_pose )
    


#def mutate_residue( pose, mutant_position, mutant_aa, pack_radius, pack_scorefxn ):
#    if pose.is_fullatom() == False:
#        IOError( 'mutate_residue only works with fullatom poses' )
#    from pyrosetta.toolbox.mutants import mutate_residue
#    test_pose = pose.clone()
#    mutate_residue( test_pose, int(mutant_position), mutant_aa, pack_radius, pack_scorefxn )
#    return test_pose

def run_ddG_of_mutation_calc( config, energy_fxn, testname, pdb, spanfile, mutation_list ): 
    
    # Initialize Pyrosetta with const options
    # option_string = "-run:constant_seed -in:ignore_unrecognized_res"
    option_string = "-in:ignore_unrecognized_res"
    
    init( extra_options=option_string )

    # Append base path to pdb, spanfile, and mutation_list
    pdb = config.benchmark_path + "targets/stability/" + testname + "/" + pdb
    spanfile = config.benchmark_path + "targets/stability/" + testname + "/" + spanfile
    mutation_list = config.benchmark_path + "targets/stability/" + testname + "/" + mutation_list 

    fa_water_bilayer_wts = 1.000
    fa_imm_elec_wts = 0.01
    f_elec_lipidlayer_wts = 0.128
    # Read database file including mutations (space delimited)
    # Expected header format: Nat Pos Mut PDb Spanfile pH exp_ddG double_mut
    with open( mutation_list, 'r' ) as f:
        content = f.readlines()
    content = [ x.strip() for x in content ]
    mutation_list = [ x.split(' ') for x in content ]

    # Create an energy function
    sfxn = ScoreFunction()
    #sfxn = get_fa_scorefxn() #this is for default score function
    if( energy_fxn == "franklin2019" ):
        option_string = option_string + " -mp:lipids:temperature 37.0 -mp:lipids:composition DLPC -mp:lipids:has_pore true -ex1 -ex2 -ex1aro -ex2aro"
        sfxn = create_score_function( energy_fxn )
        sfxn.set_weight(fa_water_to_bilayer, 1.5)

    elif( energy_fxn == "score12" or energy_fxn == "score07"):

        score_file = config.benchmark_path + "tests/python/" + energy_fxn + ".wts"
        option_string = option_string + "-ex1 -ex2 -ex1aro -ex2aro"# -restore_pre_talaris_2013_behavior"
        #-mp:lipids:temperature 37.0 -mp:lipids:composition DLPC -mp:lipids:has_pore true
        #"-score:weights" + score_file +
        sfxn = create_score_function( score_file )
    elif( energy_fxn =="franklin2021"):
        option_string = option_string + " -mp:lipids:temperature 37.0 -mp:lipids:composition DLPC -ex1 -ex2 -ex1aro -ex2aro -pH_mode False"
        #score_file = config.benchmark_path + "tests/python/" + energy_fxn + ".wts"
        sfxn = create_score_function( energy_fxn )
        #fa_w_b    fa_imm_elec f_elec_lipidlayer 
        #1.375		-0.143		0.286
        ##0.933		0.733		0.094
        #1.179		-0.015		-0.107
        #1.571		-0.025		-0.179
        #1.629	    -1.561	    0.136
        #1.484	    0.571	    -0.462
        # fa_water_to_bilayer_combined	fa_imm_elec_combined	f_elec_lipidlayer_combined
		
        # sfxn.set_weight(fa_water_to_bilayer, 0.966)
        # sfxn.set_weight(fa_imm_elec, 0.379)
        # sfxn.set_weight(f_elec_lipidlayer, 0.016)
        # sfxn.set_weight(fa_elec, 1.00)
        # sfxn.set_weight(menv_pH, 0.00)
        # 0.938	0.320	0.185
        # 1.044	0.107	0.125
        # 0.910	0.121	0.046
        # 0.863 0.001 0.152
        
        	
        
        #the pyrosetta version in the environment doesnt have the weights included. 
        sfxn.set_weight(fa_water_to_bilayer, fa_water_bilayer_wts)
        sfxn.set_weight(fa_imm_elec, fa_imm_elec_wts)
        sfxn.set_weight(f_elec_lipidlayer, f_elec_lipidlayer_wts)
        sfxn.set_weight(fa_elec, 1.00)
        sfxn.set_weight(menv_pH, 0.00)
        
    else:
        option_string = option_string + " -mp:lipids:temperature 37.0 -mp:lipids:composition DLPC -mp:lipids:has_pore true -ex1 -ex2 -ex1aro -ex2aro -pH_mode False"
        #score_file = config.benchmark_path + "tests/python/" + energy_fxn + ".wts"
        sfxn = create_score_function( energy_fxn )
        


    #sfxn = create_score_function( energy_fxn )
    ##score12 is giving error, ask Rebecca which is the right word, and how to find it. 
    ##Time being set_weight score function given in the suuple of the biophysJ paper. 
    ##sfxn = ScoreFunction()
    print(sfxn.get_weight(fa_atr))
    print("score fxn components")
    print(sfxn)

    # Set the repack radius from the option system
    repack_radius = 6.0

    # Make output directories
    outdir = config.benchmark_path + "data/" 
    if ( not os.path.isdir( outdir ) ): 
        os.system( "mkdir " + outdir )

    outdir = outdir + energy_fxn 
    if ( not os.path.isdir( outdir ) ): 
        os.system( "mkdir " + outdir )

    if(energy_fxn =="franklin2021"):
        outdir = outdir + "/ddG-of-mutation/test8_afternheavyatomcorrection_changedddG/"
        outdir = outdir + 'fa_wb_'+str(fa_water_bilayer_wts)+'_felecbilayer_'+str(f_elec_lipidlayer_wts)+'_fimm_'+str(fa_imm_elec_wts)
        print(outdir)
    else:
        outdir = outdir + "/ddG-of-mutation/"
        print(outdir)
    
    # _based_vjul8_fa_wb_0.966_felecbilayer_0.016_fimm_0.379/withpore"
    if ( not os.path.isdir( outdir ) ): 
        os.system( "mkdir " + outdir )

    targetdir = outdir + "/" + testname 
    if ( not os.path.isdir( targetdir ) ): 
        os.system( "mkdir " + targetdir ) 
    os.chdir( targetdir )

    # Setup output file
    outfile = targetdir + "/ddG_" + energy_fxn + ".dat"
    f = open( outfile, 'a' )
    f.write( "Nat Pos Mut experimental_ddG predicted_ddG class depth\n" )
    
    pose = pose_from_pdb( pdb )
    print_score_labels_to_file( pose, sfxn, targetdir+"/breakdown.sc" )
    print_score_labels_to_file( pose, sfxn, targetdir+"/breakdown_native.sc" )

    for entry in mutation_list:

        # Sanity check that the PDB file path and spanfile path exist
        if ( not os.path.isfile(pdb) or not os.path.isfile(spanfile) ):
            print(pdb + " " + spanfile)
            print("Path to PDB file or spanfile is invalid!")
            sys.exit()

        # Read PDB from table - note, must contain an absolute path
        pose = pose_from_pdb( pdb )

        # Add membrane to pose
        add_memb = rosetta.protocols.membrane.AddMembraneMover( spanfile )
        add_memb.apply( pose )

        # Setup in a topology based membrane
        init_mem_pos = rosetta.protocols.membrane.MembranePositionFromTopologyMover()
        init_mem_pos.apply( pose )

        # Calculate the native state, based on whether or not there is a shift to alanine first
        if ( entry[7] == "Y" ):
            if(testname == "C3_OmpLA_aro_ddGs"):
                repacked_native = mutate_residue_withgreenpacker( pose, int( entry[1] ), 'A', repack_radius, sfxn )
            else:
                repacked_native = mutate_residue_withgreenpacker( pose, int( entry[1] ), 'A', repack_radius, sfxn )
        else:
            native_res = pose.residue( int( entry[1] ) ).name1()

            if(testname == "C3_OmpLA_aro_ddGs"):
                repacked_native = mutate_residue_withgreenpacker( pose, int( entry[1] ), native_res, repack_radius, sfxn )
            else:
                repacked_native = mutate_residue_withgreenpacker( pose, int( entry[1] ), native_res, repack_radius, sfxn )
        
        print("before mutation:" + pose.residue( int( entry[1] ) ).name() )
        print("after mutation:" + repacked_native.residue( int( entry[1] ) ).name() )
 
        # Calculate the ddG of mutation for the given position
        ddG_of_mutation = compute_ddG( repacked_native, sfxn, int( entry[1] ), entry[2], repack_radius, targetdir, testname )
        print(ddG_of_mutation)
        print_score_labels_to_file( repacked_native, sfxn, "dummy" )

        # Calculate some additional classifications for the mutation
        depth = pose.conformation().membrane_info().residue_z_position( pose.conformation(), int( entry[1] ) )
        rsd_class = classification[ entry[2] ]

        # Write ddG data to output file
        outstr = Template( "$native $pos $mutant $exp_val $predicted $rsd_class $depth" )
        output = outstr.substitute( native=entry[0], pos=entry[1], mutant=entry[2], exp_val=entry[6], predicted=ddG_of_mutation, rsd_class=rsd_class, depth=round(depth,3))
        f.write( output + "\n" )


    f.close()

## @brief Compute ddG of mutation in a protein at specified residue and AA position
def compute_ddG( pose, sfxn, resnum, aa, repack_radius, targetdir, testname ):
    import time

    # Score the native pose and grab the native AA
    native_score = sfxn( pose )
    native_aa = pose.residue( resnum ).name1()

    # Perform the mutation at residue <resnum> to amino acid <aa> and score
    if(testname == "C3_OmpLA_aro_ddGs"):
        mutated_pose = mutate_residue_withgreenpacker( pose, resnum, aa, repack_radius, sfxn )
    else:
        mutated_pose = mutate_residue( pose, resnum, aa, repack_radius, sfxn )
    
    mutant_score = sfxn( mutated_pose )
    
    file_name = targetdir + "/output_"+ str( resnum ) + aa +"_pose.pdb"
    print("before mutation:" + pose.residue( resnum ).name() )	
    print("after mutation:" + mutated_pose.residue( resnum ).name() )
    #time.sleep(15)
    mutated_pose.dump_pdb(file_name)
    
    #with open( targetdir + "/breakdown.sc", 'a' ) as f:
    #      f.write( print("native_res_mutant") + print(sfxn) + "\n" )
    #f.close()
    
    
    #print_score_labels_to_file( mutated_pose, sfxn, targetdir+"/breakdown.sc" )
    print_ddG_breakdown( pose, mutated_pose, sfxn, resnum, aa, targetdir + "/breakdown.sc" )
    #with open( targetdir + "/breakdown_native.sc", 'a' ) as f:
    #      f.write( print("native_res_mutant") + print(sfxn) + "\n" )
    #f.close()

    #print_score_labels_to_file( pose, sfxn, targetdir + "/breakdown_native.sc" )
    print_ddG_breakdown( pose, pose, sfxn, resnum, aa, targetdir + "/breakdown_native.sc" )
    ##fn="/breakdown.sc"
    ##with open( fn, 'a' ) as f:
    ##    f.write( print("native_res_mutant") + print(sfxn) + "\n" )
    ##f.close()    
# Calculate the ddG in place
    ddG = round( mutant_score - native_score, 3 )
    return ddG

###############################################################################
#@brief Print ddG breakdown from the pose
# Extract weighted energies from the native and mutated pose. Calculate the ddG
# of each and print the component-wise ddG vlaues
def print_ddG_breakdown( native_pose, mutated_pose, sfxn, resnum, aa, fn ): 

    # Extract scores
    tmp_native = native_pose.energies().total_energies().weighted_string_of( sfxn.weights() )
    tmp_mutant = mutated_pose.energies().total_energies().weighted_string_of( sfxn.weights() )

    # Parse out scores
    array_native = list(filter( None, tmp_native.split(' ') ))
    array_mutant = list(filter( None, tmp_mutant.split(' ') ))

    # Pull out only the scores from these arrays
    native_scores = []
    for i in range( len(array_native) ): 
        if ( i % 2 != 0 ): 
            native_scores.append( float( array_native[i] ) )

    mutant_scores = []
    for i in range( len(array_mutant) ): 
        if ( i % 2 != 0 ): 
            mutant_scores.append( float( array_mutant[i] ) )

    # Make a label for the mutation
    native_res = native_pose.residue( int( resnum ) ).name1()
    mut_label = native_res + str(resnum) + aa

    # Calculate ddG of individual components
    ddGs = []
    ddGs.append( mut_label )
    for i in range( len( mutant_scores ) ): 
        ddG_component = mutant_scores[i] - native_scores[i]
        ddGs.append( round( ddG_component, 3 ) )

    ddGs_str = convert_array_to_str( ddGs )
    
    with open( fn, 'a' ) as f:
          ##f.write( print("native_res_mutant") + print(sfxn) + "\n" )
          f.write( ddGs_str + "\n" )
    f.close()

###############################################################################
#@brief Get header for ddG breakdown output
# Save the score labels, to be printed at the top of the output breakdown file
def print_score_labels_to_file( native_pose, sfxn, fn ): 

    tmp_native = native_pose.energies().total_energies().weighted_string_of( sfxn.weights() )
    array_native = list(filter( None, tmp_native.split(' ') ))
    labels = []
    labels.append( 'mutation ' ) # Append field for mutation label
    for i in range( len(array_native) ): 
        if ( i % 2 == 0 ): 
            labels.append( array_native[i].translate(':') )

    labels_str = convert_array_to_str( labels )
    print(labels_str)
    with open( fn, 'a' ) as f:
        f.write( labels_str + "\n" )
    f.close()

###############################################################################
#@brief Convert an array to a space deliminted string
# Save the score labels, to be printed at the top of the output breakdown file
def convert_array_to_str( array ): 

    linestr = ""
    for elem in array: 
        if ( linestr == "" ): 
            linestr = linestr + str( elem )
        else: 
            linestr = linestr + " " + str( elem )

    return linestr

def get_energy_components( native_pose, mutated_pose, sfxn, resnum, aa ): 

    # Extract & parse scores
    tmp_native = native_pose.energies().total_energies().weighted_string_of( sfxn.weights() )
    tmp_mutant = mutated_pose.energies().total_energies().weighted_string_of( sfxn.weights() )
    array_native = list(filter( None, tmp_native.split(' ') ))
    array_mutant = list(filter( None, tmp_mutant.split(' ') ))

    # Pull out only the scores from these arrays
    native_scores = []
    for i in range( len(array_native) ): 
        if ( i % 2 != 0 ): 
            native_scores.append( float( array_native[i] ) )

    mutant_scores = []
    for i in range( len(array_mutant) ): 
        if ( i % 2 != 0 ): 
            mutant_scores.append( float( array_mutant[i] ) )

    # Make a label for the mutation
    native_res = native_pose.residue( int( resnum ) ).name1()
    mut_label = native_res + str(resnum) + aa

    # Calculate ddG of individual components
    ddGs = []
    ddGs.append( mut_label )
    for i in range( len( mutant_scores ) ): 
        ddG_component = mutant_scores[i] - native_scores[i]
        ddGs.append( round( ddG_component, 3 ) )

    # Get labels
    labels = []
    for i in range( len(array_native) ): 
        if ( i % 2 == 0 ): 
            labels.append( array_native[i].translate(':') )

    return labels, ddGs

if __name__ == "__main__" : main(sys.argv)
