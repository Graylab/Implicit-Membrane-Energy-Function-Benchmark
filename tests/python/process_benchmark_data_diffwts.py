#!usr/bin/env python
""" Master script for processing sequence and structure benchmarks

This module analyzes data from tests that probe sequence and 
structure features. This is the second step of three for executing
and evaluating the scientific testds. 

Authors: 
Rebecca Alford <ralford3@jhu.edu> 

Example: 
$ python analyze_benchmark_data.py --energy_fxn franklin2019
--which_tests all --restore_talaris False

Arguments: 
- energy_fxn: Weights file for energy function of interest
- which_tests: Run all tests, or specify by comma-separated list

Requirements: 
- Rosetta release 246 or greater
- PyRosetta4 for Python 3.6 or 3.7
- KinkFinder
"""

import sys, os, random, read_config, hpc_util
import make_asymm_docked_complexes
import make_designed_protein_scaffolds
import make_refined_decoys
import predict_side_chain_distribution

from string import Template
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

all_tests = [ "sc-distribution", "ddG-of-insertion", "ddG-of-mutation", "ddG-of-pH-insertion", "decoy-discrimination", "tm-peptide-tiolt", "helix-kinks", "hydrophobic-length", "adsorbed-peptide-tilt-angle", "protein-protein-docking", "protein-tilt-angle", "sequence-recovery" ]
sequence_tests = [ "sc-distribution", "sequence-recovery" ]
structure_tests = [ "decoy-discrimination", "helix-kinks", "protein-protein-docking" ]

def convert_spanfile_to_string( spanfile ): 

	with open( spanfile, 'rt' ) as f: 
		spans = f.readlines()
		spans = [ x.strip() for x in spans ]

	# Parse tm regions, skipping the first four header lines
	span_info = []
	for tm_span in range(4, len(spans)): 
		single_span_data = spans[tm_span].split( ' ' )
		single_span_info = []
		for entry in single_span_data: 
			if ( entry != "" ): 
				single_span_info.append( entry )
		span_info.append( single_span_info )

	# Convert list of spans into a string that will be used as input for kink finder
	spanstr = ""
	for s in range(0, len(span_info)):
		if ( s != len(span_info)-1 ):
			spanstr = spanstr + span_info[s][0] + "-" + span_info[s][1] + " "
		else: 
			spanstr = spanstr + span_info[s][0] + "-" + span_info[s][1] 

	return spanstr

def main( args ): 
	# mp_env = False
	
 	# Read options from the command line
	parser = OptionParser(usage="usage %prog --energy_fxn franklin2019 --which_tests all" )
	parser.set_description(main.__doc__)

	parser.add_option( '--energy_fxn', '-e', action="store", help="Name of energy function weights file", )
	parser.add_option( '--which_tests', '-w', action="store", help="Specify tests run (comma separated list)", )

	(options, args) = parser.parse_args(args=args[1:])
	global Options
	Options = options
	print("after global options")

	# Check that required options have been provided
	if ( not Options.energy_fxn or not Options.which_tests ): 
		print("Missing required options --energy_fxn and/or --which_tests" )
		sys.exit()

	# Set restore variable based on energy function type
	restore = True
	if ( Options.energy_fxn == "franklin2019" or Options.energy_fxn == "ref2015" or Options.energy_fxn == "franklin2021"): 
		restore = False 

	# Read path configuration file
	config = read_config.read_config()
	print("after read config")
	# Check test categories
	test_names = []
	if ( Options.which_tests == "all" ): 
		test_names = all_tests
	else: 
		test_names = Options.which_tests.split(",")
		# check that all names are valid
		for name in test_names: 
			if name not in all_tests: 
				sys.exit( "No such test " + name + ". Exiting!" )

	wt_fa_water_to_bilayer = [0.863] 
	wt_fa_imm_elec = [0.001]
	wt_f_elec_bilayer = [0.152] 

	# Test #8: Sequence recovery calculation
	if ( "sequence-recovery" in test_names ): 
		
		
		i=0

		outdir = '/home/rsamant2/scratch16-jgray21/rsamant2/' 
		if(Options.energy_fxn=='franklin2021'):
			datadir = outdir + "data/" + Options.energy_fxn + "/"+ Options.which_tests + "/weights_from_test8_heavyatomcorrection_changedddG_disallowC"
			datadir = datadir + "/fa_wb_" + str(wt_fa_water_to_bilayer[i]) + "_felecbilayer_" + str(\
					wt_f_elec_bilayer[i]) + "_fimm_" + str(wt_fa_imm_elec[i])
		else:
			# datadir = outdir + "data/" + Options.energy_fxn + "/"+ Options.which_tests 
			#Franklin2019 for no C design
			datadir = outdir + "data/" + Options.energy_fxn + "/"+ Options.which_tests + "/weights_from_test8_heavyatomcorrection_changedddG_disallowC"
			datadir = datadir + "/fa_wb_" + str(wt_fa_water_to_bilayer[i]) + "_felecbilayer_" + str(\
					wt_f_elec_bilayer[i]) + "_fimm_" + str(wt_fa_imm_elec[i])
   
		# Make list of native and designed PDB files
		#datadir = config.benchmark_path + "data/" + Options.energy_fxn + "/sequence-recovery/highdGforDE"
		print(datadir)
		#os.chdir(datadir)
			
		#os.chdir( datadir )
		os.chdir(datadir)
		#os.system( "ls */*_0001.pdb > designed.list" )
		filename =  "design.list"
		with open( filename,'rt' ) as f: 
			contents = f.readlines()
			contents = [ x.strip() for x in contents ]
			pdbid = [ x.split("/")[0] for x in contents ]
			#only for disallowC case. temporary remove it later. 
			# print(pdbid)
			for i in range(len(pdbid)):
				if(pdbid[i]=='..'):
					pdbid[i] = contents[i].split("/")[4]
		# print(pdbid)
		with open( "natives.list", 'wt' ) as f: 
			basedir = config.benchmark_path + "targets/design/"
			for pdb in pdbid: 
				pdbpath = basedir + pdb + "/" + pdb + "_tr_ignorechain.pdb\n"
				f.write( pdbpath )
		# sys.exit()
		# Run mp_seqrecov application
		executable = config.rosetta_path + "/mp_seqrecov." + config.platform + config.compiler + config.buildenv
		print(executable)
		output_file = Options.energy_fxn + "_seqrecov.txt"
		database_file = config.rosetta_path + "/../../database"
		s = Template( " -overwrite -database $databasefile -native_pdb_list $natives -redesign_pdb_list $designed -seq_recov_filename $outfile -in:ignore_unrecognized_res -read_only_ATOM_entries")
		arguments = s.substitute( natives="natives.list", designed=filename, outfile=output_file, databasefile=database_file )
		if ( restore == False ): 
			arguments = arguments + " -mp:lipids:composition DLPC -mp:lipids:temperature 37"
			print(arguments)
   
   ##the run fails if we pass score functions, which is weird. 
   #calcultaing sequence recovery should not need any score function. 
			# if( wt_fa_imm_elec[i]<0 or wt_f_elec_bilayer[i]<0 ):
			# 	filename = 'flag_file'
			# 	#the line below should be "w" not "a" or append mode. 
			# 	with open( filename, 'w' ) as f: 
			# 		f.write(  "-set_weights fa_water_to_bilayer "+ str(wt_fa_water_to_bilayer[i]) + "\n" )
			# 		f.write(  "-set_weights f_elec_lipidlayer "+ str(wt_f_elec_bilayer[i]) + "\n" )
			# 		f.write(  "-set_weights fa_imm_elec "+ str(wt_fa_imm_elec[i]) + "\n" )
			# 		if(mp_env):
			# 			f.write(  "-set_weights fa_mpenv_smooth "+ "0.50" + "\n" )
			# 		f.close()
			# 		os.system( "chmod +x " + filename )
			# 		arguments = arguments + " @flag_file"
			# 		print(arguments)
			# else:
			# 	arguments = arguments + " -set_weights fa_water_to_bilayer " + \
			# 		str(wt_fa_water_to_bilayer[i]) + " -set_weights f_elec_lipidlayer " + str(\
			# 		wt_f_elec_bilayer[i]) + " -set_weights fa_imm_elec " + str(wt_fa_imm_elec[i])
			# 	if(mp_env):
			# 		arguments = arguments + " -set_weights fa_mpenv_smooth "+ "0.50" 

		jobname = 'process_sequence_recovery'

		print(executable + arguments)
		jobfile = hpc_util.make_jobfile(
                            datadir, jobname, executable, arguments)
                        
		# os.system( executable + arguments )

		# Process the sequence recovery data
		# os.system( "python3 process_protein_design_results.py --energy_fxn " + Options.energy_fxn + " --basedir " + config.benchmark_path + " --seqrecov_file " + output_file )
		# command_system = "python3 " + config.benchmark_path + "tests/python/process_protein_design_results.py" + " --energy_fxn " + Options.energy_fxn + " --basedir " + datadir + " --seqrecov_file " + output_file
		# print(command_system)
		# os.system( command_system )
		#   python3 process_protein_design_results.py --energy_fxn 
  		# 'franklin2021' --basedir '/home/rsamant2/scratch16-jgray21/rsamant2/data/franklin2021/sequence-recovery/weights_from_test7_heavyatomcorrection_disallowC/fa_wb_0.863_felecbilayer_0.152_fimm_0.001/' --seqrecov_file 'franklin2021_seqrecov.txt'

	# Test #9: Side chain distribution calculations
	if ( "sc-distribution" in test_names ):
		
		if(Options.energy_fxn=='franklin2021'):
			outdir = '/home/rsamant2/scratch16-jgray21/rsamant2/' 
			datadir = outdir + "data/" + Options.energy_fxn + "/"+ "sequence-recovery" + "/weights_from_test8_heavyatomcorrection_changedddG_disallowC"
			datadir = datadir + "/fa_wb_" + str(wt_fa_water_to_bilayer[0]) + "_felecbilayer_" + str(\
					wt_f_elec_bilayer[0]) + "_fimm_" + str(wt_fa_imm_elec[0]) 
			os.chdir(datadir)
  
		elif(Options.energy_fxn=='proteinmpnn'):
			outdir = '/home/rsamant2/scratch16-jgray21/rsamant2/' 
			datadir=outdir+"data/"+ "franklin2021/" + "/" + "sequence-recovery/" + "protein_mpnn/no_design_C"
			os.chdir(datadir)
		else:
			outdir = '/home/rsamant2/scratch16-jgray21/rsamant2/' 
			# for Franklin2019
			# datadir = outdir + "data/" + Options.energy_fxn + "/"+ "sequence-recovery/"
			# for franklin2019 without designing C
			datadir = outdir + "data/" + Options.energy_fxn + "/"+ "sequence-recovery/" + "/weights_from_test8_heavyatomcorrection_changedddG_disallowC"
			datadir = datadir + "/fa_wb_" + str(wt_fa_water_to_bilayer[0]) + "_felecbilayer_" + str(\
					wt_f_elec_bilayer[0]) + "_fimm_" + str(wt_fa_imm_elec[0]) 
			
			 
			os.chdir(datadir)
  		# datadir = datadir 

		##writing the native and design list for alpha and beta mp
		# filename =  "designed_best_score.list"
		# with open( filename,'rt' ) as f: 
		# 	contents = f.readlines()
		# 	contents = [ x.strip() for x in contents ]
		# 	pdbid = [ x.split("/")[0] for x in contents ]
		
		# filename_alpha = config.benchmark_path + "targets/design/" + "alpha_monomer_chains.list"
		# filename_beta = config.benchmark_path + "targets/design/" + "beta_monomer_chains.list"
	
		# with open(filename_alpha,'rt') as f:
		# 	contents = f.readlines()
		# 	contents = [ x.strip() for x in contents ]
		# 	pdbid_alpha = [ x.split("/")[0] for x in contents ]

		# with open(filename_beta,'rt') as f:
		# 	contents = f.readlines()
		# 	contents = [ x.strip() for x in contents ]
		# 	pdbid_beta = [ x.split("/")[0] for x in contents ]


		# basedir = config.benchmark_path + "targets/design/"
		# for pdb in pdbid: 
		# 	index = pdbid.index(pdb)
		# 	print(pdb)
		# 	if(pdb in pdbid_alpha):
		# 		pdbpath = basedir + pdb + "/" + pdb + "_tr_ignorechain.pdb\n"
		# 		open("natives_alpha.list", 'a+').write(pdbpath)
		# 		design_path = open("designed_best_score.list",'r').readlines()
		# 		open("design_alpha.list", 'a+').write(design_path[index])
	
		# 	elif(pdb in pdbid_beta):
		# 		pdbpath = basedir + pdb + "/" + pdb + "_tr_ignorechain.pdb\n"
		# 		open("natives_beta.list", 'a+').write(pdbpath)
		# 		design_path = open("designed_best_score.list",'r').readlines()
		# 		open("design_beta.list", 'a+').write(design_path[index])
		# 	else:
		# 		print("this pbd needs classification: ", pdb)
     
  		# Check for existence of designed and native lists
		# datadir = config.benchmark_path + "data/" + Options.energy_fxn + "/sequence-recovery/highdGforDE"
		# os.chdir( datadir )
		# if ( not os.path.isfile( "natives.list") or not os.path.isfile( "designed.list" ) ): 
		# 	os.system( "ls */*_0001.pdb > designed.list" )
		# 	with open( "designed.list",'rt' ) as f: 
		# 		contents = f.readlines()
		# 		contents = [ x.strip() for x in contents ]
		# 		pdbid = [ x.split("/")[0] for x in contents ]

		# 	with open( "natives.list", 'wt' ) as f: 
		# 		basedir = config.benchmark_path + "targets/design/"
		# 		for pdb in pdbid: 
		# 			pdbpath = basedir + pdb + "/" + pdb + "_tr_ignorechain.pdb\n"
		# 			f.write( pdbpath )

		# Run predict side chain distribution script##"/designed_best_score.list"
		
		if(not Options.energy_fxn=='proteinmpnn'):
			predict_side_chain_distribution.compute_side_chain_distribution( config, datadir + "/natives.list", datadir + "/design.list", 'overall')
   			# predict_side_chain_distribution.compute_side_chain_distribution( config, datadir + "/betanatives.list", datadir + "/beta_design.list", 'beta')
		else:
			# predict_side_chain_distribution.compute_side_chain_distribution( config, datadir + "/natives.list", datadir + "/design.list", 'overall', datadir + "/mpnn_seq_recov.dat")
			# predict_side_chain_distribution.compute_side_chain_distribution( config, datadir + "/betanatives.list", datadir + "/design.list", 'beta', datadir + "/betampnn_seq_recov.dat")
			predict_side_chain_distribution.generate_heatmap_from_file(datadir)
   
		# predict_side_chain_distribution.compute_side_chain_distribution( config, datadir + "/natives_short.list", datadir + "/design_short.list", 'overall')
		# predict_side_chain_distribution.compute_side_chain_distribution( config, datadir + "/natives_alpha.list", datadir + "/design_alpha.list", 'alpha')
		# predict_side_chain_distribution.compute_side_chain_distribution( config, datadir + "/natives_beta.list", datadir + "/design_beta.list", 'beta')
		predict_side_chain_distribution.generate_heatmap_from_file(datadir)
  		# predict_side_chain_distribution.compute_side_chain_distribution( config, datadir + "/../native_test.list", datadir + "/design_test.list" )

	if ( "decoy-discrimination" in test_names ): 

		# for each target, make a list of pdb decoys
		targets = [ "brd7", "fmr5", "ltpa", "rhod", "vatp"]
		datadir = config.benchmark_path + "data/" + Options.energy_fxn + "/decoy-discrimination/highdGforDE/"
		# first hires targets
		hires_dir = datadir + "highres"
		os.chdir( hires_dir )
		for target in targets: 
			os.chdir( target )
			os.system( "ls decoy*.pdb > decoys.list" )
			os.chdir( "../" )

		# then lowres targets
		lowres_dir = datadir + "lowres"
		os.chdir( lowres_dir )
		for target in targets: 
			os.chdir( target )
			os.system( "ls *.pdb > decoys.list" )
			os.chdir( "../" )

		# rescore each to calculate the rms and total score
		executable = config.rosetta_path + "/score_jd2" + "." + config.platform + config.compiler + config.buildenv
		s = Template( "-in:file:l decoys.list -in:file:native $native -mp:setup:spanfiles from_structure -out:file:scorefile $scorefile -in:membrane")
		for target in targets: 

			# hires
			targetdir = hires_dir + "/" + target
			os.chdir( targetdir )
			output_scores = target + "_hires.sc"
			native_pdb = config.benchmark_path + "targets/structure/D6_decoy_discrimination/highres/" + target + "/" + target + "_native.pdb"
			spanfile = config.benchmark_path + "targets/structure/D6_decoy_discrimination/highres/" + target + "/" + target + ".span"
			arguments = s.substitute( scorefile=output_scores, native=native_pdb )
			jobname = "rescore_hires_" + target
			hpc_util.submit_condor_job( targetdir, jobname, executable, arguments, 1 )

			# lowres
			targetdir = lowres_dir + "/" + target
			os.chdir( targetdir )
			output_scores = target + "_lowres.sc"
			native_pdb = config.benchmark_path + "targets/structure/D6_decoy_discrimination/lowres/" + target + "/" + target + "_native.pdb"
			spanfile = config.benchmark_path + "targets/structure/D6_decoy_discrimination/lowres/" + target + "/" + target + ".span"
			arguments = s.substitute( scorefile=output_scores, native=native_pdb )
			jobname = "rescore_lowres_" + target
			hpc_util.submit_condor_job( targetdir, jobname, executable, arguments, 1 )

	if ( "helix-kinks" in test_names ): 

		# Read list of test case IDs and PDBs
		targets_path = config.benchmark_path + "targets/structure/D5_helix_kinks/"
		list_of_targets = targets_path + "targets.list"
		with open( list_of_targets, 'rt' ) as f: 
			targets = f.readlines()
			targets = [ x.strip() for x in targets ]

		# Change directories to a data analysis dir
		outdir = config.benchmark_path + "data/" + Options.energy_fxn + "/helix-kinks/" 
		if ( not os.path.isdir( outdir ) ): 
			aya.exit( "Data generation for helix kinks was skipped! Cannot process non-existent data. Exiting!" )
		os.chdir( outdir )

		# Iterate through each target
		for target in targets:

			# Clean up the directories with refined models
			casedir = outdir + target 
			os.chdir( casedir )
			os.system( "rm *.in_progress" )
			os.system( "mkdir metadata" )
			os.system( "mv *.out metadata/" )
			os.system( "mv *.err metadata/" )
			os.system( "mv *.log metadata/" )
			os.system( "mkdir scorefiles" )
			os.system( "mv *.sc scorefiles/" )

			# Make a list of the refined models
			os.system( "ls *.pdb > models.list" )
			with open( "models.list" ) as f: 
				refined_models = f.readlines()
				refined_models = [ x.strip() for x in refined_models ]

			# Conver tthe spanfile from the targets directory to a string format
			kink_finder_script = config.benchmark_path + "external/kink-finder/Kink_Finder.py"
			spanfile = config.benchmark_path + "targets/structure/D5_helix_kinks/" + target + "/" + target + ".span"
			spanstr = convert_spanfile_to_string( spanfile )
			os.system( "mkdir kink_data" )
			i = 1
			for m in range(0, len(refined_models)): 
				kink_cmd = "python2 " + kink_finder_script + " -f " + refined_models[m] + " -l '" + spanstr + "'"
				os.system( kink_cmd )
				os.system( "mv Output/angles.csv kink_data/" + target + "_" + str(i) + ".csv" ) 
				i = i + 1

			# Run kink finder on each model
			# will need to re-format the helix definitions here and then run kink finder
			# on each model, output into a kink files directory

			# then run the process kink data script

	if ( "protein-protein-docking" in test_names ): 

		print("temp")
		# Analyze docking decoys
		make_asymm_docked_complexes.analyze_interfaces( Options.energy_fxn, config, "D2_single_TM_complexes", "protein-protein-docking", False, True )
		make_asymm_docked_complexes.analyze_interfaces( Options.energy_fxn, config, "D3_multi_TM_bound_complexes", "protein-protein-docking", False  )
		make_asymm_docked_complexes.analyze_interfaces( Options.energy_fxn, config, "D4_multi_TM_unbound_complexes", "protein-protein-docking", False  )

		# Analyze local refine decoys
		make_asymm_docked_complexes.analyze_interfaces( Options.energy_fxn, config, "D2_single_TM_complexes", "protein-protein-docking", True, True )
		make_asymm_docked_complexes.analyze_interfaces( Options.energy_fxn, config, "D3_multi_TM_bound_complexes", "protein-protein-docking", True  )
		make_asymm_docked_complexes.analyze_interfaces( Options.energy_fxn, config, "D4_multi_TM_unbound_complexes", "protein-protein-docking", True  )


if __name__ == "__main__" : main(sys.argv)

