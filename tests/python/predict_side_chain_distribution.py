#!/usr/bin/env python
###################################################################
#@file:         predict_side_chain_distributions.py                                
#@description:  Calculate the depth-dependent side chain 
#	distributions from designed PDBS also characterize the design calculations.                                                      
#@author: Rituparna Samanta
#@email: rituparna@utexas.edu                                          
#@author: Rebecca F. Alford                   
#@email: rfalford12@gmail.com
###################################################################

import random, os, sys
import hpc_util, read_config
from string import Template
from progress.bar import Bar
from pyrosetta import * 
import numpy as np
import pandas as pd
import shutil
import collections
import math
import seaborn as sns
import chart_studio.plotly as py
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from collections import OrderedDict
from pyrosetta.rosetta.core.chemical import aa_from_oneletter_code

init()

class AminoAcidInfo: 
	def __init__(self): 
		self.coords = []
		self.lipidfacing_sequence = np.zeros(20)
		self.interface_sequence = np.zeros(20)
		self.aqueous_sequence = np.zeros(20)
		self.lipid_count = 0
		self.interface_count = 0
		self.aqueous_count = 0
		self.nativeoverall_count = 0
		self.designlipid_count = 0
		self.designinterface_count = 0
		self.designaqueous_count = 0
		self.designoverall_count = 0
		

aa_label = [ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 
            'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' ]

def compute_side_chain_distribution( config, native_pdbs, designed_pdbs, suffix='',result_csv='' ): 
	"""
		A function for calculating the depth dependent amino acid distributions

		Arguments: 
			config = Path to benchmark and Rosetta executables
			native_pdbs = List of native PDBs
			designed_pdbs = List of designed PDBs
	"""

	# Read pdblist and write to src file for natives
	native_poses = read_pdbs( native_pdbs )
	get_side_chain_distribution( native_poses, "native" )

	# Read pdblist and write to src file for designed pdbs
	designed_poses = read_pdbs( designed_pdbs )
	get_side_chain_distribution( designed_poses, "design" )
	
	option_string = "-run:constant_seed -in:ignore_unrecognized_res"

	init( extra_options=option_string )
	
	if(result_csv!=''):
		# get_rotamer_distribution_from_design_seq( result_csv, native_pdbs, suffix)
		calculate_aminoacid_distribution_from_design_seq( result_csv, native_pdbs, suffix )
                
	else:
		get_rotamer_distribution( designed_pdbs, native_pdbs, suffix ) 
		calculate_aminoacid_distribution( designed_pdbs, native_pdbs, suffix )
		# calculate_hydropathy(designed_pdbs,suffix)
 
	# pore_hydration_distribution( designed_pdbs, native_pdbs )
	# pore_hydration_distribution_file()
	calculate_residue(native_pdbs)
	calculate_design_distribution( designed_pdbs )

def read_pdbs( pdbfile_list ): 

	with open( pdbfile_list, 'rt' ) as f: 
		contents = f.readlines()
		contents = [ x.strip() for x in contents ]
	f.close()

	pdbs = []
	for pdbfile in contents: 
		if ( os.path.isfile( pdbfile ) ): 
			pdb = pose_from_pdb( pdbfile )
			add_memb = pyrosetta.rosetta.protocols.membrane.AddMembraneMover( "from_structure" )
			add_memb.apply( pdb )
			pdbs.append( pdb )
		else: 
			sys.exit( "PDB File " + pdbfile + " not found!" )
	
	return pdbs

def get_side_chain_distribution( pdblist, src ): 

	print("Generating distributions for source " + src )

	# Initialize an AA info object for each amino acid
	ala = AminoAcidInfo()
	cys = AminoAcidInfo()
	asp = AminoAcidInfo()
	glu = AminoAcidInfo()
	phe = AminoAcidInfo()
	gly = AminoAcidInfo()
	his = AminoAcidInfo()
	ile = AminoAcidInfo()
	lys = AminoAcidInfo()
	leu = AminoAcidInfo()
	met = AminoAcidInfo()
	asn = AminoAcidInfo()
	pro = AminoAcidInfo()
	gln = AminoAcidInfo()
	arg = AminoAcidInfo()
	ser = AminoAcidInfo()
	thr = AminoAcidInfo()
	val = AminoAcidInfo()
	trp = AminoAcidInfo()
	tyr = AminoAcidInfo()

	# Make a progress bar
	bar = Bar('Reading amino acids in PDB files', max=len(pdblist))

	# Loop through each PDB and encode case statement
	for pdb in pdblist: 
		for i in range(1, pdb.total_residue()+1): 

			# Get the hydration value for each side chain
			hyd = pdb.conformation().membrane_info().implicit_lipids().f_hydration( pdb.residue(i).xyz("CA") )
			if ( hyd > 0.5 ):

				if ( pdb.residue(i).name1() == "A" ): 
					ala.coords.append( pdb.residue(i).xyz("CA").z )
				elif ( pdb.residue(i).name1() == "C" ): 
					cys.coords.append( pdb.residue(i).xyz("CA").z )
				elif ( pdb.residue(i).name1() == "D" ): 
					asp.coords.append( pdb.residue(i).xyz("CA").z )
				elif ( pdb.residue(i).name1() == "E" ): 
					glu.coords.append( pdb.residue(i).xyz("CA").z )
				elif ( pdb.residue(i).name1() == "F" ): 
					phe.coords.append( pdb.residue(i).xyz("CA").z )
				elif ( pdb.residue(i).name1() == "G" ): 
					gly.coords.append( pdb.residue(i).xyz("CA").z )
				elif ( pdb.residue(i).name1() == "H" ): 
					his.coords.append( pdb.residue(i).xyz("CA").z )
				elif ( pdb.residue(i).name1() == "I" ): 
					ile.coords.append( pdb.residue(i).xyz("CA").z )
				elif ( pdb.residue(i).name1() == "K" ): 
					lys.coords.append( pdb.residue(i).xyz("CA").z )
				elif ( pdb.residue(i).name1() == "L" ): 
					leu.coords.append( pdb.residue(i).xyz("CA").z )
				elif ( pdb.residue(i).name1() == "M" ):
					met.coords.append( pdb.residue(i).xyz("CA").z )
				elif ( pdb.residue(i).name1() == "N" ): 
					asn.coords.append( pdb.residue(i).xyz("CA").z )
				elif ( pdb.residue(i).name1() == "P" ): 
					pro.coords.append( pdb.residue(i).xyz("CA").z )
				elif ( pdb.residue(i).name1() == "Q" ): 
					gln.coords.append( pdb.residue(i).xyz("CA").z )
				elif ( pdb.residue(i).name1() == "R" ): 
					arg.coords.append( pdb.residue(i).xyz("CA").z )
				elif ( pdb.residue(i).name1() == "S" ): 
					ser.coords.append( pdb.residue(i).xyz("CA").z )
				elif ( pdb.residue(i).name1() == "T" ): 
					thr.coords.append( pdb.residue(i).xyz("CA").z )
				elif ( pdb.residue(i).name1() == "V" ): 
					val.coords.append( pdb.residue(i).xyz("CA").z )
				elif ( pdb.residue(i).name1() == "W" ): 
					trp.coords.append( pdb.residue(i).xyz("CA").z )
				elif ( pdb.residue(i).name1() == "Y" ): 
					tyr.coords.append( pdb.residue(i).xyz("CA").z )

		bar.next()
	bar.finish()

	# Write each of the distributions to an output data file... 
	filename = src + "_side_chain_distribution.txt"
	with open( filename, 'wt' ) as f: 

		# write file header
		f.write( "AA zcoord src" )

		# Write for each coordinate list
		aa = "A"
		print("Writing coordinates for alanine...")
		for coord in ala.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "C"
		print("Writing coordinates for cystine...")
		for coord in cys.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "D"
		print("Writing coordinates for aspartate...")
		for coord in asp.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "E"
		print("Writing coordinates for glutamate...")
		for coord in glu.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "F"
		print("Writing coordinates for phenylaalanine...")
		for coord in phe.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "G"
		print("Writing coordinates for glycine...")
		for coord in gly.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "H"
		print("Writing coordinates for histidine...")
		for coord in his.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "I"
		print("Writing coordinates for isoleucine...")
		for coord in ile.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "K"
		print("Writing coordinates for lysine...")
		for coord in lys.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "L"
		print("Writing coordinates for leucine...")
		for coord in leu.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "M"
		print("Writing coordinates for methionine...")
		for coord in met.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "N"
		print("Writing coordinates for asparagine...")
		for coord in asn.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "P"
		print("Writing coordinates for proline...")
		for coord in pro.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "Q"
		print("Writing coordinates for glutamate...")
		for coord in gln.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "R"
		print("Writing coordinates for arginine...")
		for coord in arg.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "S"
		print("Writing coordinates for serine...")
		for coord in ser.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "T"
		print("Writing coordinates for threonine...")
		for coord in thr.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "V"
		print("Writing coordinates for valine...")
		for coord in val.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "W"
		print("Writing coordinates for tryptophan...")
		for coord in trp.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "Y"
		print("Writing coordinates for tyrosine...")
		for coord in tyr.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		f.close()

def get_rotamer_distribution( pdblist, nativelist, suffix): 
	"""
		A function for calculating the depth dependent amino acid design distributions

		Arguments: 
			nativelist = List of native PDBs
			designedlist = List of designed PDBs
		Output:
			heatmap of design amino acids for each native amino acids.  
	"""

	print("Generating distributions for source " + "native and design" )

	# Initialize an AA info object for each amino acid
	ala = AminoAcidInfo()
	cys = AminoAcidInfo()
	asp = AminoAcidInfo()
	glu = AminoAcidInfo()
	phe = AminoAcidInfo()
	gly = AminoAcidInfo()
	his = AminoAcidInfo()
	ile = AminoAcidInfo()
	lys = AminoAcidInfo()
	leu = AminoAcidInfo()
	met = AminoAcidInfo()
	asn = AminoAcidInfo()
	pro = AminoAcidInfo()
	gln = AminoAcidInfo()
	arg = AminoAcidInfo()
	ser = AminoAcidInfo()
	thr = AminoAcidInfo()
	val = AminoAcidInfo()
	trp = AminoAcidInfo()
	tyr = AminoAcidInfo()

	with open( pdblist, 'rt' ) as f: 
			contents_pdb = f.readlines()
			contents_pdb = [ x.strip() for x in contents_pdb ]
	f.close()

	with open( nativelist, 'rt' ) as f: 
			contents_native = f.readlines()
			contents_native = [ x.strip() for x in contents_native ]
	f.close()

	# pdbs = []
	# for pdbfile in contents: 
	# 	if ( os.path.isfile( pdbfile ) ): 
	# 		pdb = pose_from_pdb( pdbfile )
	# 		add_memb = pyrosetta.rosetta.protocols.membrane.AddMembraneMover( "from_structure" )
	# 		add_memb.apply( pdb )

	# Make a progress bar
	bar = Bar('Reading amino acids in PDB files \n', max=len(contents_pdb))

	aa_name = np.array(['ala','cys','asp','glu','phe','gly','his','ile','lys','leu','met','asn','pro','gln','arg','ser','thr','val','trp','tyr'])
	aa_list = np.arange(1,21)

	# Loop through each PDB and encode case statement
	for pose_num in range(len(contents_pdb)): 
		pdb =[]
		native = []

		pdb = pose_from_pdb(contents_pdb[pose_num])
		native = pose_from_pdb(contents_native[pose_num])
		add_memb = pyrosetta.rosetta.protocols.membrane.AddMembraneMover( "from_structure" )
		# add_memb.apply( pdb )
		# print('mem atom added to pdb')
		# add_memb = pyrosetta.rosetta.protocols.membrane.AddMembraneMover( "from_structure" )
		add_memb.apply(native)
		# print('mem atom added to native')
		# print('native length:'+str(native.total_residue()))
		# print('pdb length:'+str(pdb.total_residue()))

		if(native.total_residue() == pdb.total_residue()):
			#loops over residues except the mem atom
			for i in range(1,native.total_residue()): 
				# print(i)
				# Get the hydration value for each side chain
				hyd = native.conformation().membrane_info().implicit_lipids().f_hydration( pdb.residue(i).xyz("CA") )
				thk = native.conformation().membrane_info().implicit_lipids().f_depth( pdb.residue(i).xyz("CA")[2] )
				if(native.conformation().membrane_info().implicit_lipids().has_pore()):
					cavity = native.conformation().membrane_info().implicit_lipids().f_cavity( pdb.residue(i).xyz("CA") )
				else:
					cavity = 0.0
				
			
				# print(hyd) 
				
				#aa_from_oneletter_code(pdb.residue(i).name1()) is from 1-20
				#the residue numbering is from 1-
				#the array numbering is from 0 though. This needs to be fixed. 
				#is there any other sum check we can do? 

				#the design and naive poses align
				#find the index
				index = np.where(aa_name == native.residue(i).name3().lower())
				# print(native.residue(i).name3().lower())
				# print(index)
				index = index[0][0]
				# print("native residue: "+native.residue(i).name3().lower())
				
				# print("aa name is: " + aa_name[index])
				# sys.exit()
				
				if(thk<0.25 and cavity<0.50):
					# print("loop 1: thk = {}, cavity = {}, hyd = {}".format(thk,cavity,hyd))
					eval(aa_name[index]).lipidfacing_sequence[int(aa_from_oneletter_code(pdb.residue(i).name1()))-1] = eval(aa_name[index]).lipidfacing_sequence[int(aa_from_oneletter_code(pdb.residue(i).name1()))-1]+1 
					eval(aa_name[index]).lipid_count = eval(aa_name[index]).lipid_count+1

				elif(thk>0.25 and thk<0.75 and cavity<0.50):
					# print("loop 2: thk = {}, cavity = {}, hyd = {}".format(thk,cavity,hyd))
					eval(aa_name[index]).interface_sequence[int(aa_from_oneletter_code(pdb.residue(i).name1()))-1] = eval(aa_name[index]).interface_sequence[int(aa_from_oneletter_code(pdb.residue(i).name1()))-1]+1 
					eval(aa_name[index]).interface_count = eval(aa_name[index]).interface_count+1

				elif(thk>0.75 or cavity>0.50):
					# print("loop 3: thk = {}, cavity = {}, hyd = {}".format(thk,cavity,hyd))		
					eval(aa_name[index]).aqueous_sequence[int(aa_from_oneletter_code(pdb.residue(i).name1()))-1] = eval(aa_name[index]).aqueous_sequence[int(aa_from_oneletter_code(pdb.residue(i).name1()))-1]+1 
					eval(aa_name[index]).aqueous_count = eval(aa_name[index]).aqueous_count+1
				
				
				# if ( native.residue(i).name1() == "A" ):
					
				# 	# print(native.residue(i).name1())
				# 	# print(int(aa_from_oneletter_code(native.residue(i).name1()))) 
				# 	if(hyd<0.25):
				# 		# print('lipid_facing')
				# 		ala.lipidfacing_sequence[int(aa_from_oneletter_code(pdb.residue(i).name1()))-1] = ala.lipidfacing_sequence[int(aa_from_oneletter_code(pdb.residue(i).name1()))-1]+1 
				# 		ala.lipid_count = ala.lipid_count+1
				# 	elif( (hyd >= 0.25) and (hyd <= 0.75) ):
				# 		# print('interface')
				# 		ala.interface_sequence[int(aa_from_oneletter_code(pdb.residue(i).name1()))-1] = ala.interface_sequence[int(aa_from_oneletter_code(pdb.residue(i).name1()))-1]+1 
				# 		ala.interface_count = ala.interface_count+1
				# 	else:
				# 		# print('aqeous')
				# 		ala.aqueous_sequence[int(aa_from_oneletter_code(pdb.residue(i).name1()))-1] = ala.aqueous_sequence[int(aa_from_oneletter_code(pdb.residue(i).name1()))-1]+1 
				# 		ala.aqueous_count = ala.aqueous_count+1
		
				
			# print('lipid_count:'+' '+ str(ala.lipid_count+cys.lipid_count+asp.lipid_count+glu.lipid_count+phe.lipid_count+gly.lipid_count+his.lipid_count+ile.lipid_count+lys.lipid_count+leu.lipid_count+met.lipid_count+asn.lipid_count+pro.lipid_count+gln.lipid_count+arg.lipid_count+ser.lipid_count+thr.lipid_count+val.lipid_count+trp.lipid_count+tyr.lipid_count))
			# print('interface_count:'+' '+ str(ala.interface_count+cys.interface_count+asp.interface_count+glu.interface_count+phe.interface_count+gly.interface_count+his.interface_count+ile.interface_count+lys.interface_count+leu.interface_count+met.interface_count+asn.interface_count+pro.interface_count+gln.interface_count+arg.interface_count+ser.interface_count+thr.interface_count+val.interface_count+trp.interface_count+tyr.interface_count))
			# print('aqueous_count:'+' '+ str(ala.aqueous_count+cys.aqueous_count+asp.aqueous_count+glu.aqueous_count+phe.aqueous_count+gly.aqueous_count+his.aqueous_count+ile.aqueous_count+lys.aqueous_count+leu.aqueous_count+met.aqueous_count+asn.aqueous_count+pro.aqueous_count+gln.aqueous_count+arg.aqueous_count+ser.aqueous_count+thr.aqueous_count+val.aqueous_count+trp.aqueous_count+tyr.aqueous_count))
			# print('total:'+' '+str(native.total_residue()))


		else:
			print(contents_pdb[pose_num])
			print(contents_native[pose_num])
			print("check the native and design path list!")
			sys.exit()
		bar.next()
	bar.finish()

	
	lipid_facing_sr = []
	interface_sr = []
	aqueous_sr = []
	for aa in aa_list:
		
		# print(eval(aa_name[aa-1]).lipidfacing_sequence/len(pdblist))
		# break
		if(aa == 1.0):
			lipid_facing_sr = eval(aa_name[aa-1]).lipidfacing_sequence/(eval(aa_name[aa-1]).lipid_count)
			interface_sr = eval(aa_name[aa-1]).interface_sequence/(eval(aa_name[aa-1]).interface_count)
			aqueous_sr = eval(aa_name[aa-1]).aqueous_sequence/(eval(aa_name[aa-1]).aqueous_count)
		else:
			lipid_facing_sr = np.vstack((lipid_facing_sr,eval(aa_name[aa-1]).lipidfacing_sequence/eval(aa_name[aa-1]).lipid_count))
			interface_sr = np.vstack((interface_sr,eval(aa_name[aa-1]).interface_sequence/eval(aa_name[aa-1]).interface_count))
			aqueous_sr = np.vstack((aqueous_sr,eval(aa_name[aa-1]).aqueous_sequence/eval(aa_name[aa-1]).aqueous_count))
	
	lipid_df = pd.DataFrame(lipid_facing_sr.round(3), columns=aa_name, index = aa_name )
	interface_df = pd.DataFrame(interface_sr.round(3), columns=aa_name, index = aa_name)
	aqueous_df = pd.DataFrame(aqueous_sr.round(3), columns=aa_name, index = aa_name)
	
	print(lipid_df)
	print('----reorder------')
 	#reordering columns in the dataframes
	#non-polar: A,I,L,M,V
	#Polar: N,Q,S,T
	#+ charge: H,K,R
	#-charge: D,E
 	#aromatic:F,Y,W
	#special: C,G,P
	residue_type = ["ala","ile","leu","met","val","asn","gln","ser","thr","his","lys","arg","asp","glu", "phe","trp","tyr","cys","gly","pro"]
	residue_type_rev = residue_type[::-1]
	lipid_df = lipid_df[["ala","ile","leu","met","val","asn","gln","ser","thr","his","lys","arg","asp","glu", "phe","trp","tyr","cys","gly","pro"]]
	interface_df = interface_df[["ala","ile","leu","met","val","asn","gln","ser","thr","his","lys","arg","asp","glu", "phe","trp","tyr","cys","gly","pro"]]
	aqueous_df = aqueous_df[["ala","ile","leu","met","val","asn","gln","ser","thr","his","lys","arg","asp","glu", "phe","trp","tyr","cys","gly","pro"]]
	
	lipid_df = lipid_df.loc[residue_type]
	interface_df = interface_df.loc[residue_type]
	aqueous_df = aqueous_df.loc[residue_type]
	
	print(lipid_df)
	print('------------------Transpose-------------------')
	lipid_df_mod = lipid_df.T
	interface_df_mod = interface_df.T
	aqueous_df_mod = aqueous_df.T

	lipid_df_mod = lipid_df_mod.loc[residue_type_rev]
	interface_df_mod = interface_df_mod.loc[residue_type_rev]
	aqueous_df_mod = aqueous_df_mod.loc[residue_type_rev]
 
	print(lipid_df_mod)
	print(interface_df_mod)
	print(aqueous_df_mod)
 
	lipid_df.to_csv('output_lipid_sr_heatmap_'+suffix+'.txt', sep=" ", index=True)
	interface_df.to_csv('output_interface_sr_heatmap_'+suffix+'.txt', sep=" ", index=True)
	aqueous_df.to_csv('output_aqueous_sr_heatmap_'+suffix+'.txt', sep=" ", index=True)
	
	get_heatmap(lipid_df_mod,'lipid_heatmap_'+suffix)
	get_heatmap(interface_df_mod,'interface_heatmap_'+suffix)
	get_heatmap(aqueous_df_mod, 'aqueous_heatmap_'+suffix)
	
 	# print(lipid_facing_sr)
	# print(interface_sr)
	# print(aqueous_sr)
 
def get_rotamer_distribution_from_design_seq( result_csv, nativelist, suffix): 
	"""
		A function for calculating the depth dependent amino acid design distributions

		Arguments: 
			nativelist = List of native PDBs
			results_csv = This is the csv file with pdb name, sequence recovery and maximum recovery sequence
							generated from proteinmpnn
		Assumption: 
			the design sequence folds into the native backbone. We have not tested this assumption yet. 
  		Output:
			heatmap of design amino acids for each native amino acids.  
	"""

	print("Generating distributions for source " + "native and design" )

	# Initialize an AA info object for each amino acid
	ala = AminoAcidInfo()
	cys = AminoAcidInfo()
	asp = AminoAcidInfo()
	glu = AminoAcidInfo()
	phe = AminoAcidInfo()
	gly = AminoAcidInfo()
	his = AminoAcidInfo()
	ile = AminoAcidInfo()
	lys = AminoAcidInfo()
	leu = AminoAcidInfo()
	met = AminoAcidInfo()
	asn = AminoAcidInfo()
	pro = AminoAcidInfo()
	gln = AminoAcidInfo()
	arg = AminoAcidInfo()
	ser = AminoAcidInfo()
	thr = AminoAcidInfo()
	val = AminoAcidInfo()
	trp = AminoAcidInfo()
	tyr = AminoAcidInfo()

	results_mpnn = pd.read_csv(result_csv,sep=' ')
	# print(results_mpnn)
	# sys.exit()
	with open( nativelist, 'rt' ) as f: 
			contents_native = f.readlines()
			contents_native = [ x.strip() for x in contents_native ]
	f.close()

	# pdbs = []
	# for pdbfile in contents: 
	# 	if ( os.path.isfile( pdbfile ) ): 
	# 		pdb = pose_from_pdb( pdbfile )
	# 		add_memb = pyrosetta.rosetta.protocols.membrane.AddMembraneMover( "from_structure" )
	# 		add_memb.apply( pdb )

	# Make a progress bar
	bar = Bar('Reading amino acids in PDB files \n', max=len(contents_native))

	aa_name = np.array(['ala','cys','asp','glu','phe','gly','his','ile','lys','leu','met','asn','pro','gln','arg','ser','thr','val','trp','tyr'])
	aa_list = np.arange(1,21)

	# Loop through each PDB and encode case statement
	for pose_num in range(len(contents_native)): 
		
		native = []

		# pdb = pose_from_pdb(contents_pdb[pose_num])
		native = pose_from_pdb(contents_native[pose_num])
		pdb_name = ((contents_native[pose_num].split('/')[-1]).split('.')[0]).split('_')[0]
		print(pdb_name)
		add_memb = pyrosetta.rosetta.protocols.membrane.AddMembraneMover( "from_structure" )
		# add_memb.apply( pdb )
		# print('mem atom added to pdb')
		# add_memb = pyrosetta.rosetta.protocols.membrane.AddMembraneMover( "from_structure" )
		add_memb.apply(native)
		# print('mem atom added to native')
		# print('native length:'+str(native.total_residue()))
		# print('pdb length:'+str(pdb.total_residue()))
		design_seq = results_mpnn[results_mpnn['pdb']==pdb_name]['max_seq'].item()
		print(design_seq)
		if(design_seq==' '):
			continue
		if(native.total_residue()-1 == len(design_seq)):
			#loops over residues except the mem atom
			for i in range(1,native.total_residue()): 
				# print(i)
				# Get the hydration value for each side chain
				hyd = native.conformation().membrane_info().implicit_lipids().f_hydration( native.residue(i).xyz("CA") )
				thk = native.conformation().membrane_info().implicit_lipids().f_depth( native.residue(i).xyz("CA")[2] )
				if(native.conformation().membrane_info().implicit_lipids().has_pore()):
					cavity = native.conformation().membrane_info().implicit_lipids().f_cavity( native.residue(i).xyz("CA") )
				else:
					cavity = 0.0
				
			
				
				#the design and naive poses align
				#find the index
				index = np.where(aa_name == native.residue(i).name3().lower())
				# print(native.residue(i).name3().lower())
				# print(index)
				index = index[0][0]
				# print("native residue: "+native.residue(i).name3().lower())
				
				# print("aa name is: " + aa_name[index])
				# sys.exit()
				
				if(thk<0.25 and cavity<0.50):
					# print("loop 1: thk = {}, cavity = {}, hyd = {}".format(thk,cavity,hyd))
					eval(aa_name[index]).lipidfacing_sequence[int(aa_from_oneletter_code(design_seq[i-1]))-1] = eval(aa_name[index]).lipidfacing_sequence[int(aa_from_oneletter_code(design_seq[i-1]))-1]+1 
					eval(aa_name[index]).lipid_count = eval(aa_name[index]).lipid_count+1

				elif(thk>0.25 and thk<0.75 and cavity<0.50):
					# print("loop 2: thk = {}, cavity = {}, hyd = {}".format(thk,cavity,hyd))
					eval(aa_name[index]).interface_sequence[int(aa_from_oneletter_code(design_seq[i-1]))-1] = eval(aa_name[index]).interface_sequence[int(aa_from_oneletter_code(design_seq[i-1]))-1]+1 
					eval(aa_name[index]).interface_count = eval(aa_name[index]).interface_count+1

				elif(thk>0.75 or cavity>0.50):
					# print("loop 3: thk = {}, cavity = {}, hyd = {}".format(thk,cavity,hyd))		
					eval(aa_name[index]).aqueous_sequence[int(aa_from_oneletter_code(design_seq[i-1]))-1] = eval(aa_name[index]).aqueous_sequence[int(aa_from_oneletter_code(design_seq[i-1]))-1]+1 
					eval(aa_name[index]).aqueous_count = eval(aa_name[index]).aqueous_count+1
				
				
				
			# print('lipid_count:'+' '+ str(ala.lipid_count+cys.lipid_count+asp.lipid_count+glu.lipid_count+phe.lipid_count+gly.lipid_count+his.lipid_count+ile.lipid_count+lys.lipid_count+leu.lipid_count+met.lipid_count+asn.lipid_count+pro.lipid_count+gln.lipid_count+arg.lipid_count+ser.lipid_count+thr.lipid_count+val.lipid_count+trp.lipid_count+tyr.lipid_count))
			# print('interface_count:'+' '+ str(ala.interface_count+cys.interface_count+asp.interface_count+glu.interface_count+phe.interface_count+gly.interface_count+his.interface_count+ile.interface_count+lys.interface_count+leu.interface_count+met.interface_count+asn.interface_count+pro.interface_count+gln.interface_count+arg.interface_count+ser.interface_count+thr.interface_count+val.interface_count+trp.interface_count+tyr.interface_count))
			# print('aqueous_count:'+' '+ str(ala.aqueous_count+cys.aqueous_count+asp.aqueous_count+glu.aqueous_count+phe.aqueous_count+gly.aqueous_count+his.aqueous_count+ile.aqueous_count+lys.aqueous_count+leu.aqueous_count+met.aqueous_count+asn.aqueous_count+pro.aqueous_count+gln.aqueous_count+arg.aqueous_count+ser.aqueous_count+thr.aqueous_count+val.aqueous_count+trp.aqueous_count+tyr.aqueous_count))
			# print('total:'+' '+str(native.total_residue()))


		else:
			# print(contents_pdb[pose_num])
			print(contents_native[pose_num])
			print("check the native and design path list!")
			sys.exit()
		bar.next()
	bar.finish()

	
	lipid_facing_sr = []
	interface_sr = []
	aqueous_sr = []
	for aa in aa_list:
		
		# print(eval(aa_name[aa-1]).lipidfacing_sequence/len(pdblist))
		# break
		if(aa == 1.0):
			lipid_facing_sr = eval(aa_name[aa-1]).lipidfacing_sequence/(eval(aa_name[aa-1]).lipid_count)
			interface_sr = eval(aa_name[aa-1]).interface_sequence/(eval(aa_name[aa-1]).interface_count)
			aqueous_sr = eval(aa_name[aa-1]).aqueous_sequence/(eval(aa_name[aa-1]).aqueous_count)
		else:
			lipid_facing_sr = np.vstack((lipid_facing_sr,eval(aa_name[aa-1]).lipidfacing_sequence/eval(aa_name[aa-1]).lipid_count))
			interface_sr = np.vstack((interface_sr,eval(aa_name[aa-1]).interface_sequence/eval(aa_name[aa-1]).interface_count))
			aqueous_sr = np.vstack((aqueous_sr,eval(aa_name[aa-1]).aqueous_sequence/eval(aa_name[aa-1]).aqueous_count))
	
	lipid_df = pd.DataFrame(lipid_facing_sr.round(3), columns=aa_name, index = aa_name )
	interface_df = pd.DataFrame(interface_sr.round(3), columns=aa_name, index = aa_name)
	aqueous_df = pd.DataFrame(aqueous_sr.round(3), columns=aa_name, index = aa_name)
	
	print(lipid_df)
	print('----reorder------')
 	#reordering columns in the dataframes
	#non-polar: A,I,L,M,V
	#Polar: N,Q,S,T
	#+ charge: H,K,R
	#-charge: D,E
 	#aromatic:F,Y,W
	#special: C,G,P
	residue_type = ["ala","ile","leu","met","val","asn","gln","ser","thr","his","lys","arg","asp","glu", "phe","trp","tyr","cys","gly","pro"]
	residue_type_rev = residue_type[::-1]
	lipid_df = lipid_df[["ala","ile","leu","met","val","asn","gln","ser","thr","his","lys","arg","asp","glu", "phe","trp","tyr","cys","gly","pro"]]
	interface_df = interface_df[["ala","ile","leu","met","val","asn","gln","ser","thr","his","lys","arg","asp","glu", "phe","trp","tyr","cys","gly","pro"]]
	aqueous_df = aqueous_df[["ala","ile","leu","met","val","asn","gln","ser","thr","his","lys","arg","asp","glu", "phe","trp","tyr","cys","gly","pro"]]
	
	lipid_df = lipid_df.loc[residue_type]
	interface_df = interface_df.loc[residue_type]
	aqueous_df = aqueous_df.loc[residue_type]
	
	print(lipid_df)
	print('------------------Transpose-------------------')
	lipid_df_mod = lipid_df.T
	interface_df_mod = interface_df.T
	aqueous_df_mod = aqueous_df.T

	lipid_df_mod = lipid_df_mod.loc[residue_type_rev]
	interface_df_mod = interface_df_mod.loc[residue_type_rev]
	aqueous_df_mod = aqueous_df_mod.loc[residue_type_rev]
 
	print(lipid_df_mod)
	print(interface_df_mod)
	print(aqueous_df_mod)
 
	lipid_df.to_csv('output_lipid_sr_heatmap_'+suffix+'.txt', sep=" ", index=True)
	interface_df.to_csv('output_interface_sr_heatmap_'+suffix+'.txt', sep=" ", index=True)
	aqueous_df.to_csv('output_aqueous_sr_heatmap_'+suffix+'.txt', sep=" ", index=True)
	
	get_heatmap(lipid_df_mod,'lipid_heatmap_'+suffix)
	get_heatmap(interface_df_mod,'interface_heatmap_'+suffix)
	get_heatmap(aqueous_df_mod, 'aqueous_heatmap_'+suffix)
	
 	# print(lipid_facing_sr)
	# print(interface_sr)
	# print(aqueous_sr)

def calculate_aminoacid_distribution( pdblist, nativelist, suffix ): 
	"""
		A function for calculating amino acid design distributions

		Arguments: 
			nativelist = List of native PDBs
			designedlist = List of designed PDBs
		Output:
			a frequency count relative to position in the membrane.   
	"""

	print("Generating distributions for source " + "native and design" )
	import seaborn as sns
	import matplotlib
	import matplotlib.pyplot as plt
	amino_acids = {'A':'ala','C':'cys','D':'asp','E':'glu',
				'F':'phe','G':'gly','H':'his','I':'ile', 
				'K':'lys','L':'leu','M':'met','N':'asn',
				'P':'pro','Q':'gln','R':'arg',
				'S':'ser','T':'thr','V':'val','W':'trp','Y':'tyr'}

	# # Initialize an AA info object for each amino acid
	# ala = AminoAcidInfo()
	# cys = AminoAcidInfo()
	# asp = AminoAcidInfo()
	# glu = AminoAcidInfo()
	# phe = AminoAcidInfo()
	# gly = AminoAcidInfo()
	# his = AminoAcidInfo()
	# ile = AminoAcidInfo()
	# lys = AminoAcidInfo()
	# leu = AminoAcidInfo()
	# met = AminoAcidInfo()
	# asn = AminoAcidInfo()
	# pro = AminoAcidInfo()
	# gln = AminoAcidInfo()
	# arg = AminoAcidInfo()
	# ser = AminoAcidInfo()
	# thr = AminoAcidInfo()
	# val = AminoAcidInfo()
	# trp = AminoAcidInfo()
	# tyr = AminoAcidInfo()

	# with open( pdblist, 'rt' ) as f: 
	# 		contents_pdb = f.readlines()
	# 		contents_pdb = [ x.strip() for x in contents_pdb ]
	# f.close()

	# with open( nativelist, 'rt' ) as f: 
	# 		contents_native = f.readlines()
	# 		contents_native = [ x.strip() for x in contents_native ]
	# f.close()

	# # pdbs = []
	# # for pdbfile in contents: 
	# # 	if ( os.path.isfile( pdbfile ) ): 
	# # 		pdb = pose_from_pdb( pdbfile )
	# # 		add_memb = pyrosetta.rosetta.protocols.membrane.AddMembraneMover( "from_structure" )
	# # 		add_memb.apply( pdb )

	# # Make a progress bar
	# bar = Bar('Reading amino acids in PDB files \n', max=len(contents_pdb))

	# aa_name = np.array(['ala','cys','asp','glu','phe','gly','his','ile','lys','leu','met','asn','pro','gln','arg','ser','thr','val','trp','tyr'])
	# aa_list = np.arange(1,21)
	aa_name_letter = np.array(['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'])
	aa_by_group = ['A','I','L','M','V','N','Q','S','T','H','K','R','D','E','F','Y','W','C','G','P']


	# # Loop through each PDB and encode case statement
	# for pose_num in range(len(contents_pdb)): 
	# 	pdb =[]
	# 	native = []

	# 	pdb = pose_from_pdb(contents_pdb[pose_num])
	# 	native = pose_from_pdb(contents_native[pose_num])
	# 	add_memb = pyrosetta.rosetta.protocols.membrane.AddMembraneMover( "from_structure" )
	# 	# add_memb.apply( pdb )
	# 	# print('mem atom added to pdb')
	# 	# add_memb = pyrosetta.rosetta.protocols.membrane.AddMembraneMover( "from_structure" )
	# 	add_memb.apply(native)
	# 	# print('mem atom added to native')
	# 	# print('native length:'+str(native.total_residue()))
	# 	# print('pdb length:'+str(pdb.total_residue()))

	# 	if(native.total_residue() == pdb.total_residue()):
	# 		#loops over residues except the mem atom
	# 		for i in range(1,native.total_residue()): 
	# 			# print(i)
	# 			# Get the hydration value for each side chain
	# 			hyd = native.conformation().membrane_info().implicit_lipids().f_hydration( pdb.residue(i).xyz("CA") )
	# 			thk = native.conformation().membrane_info().implicit_lipids().f_depth( pdb.residue(i).xyz("CA")[2] )
	# 			if(native.conformation().membrane_info().implicit_lipids().has_pore()):
	# 				cavity = native.conformation().membrane_info().implicit_lipids().f_cavity( pdb.residue(i).xyz("CA") )
	# 			else:
	# 				cavity = 0.0
				
	# 			#the design and native poses align
	# 			#find the index
	# 			index = np.where(aa_name == native.residue(i).name3().lower())
	# 			index = index[0][0]
	# 			index_design = np.where(aa_name==pdb.residue(i).name3().lower())
	# 			index_design = index_design[0][0]
    
	# 			if(thk<0.25 and cavity<0.50):
	# 				eval(aa_name[index]).lipid_count = eval(aa_name[index]).lipid_count+1
	# 				eval(aa_name[index_design]).designlipid_count = eval(aa_name[index_design]).designlipid_count + 1

	# 			elif(thk>0.25 and thk<0.75 and cavity<0.50):
	# 				eval(aa_name[index]).interface_count = eval(aa_name[index]).interface_count+1
	# 				eval(aa_name[index_design]).designinterface_count = eval(aa_name[index_design]).designinterface_count+1
	# 			elif(thk>0.75 or cavity>0.50):
							
	# 				eval(aa_name[index]).aqueous_count = eval(aa_name[index]).aqueous_count+1
	# 				eval(aa_name[index_design]).designaqueous_count = eval(aa_name[index_design]).designaqueous_count+1
				

	# 	else:
	# 		print(contents_pdb[pose_num])
	# 		print(contents_native[pose_num])
	# 		print("check the native and design path list!")
	# 		sys.exit()
	# 	bar.next()
	# bar.finish()

	
	# final=[]
 
	# for aa in aa_list:
		
	# 	# print(eval(aa_name[aa-1]).lipidfacing_sequence/len(pdblist))
	# 	# break
		
	# 	lipid_facing = (eval(aa_name[aa-1]).lipid_count)
	# 	interface = (eval(aa_name[aa-1]).interface_count)
	# 	aqueous = (eval(aa_name[aa-1]).aqueous_count)
	# 	design_lipid_facing = (eval(aa_name[aa-1]).designlipid_count)
	# 	design_interface = (eval(aa_name[aa-1]).designinterface_count)
	# 	design_aqueous = (eval(aa_name[aa-1]).designaqueous_count)
	# 	final.append(np.stack([aa_name_letter[aa-1], lipid_facing, design_lipid_facing, interface, design_interface, aqueous, design_aqueous]))
	# # print(final)
 		
	# df_final = pd.DataFrame(final,columns=['AA','native_lipid','design_lipid','native_interface','design_interface','aqueous', 'design_aqueous'])
	# df_final.to_csv('output_protein_design_distribution.dat', sep = " ", index=False)
       
	# print(df_final)
	# print('----reorder------')
 	##============================================
	
	df_final = pd.read_csv('output_protein_design_distribution.dat', sep = " ")
	print(df_final)
	print(df_final['native_lipid'].sum())
	df_final['p_native_lipid'] = round(df_final['native_lipid']/(df_final['native_lipid'].sum()),3)
	df_final['p_design_lipid'] = round(df_final['design_lipid']/(df_final['design_lipid'].sum()),3)
	df_final['p_native_interface'] = round(df_final['native_interface']/(df_final['native_interface'].sum()),3)
	df_final['p_design_interface'] = round(df_final['design_interface']/(df_final['design_interface'].sum()),3)
	df_final['p_native_aqueous'] = round(df_final['aqueous']/(df_final['aqueous'].sum()),3)
	df_final['p_design_aqueous'] = round(df_final['design_aqueous']/(df_final['design_aqueous'].sum()),3)
	df_final.to_csv('output_protein_design_distribution.dat', sep = " ", index=False)
	df_new_native = df_final[['AA','p_native_lipid','p_native_interface','p_native_aqueous']]
	df_new_native.columns=['AA','lipid','interface','aqueous']
	df_new_native['type'] = 'native'

	df_new_design = df_final[['AA','p_design_lipid','p_design_interface','p_design_aqueous']]
	df_new_design.columns=['AA','lipid','interface','aqueous']
	df_new_design['type'] = 'design'


	df_new=pd.concat([df_new_native,df_new_design])
	print(df_new)
	# print('----reorder------')

	theme = {'axes.grid': True,
		'grid.linestyle': '',
		'xtick.labelsize': 18,
		'ytick.labelsize': 18,
		"font.weight": 'regular',
		'xtick.color': 'black',
		'ytick.color': 'black',
		"axes.titlesize": 20,
		"axes.labelsize": 18
	}

	for case in ['lipid','interface','aqueous']:
		matplotlib.rcParams.update(theme)
		palette = ['darkblue', 'lightblue', 'silver', 'gray']
		ax = sns.barplot(x="AA", y=case, hue='type', data=df_new,
							palette=palette, order=aa_by_group)
		ax.set_xlabel('AA')
		ax.set_ylabel('Probability distribution')
		# plt.suptitle('RAbD Benchmark {} targets'.format(N))
		plt.tight_layout(rect=[0, 0.03, 1, 0.9])
		plt.yticks(np.arange(0, 0.5, 0.1))
		plt.legend(bbox_to_anchor=(0.5, 1.01), loc='lower center', borderaxespad=0, prop={'size': 7})
		plt.savefig('{}_distribution_franklin.png'.format(case), dpi=600, transparent=True)
		plt.close()

def calculate_aminoacid_distribution_from_design_seq( result_csv, nativelist, suffix ): 
	"""
		A function for calculating amino acid design distributions

		Arguments: 
			nativelist = List of native PDBs
			results_csv = This is the csv file with pdb name, sequence recovery and maximum recovery sequence
							generated from proteinmpnn
		Assumption: 
			the design sequence folds into the native backbone. We have not tested this assumption yet. 
  		Output:
			a frequency count relative to position in the membrane.   
	"""

	## Generating an amino acid dictionary for 3-letter to 1-letter conversions
	amino_acids = {'A':'ala','C':'cys','D':'asp','E':'glu',
				'F':'phe','G':'gly','H':'his','I':'ile', 
				'K':'lys','L':'leu','M':'met','N':'asn',
				'P':'pro','Q':'gln','R':'arg',
				'S':'ser','T':'thr','V':'val','W':'trp','Y':'tyr'}

	print("Generating distributions for source " + "native and design" )
	import seaborn as sns
	import matplotlib
	import matplotlib.pyplot as plt

	# Initialize an AA info object for each amino acid
	ala = AminoAcidInfo()
	cys = AminoAcidInfo()
	asp = AminoAcidInfo()
	glu = AminoAcidInfo()
	phe = AminoAcidInfo()
	gly = AminoAcidInfo()
	his = AminoAcidInfo()
	ile = AminoAcidInfo()
	lys = AminoAcidInfo()
	leu = AminoAcidInfo()
	met = AminoAcidInfo()
	asn = AminoAcidInfo()
	pro = AminoAcidInfo()
	gln = AminoAcidInfo()
	arg = AminoAcidInfo()
	ser = AminoAcidInfo()
	thr = AminoAcidInfo()
	val = AminoAcidInfo()
	trp = AminoAcidInfo()
	tyr = AminoAcidInfo()

	results_mpnn = pd.read_csv(result_csv, sep=' ')

	with open( nativelist, 'rt' ) as f: 
			contents_native = f.readlines()
			contents_native = [ x.strip() for x in contents_native ]
	f.close()

	# pdbs = []
	# for pdbfile in contents: 
	# 	if ( os.path.isfile( pdbfile ) ): 
	# 		pdb = pose_from_pdb( pdbfile )
	# 		add_memb = pyrosetta.rosetta.protocols.membrane.AddMembraneMover( "from_structure" )
	# 		add_memb.apply( pdb )

	# Make a progress bar
	bar = Bar('Reading amino acids in PDB files \n', max=len(contents_native))

	aa_name = np.array(['ala','cys','asp','glu','phe','gly','his','ile','lys','leu','met','asn','pro','gln','arg','ser','thr','val','trp','tyr'])
	aa_name_letter = np.array(['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'])
	aa_by_group = ['A','I','L','M','V','N','Q','S','T','H','K','R','D','E','F','Y','W','C','G','P']
	aa_list = np.arange(1,21)

	# Loop through each PDB and encode case statement
	# for pose_num in range(len(contents_native)): 
	
	native = []

	pdb_name = ((contents_native[pose_num].split('/')[-1]).split('.')[0]).split('_')[0]
	print(pdb_name)

	native = pose_from_pdb(contents_native[pose_num])
	add_memb = pyrosetta.rosetta.protocols.membrane.AddMembraneMover( "from_structure" )
	# add_memb.apply( pdb )
	# print('mem atom added to pdb')
	# add_memb = pyrosetta.rosetta.protocols.membrane.AddMembraneMover( "from_structure" )
	add_memb.apply(native)
	file_output_aqueous = open(pdb_name+ "aqueous_res.list", "w")
	file_output_interface = open(pdb_name+ "interface_res.list","w")
	file_output_lipid = open(pdb_name+ "lipid_res.list","w")

	design_seq = results_mpnn[results_mpnn['pdb']==pdb_name]['max_seq'].item()

	if(native.total_residue()-1 == len(design_seq)):
		#loops over residues except the mem atom
		for i in range(1,native.total_residue()): 
			# print(i)
			# Get the hydration value for each side chain
			hyd = native.conformation().membrane_info().implicit_lipids().f_hydration( native.residue(i).xyz("CA") )
			thk = native.conformation().membrane_info().implicit_lipids().f_depth( native.residue(i).xyz("CA")[2] )
			if(native.conformation().membrane_info().implicit_lipids().has_pore()):
				cavity = native.conformation().membrane_info().implicit_lipids().f_cavity( native.residue(i).xyz("CA") )
			else:
				cavity = 0.0
			
			#the design and native poses align
			#find the index
			index = np.where(aa_name == native.residue(i).name3().lower())
			index = index[0][0]
			index_design = np.where(aa_name == amino_acids[design_seq[i-1]])
			index_design = index_design[0][0]

			if(thk<0.25 and cavity<0.50):
				eval(aa_name[index]).lipid_count = eval(aa_name[index]).lipid_count+1
				eval(aa_name[index_design]).designlipid_count = eval(aa_name[index_design]).designlipid_count + 1
				file_output_lipid.write('{}'.format(i))
				file_output_lipid.write("\n")
			elif(thk>0.25 and thk<0.75 and cavity<0.50):
				eval(aa_name[index]).interface_count = eval(aa_name[index]).interface_count+1
				eval(aa_name[index_design]).designinterface_count = eval(aa_name[index_design]).designinterface_count+1
				file_output_interface.write('{}'.format(i))
				file_output_interface.write("\n")
			elif(thk>0.75 or cavity>0.50):		
				eval(aa_name[index]).aqueous_count = eval(aa_name[index]).aqueous_count+1
				eval(aa_name[index_design]).designaqueous_count = eval(aa_name[index_design]).designaqueous_count+1
				file_output_aqueous.write('{}'.format(i))
				file_output_aqueous.write("\n")

	else:
		# print(contents_pdb[pose_num])
		print(contents_native[pose_num])
		print("check the native and design path list!")
		sys.exit()
	bar.next()
	file_output_lipid.close()
	file_output_interface.close()
	file_output_aqueous.close()

	bar.finish()

	
	final=[]
 
	for aa in aa_list:
		
		# print(eval(aa_name[aa-1]).lipidfacing_sequence/len(pdblist))
		# break
		
		lipid_facing = (eval(aa_name[aa-1]).lipid_count)
		interface = (eval(aa_name[aa-1]).interface_count)
		aqueous = (eval(aa_name[aa-1]).aqueous_count)
		design_lipid_facing = (eval(aa_name[aa-1]).designlipid_count)
		design_interface = (eval(aa_name[aa-1]).designinterface_count)
		design_aqueous = (eval(aa_name[aa-1]).designaqueous_count)
		final.append(np.stack([aa_name_letter[aa-1], lipid_facing, design_lipid_facing, interface, design_interface, aqueous, design_aqueous]))

 	# print(final)
 		
	df_final = pd.DataFrame(final,columns=['AA','native_lipid','design_lipid','native_interface','design_interface','aqueous', 'design_aqueous'])
	df_final.to_csv('output_protein_design_distribution.dat', sep = " ", index=False)
	##============================================
	df_final = pd.read_csv('output_protein_design_distribution.dat', sep = " ")
	print(df_final)
	print(df_final['native_lipid'].sum())
	df_final['p_native_lipid'] = round(df_final['native_lipid']/(df_final['native_lipid'].sum()),3)
	df_final['p_design_lipid'] = round(df_final['design_lipid']/(df_final['design_lipid'].sum()),3)
	df_final['p_native_interface'] = round(df_final['native_interface']/(df_final['native_interface'].sum()),3)
	df_final['p_design_interface'] = round(df_final['design_interface']/(df_final['design_interface'].sum()),3)
	df_final['p_native_aqueous'] = round(df_final['aqueous']/(df_final['aqueous'].sum()),3)
	df_final['p_design_aqueous'] = round(df_final['design_aqueous']/(df_final['design_aqueous'].sum()),3)
	df_final.to_csv('output_protein_design_distribution.dat', sep = " ", index=False)
	df_new_native = df_final[['AA','p_native_lipid','p_native_interface','p_native_aqueous']]
	df_new_native.columns=['AA','lipid','interface','aqueous']
	df_new_native['type'] = 'native'

	df_new_design = df_final[['AA','p_design_lipid','p_design_interface','p_design_aqueous']]
	df_new_design.columns=['AA','lipid','interface','aqueous']
	df_new_design['type'] = 'design'


	df_new=pd.concat([df_new_native,df_new_design])
	print(df_new)
	# print('----reorder------')

	theme = {'axes.grid': True,
		'grid.linestyle': '',
		'xtick.labelsize': 18,
		'ytick.labelsize': 18,
		"font.weight": 'regular',
		'xtick.color': 'black',
		'ytick.color': 'black',
		"axes.titlesize": 20,
		"axes.labelsize": 18
	}

	for case in ['lipid','interface','aqueous']:
		matplotlib.rcParams.update(theme)
		palette = ['darkblue', 'lightblue', 'silver', 'gray']
		ax = sns.barplot(x="AA", y=case, hue='type', data=df_new,
							palette=palette, order=aa_by_group)
		ax.set_xlabel('AA')
		ax.set_ylabel('Probability distribution')
		# plt.suptitle('RAbD Benchmark {} targets'.format(N))
		plt.tight_layout(rect=[0, 0.03, 1, 0.9])
		plt.yticks(np.arange(0, 0.5, 0.1))
		plt.legend(bbox_to_anchor=(0.5, 1.01), loc='lower center', borderaxespad=0, prop={'size': 7})
		plt.savefig('{}_distribution_mpnn.png'.format(case), dpi=600, transparent=True)
		plt.close()
 
	
def calculate_overallaminoacid_distribution( pdblist, nativelist, suffix ): 
	"""
		A function for calculating overall amino acid design distributions without any filter

		Arguments: 
			nativelist = List of native PDBs
			designedlist = List of designed PDBs
		Output:
			a frequency count relative to position in the membrane.   
	"""

	print("Generating distributions for source " + "native and design" )

	# Initialize an AA info object for each amino acid
	ala = AminoAcidInfo()
	cys = AminoAcidInfo()
	asp = AminoAcidInfo()
	glu = AminoAcidInfo()
	phe = AminoAcidInfo()
	gly = AminoAcidInfo()
	his = AminoAcidInfo()
	ile = AminoAcidInfo()
	lys = AminoAcidInfo()
	leu = AminoAcidInfo()
	met = AminoAcidInfo()
	asn = AminoAcidInfo()
	pro = AminoAcidInfo()
	gln = AminoAcidInfo()
	arg = AminoAcidInfo()
	ser = AminoAcidInfo()
	thr = AminoAcidInfo()
	val = AminoAcidInfo()
	trp = AminoAcidInfo()
	tyr = AminoAcidInfo()

	with open( pdblist, 'rt' ) as f: 
			contents_pdb = f.readlines()
			contents_pdb = [ x.strip() for x in contents_pdb ]
	f.close()

	with open( nativelist, 'rt' ) as f: 
			contents_native = f.readlines()
			contents_native = [ x.strip() for x in contents_native ]
	f.close()

	# pdbs = []
	# for pdbfile in contents: 
	# 	if ( os.path.isfile( pdbfile ) ): 
	# 		pdb = pose_from_pdb( pdbfile )
	# 		add_memb = pyrosetta.rosetta.protocols.membrane.AddMembraneMover( "from_structure" )
	# 		add_memb.apply( pdb )

	# Make a progress bar
	bar = Bar('Reading amino acids in PDB files \n', max=len(contents_pdb))

	aa_name = np.array(['ala','cys','asp','glu','phe','gly','his','ile','lys','leu','met','asn','pro','gln','arg','ser','thr','val','trp','tyr'])
	aa_list = np.arange(1,21)

	# Loop through each PDB and encode case statement
	for pose_num in range(len(contents_pdb)): 
		pdb =[]
		native = []

		pdb = pose_from_pdb(contents_pdb[pose_num])
		native = pose_from_pdb(contents_native[pose_num])
		add_memb = pyrosetta.rosetta.protocols.membrane.AddMembraneMover( "from_structure" )
		add_memb.apply(native)
	
		if(native.total_residue() == pdb.total_residue()):
			#loops over residues except the mem atom
			for i in range(1,native.total_residue()): 
				
				#the design and native poses align
				#find the index
				index = np.where(aa_name == native.residue(i).name3().lower())
				index = index[0][0]
				index_design = np.where(aa_name==pdb.residue(i).name3().lower())
				index_design = index_design[0][0]
    
				eval(aa_name[index]).nativeoverall_count = eval(aa_name[index]).nativeoverall_count+1
				eval(aa_name[index_design]).designoverall_count = eval(aa_name[index_design]).designoverall_count + 1

				
				

		else:
			print(contents_pdb[pose_num])
			print(contents_native[pose_num])
			print("check the native and design path list!")
			sys.exit()
		bar.next()
	bar.finish()

	
	final=[]
 
	for aa in aa_list:
		
		# print(eval(aa_name[aa-1]).lipidfacing_sequence/len(pdblist))
		# break
		
		native_count = (eval(aa_name[aa-1]).nativeoverall_count)
		design_count = (eval(aa_name[aa-1]).designoverall_count)
		
		final.append(np.stack([aa, native_count, design_count]))
	# print(final)
 		
	df_final = pd.DataFrame(final,columns=['AA','native_lipid','design_lipid','native_interface','design_interface','aqueous', 'design_aqueous'])
	df_final.to_csv('output_protein_design_distribution.dat', sep = " ", index=False)
       
	# print(df_final)
	# print('----reorder------')
 	
	
 
def get_heatmap( df, output_filename ):
    # ["ala","ile","leu","met","val","asn","gln","ser","thr","his","lys","arg","asp","glu", "phe","trp","tyr","cys","gly","pro"]
	aa_label = ["A","I","L","M","V","N","Q","S","T","H","K","R","D","E","F","W","Y","C","G","P"]
	res_labels = aa_label[::-1]
	fig, ax = plt.subplots(figsize=(10,10))
	# cmap = sns.color_palette("RdBu_r", 100)
	cmap = sns.color_palette("rocket_r", 100)
	# cmap = sns.color_palette("Blues", 100)
    
	ax = sns.heatmap(df,vmin=0,vmax=1.0,robust = True,linewidth = 0.5,cmap=cmap,center = 0.50, cbar=True,cbar_kws ={"label":'% Sequence recovery', "shrink":0.80},linecolor='black',annot = True,annot_kws={"size":5, "fontweight":250})
	# ax = sns.heatmap(df,vmin=0,vmax=1.0,robust = True,linewidth = 0.5,cmap=cmap,center = 0.50, cbar=True,cbar_kws ={"label":'Perplexity', "shrink":0.80},linecolor='black',annot = True,annot_kws={"size":5, "fontweight":250})
	
 	# cbar.set_label('$%$ Sequence recovery', fontsize=50)
	ax.set_aspect("equal")
	ax.set_xticklabels(aa_label, fontsize=20, fontname="Arial")
	ax.set_yticklabels(res_labels,fontsize=20, fontname="Arial")
	ax.axhline(y=0, color='k',linewidth=2)
	# ax.axhline(y=data_corr.shape[1], color='k',linewidth=1)
	ax.axvline(x=0, color='k',linewidth=2)
	# ax.axvline(x=data_corr.shape[0], color='k',linewidth=1)
	# plt.colorbar()
	axis_font = {'fontname':'Arial', 'size':'20'}
	plt.xticks(rotation=45)
	plt.yticks(rotation=0)
	plt.xlabel('Native residues', **axis_font)
	plt.ylabel('Design residues', **axis_font)
	ax.figure.axes 
	ax.tick_params(axis='both', which='major', pad=15, width=2)
	plt.savefig( output_filename+ ".png", bbox_inches="tight", transparent=True, dpi=300, format="png" )
 
def calculate_design_distribution( pdblist ):
	with open( pdblist, 'rt' ) as f: 
		contents_pdb = f.readlines()
		contents_pdb = [ x.strip() for x in contents_pdb ]
	f.close()

	# Make a progress bar
	bar = Bar('Reading amino acids in PDB files \n', max=len(contents_pdb))

	polar = 0
	nonpolar = 0
	charged = 0
	# ncharged = 0
	special = 0
	aromatic = 0
	total = 0

	# Loop through each PDB and encode case statement
	for pose_num in range(len(contents_pdb)-1): 
		pdb =[]


		pdb = pose_from_pdb(contents_pdb[pose_num])

		for i in range(1,pdb.total_residue()+1): 
			# print(i)

			#non-polar: A,I,L,M,V
			#Polar: N,Q,S,T
			#+ charge: H,K,R
			#-charge: D,E
			#aromatic:F,Y,W
			#special: C,G,P
			#the design and naive poses align
			if ( pdb.residue(i).name1() in ["A","I","L","M","V"] ):
				nonpolar = nonpolar+1
				total=total+1
			elif(pdb.residue(i).name1() in ["N","Q","S","T"]):
				polar = polar+1
				total=total+1
			elif(pdb.residue(i).name1() in ["H","K","R","D","E"]):
				charged=charged+1
				total=total+1
			elif(pdb.residue(i).name1() in ["F","W","Y"]):
				aromatic=aromatic+1
				total=total+1
			else:
				special=special+1
				total=total+1
		
		bar.next()
	bar.finish()
 
	category = ["non-polar","polar","charged","aromatic","special","total"]
	design = [nonpolar,polar,charged,aromatic,special,total]
	print(design)
	
	design_distribution = pd.DataFrame([category, design], columns=["category","design"])
	design_distribution.to_csv('output_design_distribution.txt', sep=" ", index=False)
 
def generate_heatmap_from_file( path ):
    
	residue_type = ["ala","ile","leu","met","val","asn","gln","ser","thr","his","lys","arg","asp","glu", "phe","trp","tyr","cys","gly","pro"]
	residue_type_rev = residue_type[::-1]
	
	lipid_df = pd.read_csv(path+'/output_lipid_sr_heatmap_overall.txt', delimiter = " ", index_col = 0, header=0)
	interface_df = pd.read_csv(path+'/output_interface_sr_heatmap_overall.txt', delimiter = " ", index_col=0, header = 0)
	aqueous_df = pd.read_csv(path+'/output_aqueous_sr_heatmap_overall.txt', delimiter =" ", index_col=0, header = 0)
	print(lipid_df)
	print('------------------Transpose-------------------')
	lipid_df_mod = lipid_df.T
	interface_df_mod = interface_df.T
	aqueous_df_mod = aqueous_df.T

	lipid_df_mod = lipid_df_mod.loc[residue_type_rev]
	interface_df_mod = interface_df_mod.loc[residue_type_rev]
	aqueous_df_mod = aqueous_df_mod.loc[residue_type_rev]

	print(lipid_df_mod)
	print(interface_df_mod)
	print(aqueous_df_mod)
	
	# for i in lipid_df_mod.columns:
	# 	for j in lipid_df_mod.rows:
	# 		print(i,j)
	# 		# lipid_df_mod[i][j] = -lipid_df_mod[i][j]*math.log2(lipid_df_mod[i][j])
	
	# lipid_df_mod.to_csv(path+'/output_lipid_sr_heatmap_perplex.txt',sep=" ", index=False)

	get_heatmap(lipid_df_mod,'lipid_heatmap')
	get_heatmap(interface_df_mod,'interface_heatmap')
	get_heatmap(aqueous_df_mod, 'aqueous_heatmap')


def pore_hydration_distribution( pdblist, nativelist ):
    
	print("Generating distributions for source " + "native and design" )

	with open( pdblist, 'rt' ) as f: 
			contents_pdb = f.readlines()
			contents_pdb = [ x.strip() for x in contents_pdb ]
	f.close()
 
	with open( nativelist, 'rt' ) as f: 
			contents_native = f.readlines()
			contents_native = [ x.strip() for x in contents_native ]
	f.close()

 	
	# Make a progress bar
	bar = Bar('Reading amino acids in PDB files \n', max=len(contents_pdb))
	
	np_f_hyd = []
	np_f_thk = []
	np_f_cav = []
	# Loop through each PDB and encode case statement
	for pose_num in range(len(contents_pdb)): 
		pdb =[]
		native = []
    
		pdb = pose_from_pdb(contents_pdb[pose_num])
		native = pose_from_pdb(contents_native[pose_num])
		add_memb = pyrosetta.rosetta.protocols.membrane.AddMembraneMover( "from_structure" )
		# add_memb.apply( pdb )
		# print('mem atom added to pdb')
		# add_memb = pyrosetta.rosetta.protocols.membrane.AddMembraneMover( "from_structure" )
		add_memb.apply(native)
		# print('mem atom added to native')
		# print('native length:'+str(native.total_residue()))
		# print('pdb length:'+str(pdb.total_residue()))
		
		if(native.total_residue() == pdb.total_residue()):
			
			for i in range(1,native.total_residue()+1): 
				# print(i)
				# Get the hydration value for each side chain
				hyd = native.conformation().membrane_info().implicit_lipids().f_hydration( pdb.residue(i).xyz("CA") )
				thk = native.conformation().membrane_info().implicit_lipids().f_depth( pdb.residue(i).xyz("CA")[2] )
				# if(thk != 1.0):
				# 	cavity = (hyd - thk)/(1-thk)
				# else:
				# 	cavity = 0.0
       			# cavity = native.conformation().membrane_info().implicit_lipids().f_cavity( pdb.residue(i).xyz("CA") )
				# print(not native.conformation().membrane_info().implicit_lipids().has_pore())
				if(native.conformation().membrane_info().implicit_lipids().has_pore()):
					cavity = native.conformation().membrane_info().implicit_lipids().f_cavity( pdb.residue(i).xyz("CA") )
				else:
					cavity = 0.0
    			
       
       			# a = pdb.residue(i).xyz("CA")
				# print("a: " + str(a))
				# print("z is:" + str(a[2]))
				# print("z is:" + str(a.z()))
				# sys.exit()
				# thk = native.conformation().membrane_info().implicit_lipids().f_depth( pdb.residue(i).xyz("CA")[2] )
				# cavity = native.conformation().membrane_info().implicit_lipids().f_cavity( pdb.residue(i).xyz("CA") )
				# if(pose_num==0 and i==1):
				# 	np_f_thk = thk
				# 	np_f_hyd = hyd
				# 	np_f_cav = cavity
				# else:
				# 	np_f_thk = np.vstack((np_f_thk, thk))
				# 	np_f_hyd = np.vstack((np_f_hyd, hyd))
				# 	np_f_cav = np.vstack((np_f_cav, cavity))
				np_f_thk.append(thk)
				np_f_hyd.append(hyd)
				np_f_cav.append(cavity)
			
	# print(np_f_cav)
	# print(np_f_thk)
	# print(np_f_hyd)
	
 
	final_ar = np.array([np_f_hyd, np_f_thk, np_f_cav])
	# print(final_ar)
	property_distribution = pd.DataFrame(final_ar)
	property_distribution_mod = property_distribution.T
	property_distribution_mod.columns=["hydration","thickness","cavity"]
	print(property_distribution_mod)
	# sys.exit()
	property_distribution_mod.to_csv('membrane_property_distribution.txt', sep=" ", index=False)
	print("writing distribution file")
	#extract the df which have range of f_hyd
	property_distribution_mod = pd.read_csv('membrane_property_distribution.txt', delimiter = " ", header=None, skiprows=1)
	property_distribution_mod.columns=["hydration","thickness","cavity"]
	print(property_distribution_mod)
 
	property_aqeous = property_distribution_mod[property_distribution_mod["thickness"]<0.5]
	bins = np.linspace(0, 1.2, 20)

	fig, ax = plt.subplots(figsize=(10,10))
	plt.hist(property_aqeous["thickness"], bins, alpha=0.5, label='thickness')
	plt.hist(property_aqeous["cavity"], bins, alpha=0.5, label='cavity')
	plt.hist(property_aqeous["hydration"], bins, alpha=0.5, label='hydration')
	plt.legend(loc='upper right')
	print("plotting file")
	output_filename = 'membrane_property_distribution_aqua_thk05'
	plt.savefig( output_filename+ ".png", bbox_inches="tight", transparent=True, dpi=300, format="png" )

def pore_hydration_distribution_file():
    
	property_distribution_mod = pd.read_csv('membrane_property_distribution.txt', delimiter = " ", header=None, skiprows=1)
	property_distribution_mod.columns=["hydration","thickness","cavity"]
	print(property_distribution_mod)
	# sys.exit()
	property_aqeous_temp = property_distribution_mod[property_distribution_mod["thickness"]>0.7]
	property_aqeous = property_aqeous_temp
	#  [property_aqeous_temp["thickness"]<0.7]
	
	bins = np.linspace(0, 1.2, 20)

	fig, ax = plt.subplots(3, figsize=(30,10))
	ax[0].hist(property_aqeous["thickness"], bins, alpha=0.5, label='thickness')
	ax[0].set_title('thickness')
	ax[1].hist(property_aqeous["cavity"], bins, alpha=0.5, label='cavity')
	ax[1].set_title('cavity')
	ax[2].hist(property_aqeous["hydration"], bins, alpha=0.5, label='hydration')
	ax[2].set_title('hydration')
 	# plt.hist(property_aqeous["thickness"], bins, alpha=0.5, label='thickness')
	# plt.hist(property_aqeous["cavity"], bins, alpha=0.5, label='cavity')
	# plt.hist(property_aqeous["hydration"], bins, alpha=0.5, label='hydration')
	# plt.legend(loc='upper right')
	print("plotting file")
	output_filename = 'membrane_property_distribution_thkgt07'
	plt.savefig( output_filename+ ".png", bbox_inches="tight", transparent=True, dpi=300, format="png" )
 
def calculate_residue(nativelist):
	with open( nativelist, 'rt' ) as f: 
		contents_native = f.readlines()
		contents_native = [ x.strip() for x in contents_native ]
	f.close()

 	
	# Make a progress bar
	bar = Bar('Reading amino acids in PDB files \n', max=len(contents_native))
	
	filename = "pdblist_cystein_distribution.txt"
	f = open(filename, 'w')
 	# Loop through each PDB and encode case statement
	for pose_num in range(len(contents_native)): 
		
		native = []
    
		native = pose_from_pdb(contents_native[pose_num])
		
		count = 0
		for i in range(1,native.total_residue()): 
			if(native.residue(i).name1() in ["C"]):
				count = count+1
		# print("pdb:{} #cys:{}".format(contents_native[pose_num], count))
  
		if(count>0):
			f.write( f'{contents_native[pose_num]} {count}\n' )
    
	f.close()


def calculate_hydropathy(pdb_list, suffix=['']):
	"""
		Based on Marissas code, Anastassias lab. 
	"""
 	#@title function to compute the GRAVY hydropathy score of the surface positions
	hydro_scale = {"A":1.800, "C":2.500, "D":-3.500, "E":-3.500, "F":2.800, "G":-0.400, "H":-3.200, "I":4.500, "K":-3.900, "L":3.800, "M":1.900, "N":-3.500, "P":-1.600, "Q":-3.500, "R":-4.500, "S":-0.800, "T":-0.700, "V":4.200, "W":-0.900, "Y":-1.300}

	with open( pdb_list, 'rt' ) as f: 
				contents_pdb = f.readlines()
				contents_pdb = [ x.strip() for x in contents_pdb ]
	f.close()
	pdb_names = []
	avg_interface_hp = []
	avg_lipid_hp = []
	avg_aqueous_hp = []
		
	for pose_num in range(len(contents_pdb)):
		pose = pose_from_pdb(contents_pdb[pose_num])
		pose_seq = list(pose.sequence())
		#print(pose_seq)
		#sys.exit()
		pdb_name = ((contents_pdb[pose_num].split('/')[-1]).split('.')[0]).split('_')[0]
		print(pdb_name)
		with open("/scratch16/jgray21/rsamant2/data/franklin2021/sequence-recovery/protein_mpnn/no_design_C/" + pdb_name+'aqueous_res.list') as f:
			aqueous_res = f.readlines()
			aqueous_res = [x.strip() for x in aqueous_res]
			aqueous_res = [int(x) for x in aqueous_res]
		f.close()
                
		with open("/scratch16/jgray21/rsamant2/data/franklin2021/sequence-recovery/protein_mpnn/no_design_C/" + pdb_name+'interface_res.list') as f:
			interface_res = f.readlines()
			interface_res = [x.strip() for x in interface_res]
			interface_res = [int(x) for x in interface_res]
		f.close()
			
		with open("/scratch16/jgray21/rsamant2/data/franklin2021/sequence-recovery/protein_mpnn/no_design_C/" + pdb_name+'lipid_res.list') as f:
			lipid_res = f.readlines()
			lipid_res = [x.strip() for x in lipid_res]
			lipid_res = [int(x) for x in lipid_res]
		f.close()
		aqueous_hp = 0.0
		lipid_hp = 0.0
		interface_hp = 0.0
  
		for position in range(len(pose_seq)):
			if(position in aqueous_res):
				aqueous_hp += hydro_scale[pose_seq[position - 1]] 			
			elif(position in lipid_res):
				lipid_hp += hydro_scale[pose_seq[position - 1]]
			elif(position in interface_res):
				interface_hp += hydro_scale[pose_seq[position - 1]]
		# print(aqueous_hp)
		# print(lipid_hp)
		# print(interface_hp)
		aqueous_hp = aqueous_hp/len(aqueous_res)
		lipid_hp = lipid_hp/len(lipid_res)
		interface_hp = interface_hp/len(interface_res)
		# print(aqueous_hp)
		# print(lipid_hp)
		# print(interface_hp)
		# sys.exit()
		pdb_names.append(pdb_name)
		avg_interface_hp.append(round(interface_hp,3))
		avg_lipid_hp.append(round(lipid_hp,3))
		avg_aqueous_hp.append(round(aqueous_hp,3))
  
	final = np.stack(( pdb_names, avg_lipid_hp, avg_interface_hp, avg_aqueous_hp ), axis=1)
	df = pd.DataFrame(final,columns=['pdb','avg_lipid_hydropathy', 'avg_interface_hydropathy', 'avg_aqueous_hydropathy'])

	# print(len(df[df['pdb']=='2KSD']['max_seq'].item()))
	# design_seq = df[df['pdb']=='2KSD']['max_seq'].item()
	# print(design_seq[0])
	output_file = './franklin_hydropathy'+suffix+'.dat'
	df.to_csv(output_file,sep=' ', index=False)  
    
