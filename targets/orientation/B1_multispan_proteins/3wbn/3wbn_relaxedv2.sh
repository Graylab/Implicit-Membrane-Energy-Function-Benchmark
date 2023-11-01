#!/bin/bash
/home/rsamant2/Softwares/Rosetta/main/source/bin/relax.linuxgccrelease -database /home/rsamant2/Softwares/Rosetta/main/source/bin/../../database -relax:constrain_relax_to_start_coords -relax:bb_move false -relax:chi_move true -in:file:s /home/rsamant2/Softwares/Implicit-Membrane-Energy-Function-Benchmark-Electrostatics/Implicit-Membrane-Energy-Function-Benchmark-master/targets/orientation/B1_multispan_proteins/3wbn/3wbn_ignorechain.pdb -ex1 -ex2 -ex1aro -ex2aro -overwrite -nstruct 100 -out:path:pdb /home/rsamant2/Softwares/Implicit-Membrane-Energy-Function-Benchmark-Electrostatics/Implicit-Membrane-Energy-Function-Benchmark-master/targets/orientation/B1_multispan_proteins/3wbn/relaxed_structure_v2/ -out:path:score /home/rsamant2/Softwares/Implicit-Membrane-Energy-Function-Benchmark-Electrostatics/Implicit-Membrane-Energy-Function-Benchmark-master/targets/orientation/B1_multispan_proteins/3wbn/relaxed_structure_v2/
