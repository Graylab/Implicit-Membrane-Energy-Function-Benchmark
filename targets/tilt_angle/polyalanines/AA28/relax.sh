~/Softwares/Rosetta/main/source/bin/relax.linuxgccrelease -database ~/Softwares/Rosetta/main/database/ -relax:constrain_relax_to_start_coords -in:file:s AA28_0001_0001.pdb -ex1 -ex2 -ex1aro -ex2aro -overwrite -nstruct 6 -out:path:pdb ./ -out:path:score ./ -relax:ramp_constraints false
