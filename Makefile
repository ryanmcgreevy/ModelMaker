.SILENT:

VMFILES = CCColor.tcl  clash_selection.tcl  find_selection.tcl  helix-turn-cc.tcl  \
makepsf.tcl  modelmaker.tcl  pkgIndex.tcl  rosetta_input_generator.tcl  \
rosetta_scoring_wrapper.tcl  rosetta_vmd.tcl  ss_analysis.tcl  utilities.tcl \
align2d.py  clusterK.py  model-single.py \
README.md


VMVERSION = 0.1
DIR = $(PLUGINDIR)/noarch/tcl/modelmaker$(VMVERSION)

bins:
win32bins:
dynlibs:
staticlibs:
win32staticlibs:

distrib:
	@echo "Copying modelmaker $(VMVERSION) files to $(DIR)"
	mkdir -p $(DIR) 
	cp $(VMFILES) $(DIR) 

	
