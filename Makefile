.SILENT:

VMFILES = pkgIndex.tcl modelmaker.tcl CColor.tcl README.md analysis_mdff.tcl auto_mdff.tcl clash_selection.tcl find_selection.tcl \
helix-turn-cc.tcl makepsf.tcl pdb_density_selector.tcl rosetta_input_generator.tcl \
rosetta_scoring_wrapper.tcl rosetta_vmd.tcl ss_analysis.tcl utilities.tcl

VMVERSION = 0.1
DIR = $(PLUGINDIR)/noarch/tcl/mdodelmaker$(VMVERSION)

bins:
win32bins:
dynlibs:
staticlibs:
win32staticlibs:

distrib:
	@echo "Copying modelmaker $(VMVERSION) files to $(DIR)"
	mkdir -p $(DIR) 
	cp $(VMFILES) $(DIR) 

	
