.SILENT:

VMFILES = pkgIndex.tcl modelmaker.tcl CColor.tcl README_requirements.txt analysis_mdff.tcl auto_mdff.tcl check_clusterjob.tcl clash_selection.tcl find_selection.tcl \
helix-turn-cc.tcl make_folders_rosetta.sh makepsf.tcl package_example.tcl package_read_example.tcl pdb_density_selector.tcl rosetta_input_generator.tcl rosetta_scoring_new.tcl \
rosetta_scoring_wrapper.tcl rosetta_vmd.tcl run_situs.sh ss_analysis.tcl sstemp.tcl utilities.tcl

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

	
