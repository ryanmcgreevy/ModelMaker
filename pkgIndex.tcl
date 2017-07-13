package ifneeded modelmaker 0.1 "set env(RosettaVMDDIR) [list $dir]; [list source [file join $dir modelmaker.tcl]]"

package ifneeded RosettaVMD 0.1 "set env(RosettaVMDDIR) [list $dir]; [list source [file join $dir RosettaVMD/rosetta_vmd.tcl]]"
package ifneeded RosettaInputGenerator 0.1 [list source [file join $dir RosettaVMD/rosetta_input_generator.tcl]]
package ifneeded MakePsf 0.1 [list source [file join $dir RosettaVMD/makepsf.tcl]]
package ifneeded AutoMDFF 0.1 [list source [file join $dir RosettaVMD/auto_mdff.tcl]]
package ifneeded RosettaScoring 0.1 [list source [file join $dir RosettaVMD/rosetta_scoring_wrapper.tcl]]
package ifneeded AnalysisMDFF 0.1 [list source [file join $dir RosettaVMD/analysis_mdff.tcl]]
package ifneeded CCColor 0.1 [list source [file join $dir RosettaVMD/CCColor.tcl]]
package ifneeded FindSelection 0.1 [list source [file join $dir RosettaVMD/find_selection.tcl]]
package ifneeded DensitySelector 0.1 [list source [file join $dir RosettaVMD/pdb_density_selector.tcl]]
package ifneeded CheckCluster 0.1 [list source [file join $dir RosettaVMD/check_clusterjob.tcl]]
package ifneeded RosettaUtilities 0.1 [list source [file join $dir RosettaVMD/utilities.tcl]]
package ifneeded SSAnalysis 0.1 [list source [file join $dir RosettaVMD/ss_analysis.tcl]]
