package ifneeded modelmaker 0.1 "set env(RosettaVMDDIR) [list $dir]; [list source [file join $dir modelmaker.tcl]]"
package ifneeded RosettaVMD 0.1 "set env(RosettaVMDDIR) [list $dir]; [list source [file join $dir rosetta_vmd.tcl]]"
package ifneeded RosettaInputGenerator 0.1 [list source [file join $dir rosetta_input_generator.tcl]]
package ifneeded MakePsf 0.1 [list source [file join $dir makepsf.tcl]]
package ifneeded RosettaScoring 0.1 [list source [file join $dir rosetta_scoring_wrapper.tcl]]
package ifneeded CCColor 0.1 [list source [file join $dir CCColor.tcl]]
package ifneeded FindSelection 0.1 [list source [file join $dir find_selection.tcl]]
package ifneeded RosettaUtilities 0.1 [list source [file join $dir utilities.tcl]]
package ifneeded SSAnalysis 0.1 [list source [file join $dir ss_analysis.tcl]]
