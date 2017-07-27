####################################################
# VMD Rosetta Scripts Generator
# Maximilian Scheurer, March 2016
####################################################


# JOB CONTROLLER
#### DO NOT CHANGE ANYTTHING HERE
package require mdff

package require RosettaInputGenerator
package require AutoMDFF
package require RosettaScoring
package require CCColor
package require FindSelection
package require CheckCluster
package require RosettaUtilities
package require SSAnalysis

## RosettaVMD namespace
# main namespace for RosettaVMD package
namespace eval ::RosettaVMD {
    namespace export start_rosetta_refine ;#< refinement protocol
    namespace export start_rosetta_refine_sidechains_density ;#< sidechain refinement protocol
    namespace export start_rosetta_abinitio ;# ab-initio procotol
    namespace export analyze_abinitio
    namespace export start_mdff_run
    namespace export helix_reg
    namespace export readscorefile
	namespace export make_dx_file
	namespace export make_mrc_file
	namespace export smooth_density
	namespace export write_phenixpdb

	# Set up Variable
	set version 0.1
	set packageDescription "RosettaVMD Plugin"

    # Variable for the path of the script
    variable home [file join [pwd] [file dirname [info script]]]
}
package provide RosettaVMD $RosettaVMD::version

proc start_rosetta_refine {jobname mol selections anchor cartesian mapname mapresolution score_dens bestN nstruct {cluster 0} {nPerTask 5} {scoreOnly 0} args} \
{
  set username $::MODELMAKER::settings(username)
	# prepare configuration
	set selection_length [llength $selections]
	set find_cfg []
	for {set i 0} {$i < $selection_length} {incr i} {
		# allow backbone and sidechain movement for selections
		lappend find_cfg {1 1}
	}
	# MOL selections config {offset 4}
	# offset set to 5 for cartesianSampler
	if {$cartesian} {
		set find_sel [find_selection full_length_model/$mol $selections $find_cfg 5]
	} else {
		set find_sel [find_selection full_length_model/$mol $selections $find_cfg 0]
	}
	set spans [lindex $find_sel 0]
	set exclude [lindex $find_sel 1]
	set chains [lindex $find_sel 2]
	set constraints [lindex $find_sel 3]

	set ros_config [list $chains $spans $exclude [list $anchor $constraints]]

	######################
	# Cleanup input file
	exec sed -i -e {s/HSD/HIS/g} full_length_model/$mol.pdb
	exec sed -i -e {s/HSE/HIS/g} full_length_model/$mol.pdb
	exec sed -i -e {s/HSP/HIS/g} full_length_model/$mol.pdb
	######################
	puts "Rosetta started"
	# refine_with_rosetta {jobname MOL mapname res score_dens nstruct cluster nPerTask configuration {cartesianSampler 0}}
	refine_with_rosetta $jobname $mol.pdb $mapname $mapresolution $score_dens $nstruct $cluster $nPerTask $ros_config $cartesian
	exec chmod +x $jobname.sh
	exec mv $jobname.sh rosetta_output_$jobname/
	exec cp $mapname.mrc rosetta_input_$jobname/
	cd rosetta_output_$jobname

	if {!$scoreOnly} {
		exec mkdir -p sc_out
		exec mkdir -p pdb_out
		exec mkdir -p OUTPUT_FILES
		set output [exec "[pwd]/$jobname.sh" "$jobname" "$mol.pdb" >> rosetta_log_$jobname.log &]
		set current [exec ls -1v pdb_out | wc -l]
		while {$current < $nstruct} {
			set n 5
			puts "Files are not yet available."
			puts "Current number: [exec ls -1v pdb_out | wc -l] - [expr double($current)/($nstruct) * 100.0] %"
			after [expr {int($n * 1000)}]
			set current [exec ls -1v pdb_out | wc -l]
			if {$cluster} {
				set logfile [open "rosetta_log_$jobname.log" r]
				set dt [read $logfile]
				close $logfile
				set lns [split $dt "\n"]
				set infoline [lindex $lns 0]
				set res [regexp {([0-9]+)} $infoline jobid]
				set tasks [expr int(ceil(double($nstruct)/double($nPerTask)))]
				set status [check_clusterjob $username $jobid $tasks]
				if {$status == 0} {
					break
				}
			}
		}
		puts $output
		puts "Rosetta finished"
	}

	# Scoring
	#exec mv {*}[glob *.sc] sc_out/
	#exec mv {*}[glob *.pdb] pdb_out/

	# MOL max_structures
	if {!$cluster} {
		puts "Scoring normal run."
		score_refinement ${jobname}_$mol $bestN
	} else {
		puts "Scoring cluster run."
		score_refinement_cluster $jobname ${jobname}_$mol $bestN
	}


	cd ..
}

proc start_rosetta_refine_sidechains_density {jobname mol selections anchor mapname mapresolution score_dens bestN nstruct {cluster 0} {nPerTask 5} {scoreOnly 0} args} \
{
	# prepare configuration
	set username $::MODELMAKER::settings(username)
	set selection_length [llength $selections]
	set find_cfg []
	for {set i 0} {$i < $selection_length} {incr i} {
		# allow backbone and sidechain movement for selections
		lappend find_cfg {1 0}
	}
	# MOL selections config {offset 4}
	set find_sel [find_selection full_length_model/$mol $selections $find_cfg 0]
	set spans [lindex $find_sel 0]
	set exclude [lindex $find_sel 1]
	set chains [lindex $find_sel 2]
	set constraints [lindex $find_sel 3]

	set ros_config [list $chains $spans $exclude [list $anchor $constraints]]
	######################
	# Cleanup input file
	exec sed -i -e {s/HSD/HIS/g} full_length_model/$mol.pdb
	exec sed -i -e {s/HSE/HIS/g} full_length_model/$mol.pdb
	exec sed -i -e {s/HSP/HIS/g} full_length_model/$mol.pdb
	######################
	puts "Rosetta sidechain started"
	refine_with_rosetta $jobname $mol.pdb $mapname $mapresolution $score_dens $nstruct $cluster $nPerTask $ros_config
	# refine_sidechains_rosetta
	exec chmod +x $jobname.sh
	exec mv $jobname.sh rosetta_output_$jobname/
	exec cp $mapname.mrc rosetta_input_$jobname/
	cd rosetta_output_$jobname

	if {!$scoreOnly} {
		exec mkdir -p sc_out
		exec mkdir -p pdb_out
		exec mkdir -p OUTPUT_FILES
		set output [exec "[pwd]/$jobname.sh" "$jobname" "$mol.pdb" >> rosetta_log_$jobname.log &]
		set current [exec ls -1v pdb_out | wc -l]
		while {$current < $nstruct} {
			set n 5
			puts "Files are not yet available."
			puts "Current number: [exec ls -1v pdb_out | wc -l] - [expr double($current)/($nstruct) * 100.0] %"
			after [expr {int($n * 1000)}]
			set current [exec ls -1v pdb_out | wc -l]
			if {$cluster} {
				set logfile [open "rosetta_log_$jobname.log" r]
				set dt [read $logfile]
				close $logfile
				set lns [split $dt "\n"]
				set infoline [lindex $lns 0]
				set res [regexp {([0-9]+)} $infoline jobid]
				set tasks [expr int(ceil(double($nstruct)/double($nPerTask)))]
				set status [check_clusterjob $username $jobid $tasks]
				if {$status == 0} {
					break
				}
			}
		}
		puts $output
		puts "Rosetta finished"
	}

	# MOL max_structures
	puts [pwd]
	if {!$cluster} {
		puts "Scoring normal run."
		score_refinement ${jobname}_$mol $bestN
	} else {
		puts "Scoring cluster run."
		score_refinement_cluster $jobname ${jobname}_$mol $bestN
	}


	cd ..
}



#refinement without density (fast relax)
#sidechains_only: puts backbone constraints
proc start_rosetta_basic_refine {jobname mol selections anchor sidechains_only bestN nstruct {cluster 0} {nPerTask 5} {scoreOnly 0} args} \
{
	set username $::MODELMAKER::settings(username)
	# prepare configuration
	set selection_length [llength $selections]
	set find_cfg []
	for {set i 0} {$i < $selection_length} {incr i} {
		# allow backbone and sidechain movement for selections
		if {$sidechains_only} {
			lappend find_cfg {1 0}
		} else {
			lappend find_cfg {1 1}
		}
	}
	# MOL selections config {offset 4}
	set find_sel [find_selection full_length_model/$mol $selections $find_cfg 0]
	set spans [lindex $find_sel 0]
	set exclude [lindex $find_sel 1]
	set chains [lindex $find_sel 2]
	set constraints [lindex $find_sel 3]

	set ros_config [list $chains $spans $exclude [list $anchor $constraints]]

	######################
	# Cleanup input file
	exec sed -i -e {s/HSD/HIS/g} full_length_model/$mol.pdb
	exec sed -i -e {s/HSE/HIS/g} full_length_model/$mol.pdb
	exec sed -i -e {s/HSP/HIS/g} full_length_model/$mol.pdb
	######################
	puts "Rosetta basic refinement (without density) started"
	# rosetta_basic_refinement {jobname MOL nstruct cluster nPerTask configuration}
	rosetta_basic_refinement $jobname $mol.pdb $nstruct $cluster $nPerTask $ros_config

	exec chmod +x $jobname.sh
	exec mv $jobname.sh rosetta_output_$jobname/
	# exec cp $mapname.mrc rosetta_input_$jobname/
	cd rosetta_output_$jobname

	if {!$scoreOnly} {
		exec mkdir -p sc_out
		exec mkdir -p pdb_out
		exec mkdir -p OUTPUT_FILES
		set output [exec "[pwd]/$jobname.sh" "$jobname" "$mol.pdb" >> rosetta_log_$jobname.log &]
		set current [exec ls -1v pdb_out | wc -l]
		while {$current < $nstruct} {
			set n 5
			puts "Files are not yet available."
			puts "Current number: [exec ls -1v pdb_out | wc -l] - [expr double($current)/($nstruct) * 100.0] %"
			after [expr {int($n * 1000)}]
			set current [exec ls -1v pdb_out | wc -l]
			if {$cluster} {
				set logfile [open "rosetta_log_$jobname.log" r]
				set dt [read $logfile]
				close $logfile
				set lns [split $dt "\n"]
				set infoline [lindex $lns 0]
				set res [regexp {([0-9]+)} $infoline jobid]
				set tasks [expr int(ceil(double($nstruct)/double($nPerTask)))]
				set status [check_clusterjob $username $jobid $tasks]
				if {$status == 0} {
					break
				}
			}
		}
		puts $output
		puts "Rosetta finished"
	}

	# Scoring
	#exec mv {*}[glob *.sc] sc_out/
	#exec mv {*}[glob *.pdb] pdb_out/

	# MOL max_structures
	if {!$cluster} {
		puts "Scoring normal run."
		score_refinement ${jobname}_$mol $bestN
	} else {
		puts "Scoring cluster run."
		score_refinement_cluster $jobname ${jobname}_$mol $bestN
	}


	cd ..
}



proc start_rosetta_abinitio {jobname mol selections anchor fragfiles nstruct {cluster 0} {nPerTask 5} {testrun 0} args} \
{
	# make sure we are in the correct directory
#	cd $::MODELMAKER::workdir
	# prepare configuration
	set selection_length [llength $selections]
	set find_cfg []
	for {set i 0} {$i < $selection_length} {incr i} {
		# allow backbone and sidechain movement for selections
		lappend find_cfg {1 1}
	}

	# Cleanup input file
	exec sed -i -e {s/HSD/HIS/g} $::MODELMAKER::workdir/setup-$jobname/$mol.pdb
	exec sed -i -e {s/HSE/HIS/g} $::MODELMAKER::workdir/setup-$jobname/$mol.pdb
	exec sed -i -e {s/HSP/HIS/g} $::MODELMAKER::workdir/setup-$jobname/$mol.pdb


	# MOL selections config {offset 4}
	set find_sel [find_selection $::MODELMAKER::workdir/setup-$jobname/$mol $selections $find_cfg 0]
	set spans [lindex $find_sel 0]
	set exclude [lindex $find_sel 1]
	set chains [lindex $find_sel 2]
	set constraints [lindex $find_sel 3]

	set ros_config [list $chains $spans $exclude [list $anchor $constraints]]

	######
	#find chain identifiers for folding
	set seltexts []
	mol delete all
	foreach seltext $selections {
		lappend seltexts $seltext
	}
	set searchmol [mol new $::MODELMAKER::workdir/setup-$jobname/$mol.pdb]

	set chain_idents {}
	foreach findseltext $seltexts {
		set sel_chainfind [atomselect $searchmol "$findseltext"]
		set occuring_chains [lsort -unique [$sel_chainfind get chain]]
		foreach c $occuring_chains {
			lappend chain_idents $c
		}
	}
	mol delete all
	######################

	set username $::MODELMAKER::settings(username)
	puts "Rosetta abinitio started."
# 	rosetta_abinitio {jobname MOL fragfiles nstruct cluster nPerTask test configuration}
	#cd $::MODELMAKER::workdir/run-$jobname
	rosetta_abinitio $jobname $mol.pdb $fragfiles $nstruct $cluster $nPerTask $testrun $ros_config $chain_idents
	exec chmod +x $::MODELMAKER::workdir/run-$jobname/$jobname.sh

	#exec cp $mapname.mrc rosetta_input_$jobname/

	file mkdir $::MODELMAKER::workdir/run-$jobname/sc_out
	file mkdir $::MODELMAKER::workdir/run-$jobname/pdb_out
	file mkdir $::MODELMAKER::workdir/run-$jobname/OUTPUT_FILES
	set output [exec "$::MODELMAKER::workdir/run-$jobname/$jobname.sh" "$jobname" "$mol.pdb" >> $::MODELMAKER::workdir/run-$jobname/rosetta_log_$jobname.log &]
  set current [llength [glob -nocomplain $::MODELMAKER::workdir/run-$jobname/*.pdb ] ]
	while {$current < $nstruct} {
		set n 5
		puts "Files are not yet available."
		puts "Current number: $current - [expr double($current)/($nstruct) * 100.0] %"
		after [expr {int($n * 1000)}]
		set current [llength [glob -nocomplain $::MODELMAKER::workdir/run-$jobname/*.pdb ] ]
		if {$cluster} {
			set logfile [open "rosetta_log_$jobname.log" r]
			set dt [read $logfile]
			close $logfile
			set lns [split $dt "\n"]
			set infoline [lindex $lns 0]
			set res [regexp {([0-9]+)} $infoline jobid]
			set tasks [expr int(ceil(double($nstruct)/double($nPerTask)))]
			set status [check_clusterjob $username $jobid $tasks]
			if {$status == 0} {
				break
			}
		}
	}
#	puts "pdbs: [glob *.pdb]"
#  puts "score: [glob -nocomplain *.sc] pwd: [pwd]"
  if { [glob -nocomplain $::MODELMAKER::workdir/run-$jobname/*.sc] != ""} {
	  file rename {*}[glob -nocomplain $::MODELMAKER::workdir/run-$jobname/*.sc] $::MODELMAKER::workdir/run-$jobname/sc_out/
  }
	file rename {*}[glob -nocomplain $::MODELMAKER::workdir/run-$jobname/*.pdb] $::MODELMAKER::workdir/run-$jobname/pdb_out/

  puts $output
	puts "Rosetta abinitio finished."
}


proc start_rosetta_insertion {jobname mol selections fragfiles fasta nstruct {cluster 0} {nPerTask 5} args} \
{
	# prepare configuration
#	cd $::MODELMAKER::workdir
	set selection_length [llength $selections]
	set find_cfg []
	for {set i 0} {$i < $selection_length} {incr i} {
		# allow backbone and sidechain movement for selections
		lappend find_cfg {1 1}
	}

	# Cleanup input file
	exec sed -i -e {s/HSD/HIS/g} $::MODELMAKER::workdir/setup-$jobname/$mol.pdb
	exec sed -i -e {s/HSE/HIS/g} $::MODELMAKER::workdir/setup-$jobname/$mol.pdb
	exec sed -i -e {s/HSP/HIS/g} $::MODELMAKER::workdir/setup-$jobname/$mol.pdb


	# MOL selections config {offset 4}
	set find_sel [find_selection $::MODELMAKER::workdir/setup-$jobname/$mol $selections $find_cfg 0]
	set spans [lindex $find_sel 0]
	set exclude [lindex $find_sel 1]
	set chains [lindex $find_sel 2]

	set ros_config [list $chains $spans $exclude]

	#cd $::MODELMAKER::workdir/run-$jobname
	puts "Rosetta insertion folding started."
	rosetta_insertion $jobname $mol $fragfiles $fasta $nstruct $cluster $nPerTask $ros_config
	exec chmod +x $::MODELMAKER::workdir/run-$jobname/$jobname.sh

	file mkdir $::MODELMAKER::workdir/run-$jobname/sc_out
	file mkdir $::MODELMAKER::workdir/run-$jobname/pdb_out
	file mkdir $::MODELMAKER::workdir/run-$jobname/OUTPUT_FILES
	set output [exec "$::MODELMAKER::workdir/run-$jobname/$jobname.sh" >> $::MODELMAKER::workdir/run-$jobname/rosetta_log_$jobname.log &]
	set current [llength [glob -nocomplain $::MODELMAKER::workdir/run-$jobname/pdb_out/$jobname*.pdb ] ]
	while {$current < $nstruct} {
		set n 20
		puts "Files are not yet available."
		puts "Current number: $current - [expr double($current)/($nstruct) * 100.0] %"
		after [expr {int($n * 1000)}]
		set current [llength [glob -nocomplain $::MODELMAKER::workdir/run-$jobname/pdb_out/$jobname*.pdb ] ]
		if {$cluster} {
			set logfile [open "rosetta_log_$jobname.log" r]
			set dt [read $logfile]
			close $logfile
			set lns [split $dt "\n"]
			set infoline [lindex $lns 0]
			set res [regexp {([0-9]+)} $infoline jobid]
			set tasks [expr int(ceil(double($nstruct)/double($nPerTask)))]
			set status [check_clusterjob $username $jobid $tasks]
			if {$status == 0} {
				break
			}
		}
	}
  file mkdir $::MODELMAKER::workdir/run-$jobname/intermediates
	file rename {*}[glob $::MODELMAKER::workdir/run-$jobname/pdb_out/loops_closed*.pdb] $::MODELMAKER::workdir/run-$jobname/intermediates/
#  file rename {*}[glob $::MODELMAKER::workdir/run-$jobname/*.sc] $::MODELMAKER::workdir/run-$jobname/sc_out/
#	file rename {*}[glob $::MODELMAKER::workdir/run-$jobname/*.pdb] $::MODELMAKER::workdir/run-$jobname/pdb_out/

	puts $output
	puts "Rosetta insertion folding finished."
	#cd $::MODELMAKER::workdir
}

proc analyze_abinitio {jobname mol template bestN nstruct cluster align_template align_rosetta analysis_components {insertion 0} args} \
{
	global vmdexe
	global packagePath
	cd $::MODELMAKER::workdir/run-$jobname

	puts "Ab-initio analysis started."

	# ALIGNMENT
	puts "current directory: [pwd]"
	file mkdir pdb_out_aligned

	if {!$cluster} {
		if {$insertion != 0} {
			align_rosetta_local 1 $nstruct ${jobname}_$mol $template $insertion $align_template $align_rosetta
		} else {
			align_rosetta_local 1 $nstruct ${jobname}_$mol $template $mol $align_template $align_rosetta
		}
	} else {
		if {$insertion != 0} {
			align_rosetta_cluster 1 $nstruct ${jobname}_$mol $template $insertion $align_template $align_rosetta
		} else {
			align_rosetta_cluster 1 $nstruct ${jobname}_$mol $template $mol $align_template $align_rosetta
		}
	}

	#SCORING MOL max_structures cluster?
	# extra variable to run scoring on cluster folded insertion runs
	set extra 0
	if {$insertion != 0} {
		set extra $mol
	}

	score_abinitio $jobname ${jobname}_$mol $bestN $cluster $extra

	# ANALYSIS
	file mkdir analysis
	file copy ${jobname}_${mol}_rosetta_scoring_min_$bestN.dcd analysis
	file copy ${jobname}_${mol}_rosetta_scoring_min_$bestN.pdb analysis
	cd analysis

	set prefix ${jobname}_${mol}

	foreach ana $analysis_components {
		set type [lindex $ana 0]
		switch -exact -- $type {
			cluster {
				set start [lindex $ana 1]
				set end [lindex $ana 2]
				set cluster_number [lindex $ana 4]
				exec mkdir -p cluster_${start}_$end
				exec cp ${jobname}_${mol}_rosetta_scoring_min_$bestN.pdb cluster_${start}_$end
				cd cluster_${start}_${end}
				#CLUSTER INPUT
				make_cluster_input $prefix $start $end $bestN
				# run_vmd_clustering $mol $start $end $bestN 0.25 $cluster_number
				# TODO: refine clustering settings!
				if {$cluster} {
					run_clustering $prefix $start $end $bestN
				} else {
					run_rosetta_clustering $mol $start $end $bestN -1 $cluster_number
				}
				# run_rosetta_clustering $mol $start $end $bestN -1 $cluster_number
				# run_clustering $prefix $start $end $bestN
				cd $::MODELMAKER::workdir/run-$jobname/analysis
			}
			ss {
				set start [lindex $ana 1]
				set end [lindex $ana 2]
				set chain [lindex $ana 3]
				exec mkdir -p ss_${start}_$end
				exec cp ${jobname}_${mol}_rosetta_scoring_min_$bestN.pdb ss_${start}_$end
				cd ss_${start}_${end}

				set MOL $prefix
				set max_structures $bestN
				set resid_start $start
				set resid_stop $end
				set chain_id $chain

				set ss_file [open "run_ss.tcl" w]
				puts $ss_file "lappend ::auto_path \"$packagePath\""
				puts $ss_file "package require SSAnalysis"
				puts $ss_file "ss_analysis $MOL $max_structures $resid_start $resid_stop $chain_id"
				puts $ss_file "quit"
				close $ss_file
				puts "running SS analysis"
				exec >&@stdout $vmdexe -dispdev text -e run_ss.tcl
				evaluate_ss_analysis $bestN $start $end
				cd $::MODELMAKER::workdir/run-$jobname/analysis
			}
			default {
				cd $::MODELMAKER::workdir
				puts "Wrong input. Exiting."
			}
		}
	}

	cd $::MODELMAKER::workdir
}

# MDFF only starts with pdb created by rosetta
proc start_mdff_run {jobname mol mapname fixedselection gscale minSteps num res bestN {score 1} {cascade 0} {config 0} {gridconfig 0} {pdbfolder 0} args} \
{
	package require AnalysisMDFF
	global ch_seg
	global mutations
	global parfiles
	global topdir
	global topfiles
	global dcdfreq
	global path
	global namdArgs
	global rosettapath
	global rosettaDBpath
	global platform
	set scores {}

	### normal MDFF rosetta run
	if {$bestN > 0 && !$cascade} {
		for {set i 1} {$i <= $bestN} {incr i} {
			set folder mdff_${jobname}_$i
			set prefix ${jobname}_${mol}_best$i
			exec mkdir -p $folder

			exec cp rosetta_output_$jobname/$prefix.pdb $folder/$prefix.pdb

			exec cp $mapname.dx $folder/
			cd $folder

			#arguments: jobname MOL mapname fixedselection gscale minSteps num ch_seg topdir topfile parfile
			auto_mdff_init ${jobname}_$i $prefix $mapname $fixedselection $gscale $minSteps $num $ch_seg $mutations $topdir $topfiles $parfiles

			exec sed -i -e "s/dcdfreq.*/dcdfreq\ ${dcdfreq}/g" mdff_template.namd
			puts "Starting NAMD with job $jobname-$i"
			exec $path/namd2 $namdArgs mdff_${jobname}_$i-step1.namd > mdff_${jobname}_$i-step1.log
			puts "NAMD finished"

			puts "Analysis started."
			analyze_mdff $prefix mdff_${jobname}_$i-step1 $mapname $res $ch_seg
			puts "Analysis finished."

			set outname mdff_${jobname}_$i-step1.dcd-last

			exec sed -i -e {s/HSD/HIS/g} $outname.pdb

			#set beta & occupancy to 1.0
			mol delete all
			mol new $outname.pdb
			[atomselect top all] set beta 1.0
			[atomselect top all] set occupancy 1.0
			[atomselect top all] writepdb $outname.pdb
			mol delete all

			if {$score != 0} {
				exec $rosettapath/score.$platform -database $rosettaDBpath -out:file:scorefile ${jobname}_$i-mdff.sc -ignore_zero_occupancy false -in:file:s $outname.pdb > ${jobname}_$i-scoring.log
				lappend scores [list $i [readscorefile ${jobname}_$i-mdff.sc]]
			}
			exec cp $outname.pdb ../full_length_model/$prefix-mdff.pdb
			cd ..
		}
	} elseif {$bestN == 0} {
		global inputfolder
		if {$pdbfolder == 0} {
			set pdbfolder $inputfolder
		}
		# for normal MDFF run or cascade mdff run, input files should be in inputfolder
		puts [pwd]
		set folder mdff_${jobname}
		set prefix ${jobname}_${mol}
		exec mkdir -p $folder
		exec cp $pdbfolder/$mol.pdb $folder/$prefix.pdb
		exec cp $inputfolder/$mapname.dx $folder/
		cd $folder
		puts [pwd]
		if {$cascade} {
			for {set i 1} {$i <= [llength $config]} {incr i} {
				set smoothr [lindex $config [expr $i-1] 1]
				puts "Smooth radius $smoothr"
				exec mkdir -p $inputfolder/density_grids
				if {$smoothr == 0} {
						mdff griddx -i $inputfolder/$mapname.dx -o $inputfolder/density_grids/griddx$i.dx
						#exec /bin/cp -u $inputfolder/density_grids/griddx0.dx $inputfolder/density_grids/griddx$i.dx
					} else {
						volutil -smooth $smoothr $inputfolder/$mapname.dx -o $inputfolder/${mapname}_smooth_$smoothr.dx
						mdff griddx -i $inputfolder/${mapname}_smooth_$smoothr.dx -o $inputfolder/density_grids/griddx$i.dx
					}
			}
			auto_mdff_casc_init $jobname $prefix $mapname $fixedselection $config $minSteps $dcdfreq $ch_seg $mutations $topdir $topfiles $parfiles $inputfolder $gridconfig

			for {set i 1} {$i <= [llength $config]} {incr i} {
				puts "Running NAMD step $i"
				exec $path/namd2 $namdArgs mdff_${prefix}_step$i.namd > mdff_${prefix}_step$i.log
				puts "Finished NAMD step $i"
			}
			puts "Analysis started."
			analyze_mdff $prefix mdff_${prefix}_step[llength $config] $mapname $res $ch_seg
			puts "Analysis finished."
		} else {
			# puts "please specify a configuration!!!"
			puts "Normal MDFF run, without scoring or cascade"
			exec cp $inputfolder/$mol.pdb .

			exec cp $inputfolder/$mapname.dx .

			#arguments: jobname MOL mapname fixedselection gscale minSteps num ch_seg topdir topfile parfile
			auto_mdff_init ${jobname} $prefix $mapname $fixedselection $gscale $minSteps $num $ch_seg $mutations $topdir $topfiles $parfiles

			exec sed -i -e "s/dcdfreq.*/dcdfreq\ ${dcdfreq}/g" mdff_template.namd
			puts "Starting NAMD with job $jobname"
			exec $path/namd2 $namdArgs mdff_${jobname}-step1.namd > mdff_${jobname}-step1.log
			puts "NAMD finished"

			puts "Analysis started."
			analyze_mdff $prefix mdff_${jobname}-step1 $mapname $res $ch_seg
			puts "Analysis finished."

			set outname mdff_${jobname}-step1.dcd-last

			exec sed -i -e {s/HSD/HIS/g} $outname.pdb
		}

		cd ..
	} elseif {$bestN > 0 && $cascade} {
		global inputfolder
		if {$pdbfolder == 0} {
			set pdbfolder $inputfolder
		}
		for {set i 1} {$i <= $bestN} {incr i} {
			set folder mdff_${jobname}_$i
			set prefix ${jobname}_${mol}_best$i
			exec mkdir -p $folder

			exec cp rosetta_output_$jobname/$prefix.pdb $folder/$prefix.pdb

			exec cp $mapname.dx $folder/
			cd $folder

			# TODO: check iteration variables!
			#arguments: jobname MOL mapname fixedselection gscale minSteps num ch_seg topdir topfile parfile
			for {set j 1} {$j <= [llength $config]} {incr j} {
				set smoothr [lindex $config [expr $j-1] 1]
				puts "Smooth radius $smoothr"
				exec mkdir -p $inputfolder/density_grids
				if {$smoothr == 0} {
						mdff griddx -i $mapname.dx -o $inputfolder/density_grids/griddx$j.dx
						#exec /bin/cp -u $inputfolder/density_grids/griddx0.dx $inputfolder/density_grids/griddx$i.dx
					} else {
						volutil -smooth $smoothr $mapname.dx -o $inputfolder/${mapname}_smooth_$smoothr.dx
						mdff griddx -i $inputfolder/${mapname}_smooth_$smoothr.dx -o $inputfolder/density_grids/griddx$j.dx
					}
			}
			auto_mdff_casc_init $jobname $prefix $mapname $fixedselection $config $minSteps $dcdfreq $ch_seg $mutations $topdir $topfiles $parfiles $inputfolder $gridconfig

			for {set k 1} {$k <= [llength $config]} {incr k} {
				puts "Running NAMD step $k"
				exec $path/namd2 $namdArgs mdff_${prefix}_step$k.namd > mdff_${prefix}_step$k.log
				puts "Finished NAMD step $k"
			}
			puts "Analysis started."
			analyze_mdff $prefix mdff_${prefix}_step[llength $config] $mapname $res $ch_seg
			puts "Analysis finished."

			set outname mdff_${prefix}_step[llength $config].dcd-last

			exec sed -i -e {s/HSD/HIS/g} $outname.pdb
			exec sed -i -e {s/HSE/HIS/g} $outname.pdb
			exec sed -i -e {s/HSP/HIS/g} $outname.pdb

			#set beta & occupancy to 1.0
			mol delete all
			mol new $outname.pdb
			[atomselect top all] set beta 1.0
			[atomselect top all] set occupancy 1.0
			[atomselect top all] writepdb $outname.pdb
			mol delete all

			if {$score != 0} {
				exec $rosettapath/score.$platform -database $rosettaDBpath -out:file:scorefile ${jobname}_$i-mdff.sc -ignore_zero_occupancy false -in:file:s $outname.pdb > ${jobname}_$i-scoring.log
				lappend scores [list $i [readscorefile ${jobname}_$i-mdff.sc]]
			}
			exec cp $outname.pdb ../full_length_model/$prefix-mdff.pdb
			cd ..
		}
	}
	return $scores
}

proc helix_reg {run MOL seltext angleStep xLowerLimit xUpperLimit xStep Density Res cutoff Threshold inputfolder config {writetofulllengthmodel 0}} \
{
	global ch_seg
	global topdir
	global topfiles
	exec mkdir -p register_helix_$run
	# exec cp regression.py $inputfolder
	# exec cp run_python_reg.sh $inputfolder
	cd register_helix_$run
	# exec chmod +x $inputfolder/run_python_reg.sh

	mol new $inputfolder/$MOL.pdb
	set sel [atomselect top "$seltext"]
	$sel writepdb $inputfolder/$MOL-$run.pdb
	# volmap mask $sel -o $inputfolder/mask.dx -cutoff $cutoff
	# volutil -mult $inputfolder/$Density.dx $inputfolder/mask.dx -o $inputfolder/$Density-$run.dx
	mol delete all
	set prefix $MOL
	set densprefix $Density-$run

	set current [pwd]
	set files [helixRotateMove $MOL $seltext $angleStep $xLowerLimit $xUpperLimit $xStep $Density $Res $Threshold $inputfolder]
	puts $files
	set cc_list []
	foreach pdb $files {
		cd ${current}/$pdb
		start_mdff_run $run $pdb human_0308_4.1_density-test "none" 0.3 300 0 4.1 0 0 1 $config [pwd]
		set f [open "${current}/$pdb/mdff_$run/last_ccc.txt" r]
		set d [read $f]
		close $f
		set cc [lindex [split $d "\n"] 0]
		lappend cc_list $cc
	}
	cd ..

	puts $cc_list
	set best [::tcl::mathfunc::max {*}$cc_list]
	set best_idx [lsearch -exact $cc_list $best]
	puts "best helix is [lindex $files $best_idx] with a cc of: $best"
	set best_pdb [lindex $files $best_idx]

	exec cp $best_pdb/mdff_$run/mdff_${run}_${best_pdb}_step[llength $config].dcd-last.pdb ../full_length_model/


	cd ..
	return "mdff_${run}_${best_pdb}_step[llength $config].dcd-last"
}

proc readscorefile {fname} \
{
	set f [open $fname r]
	set data [read $f]
	close $f
	set lines [split $data "\n"]
	set score [lindex $lines 1 1]
	return $score
}

## make_dx_file
# convert mrcfile to dx file
proc make_dx_file {mrcfile} \
{
	mdff griddx -i $mrcfile.mrc -o ${mrcfile}_potential.dx
	mdff griddx -i ${mrcfile}_potential.dx -o ${mrcfile}.dx
}

proc make_mrc_file {dxfile} \
{
	global packagePath
	volutil ${dxfile}.dx -o ${dxfile}.situs
	exec $packagePath/run_situs.sh ${dxfile}.situs ${dxfile}.mrc
}

proc smooth_density {smoothr mapname} \
{
	volutil -smooth $smoothr $mapname.dx -o ${mapname}_smooth_$smoothr.dx
}


proc write_phenixpdb {filename seltext} {
  set ml [mol new $filename.pdb]
  set sel [atomselect $ml "$seltext"]
  $sel set occupancy 1
  $sel writepdb $filename-phenix.pdb

  set frpdb [open $filename "r"]
  set spdb [read $frpdb]
  close $frpdb
  set fwpdb [open $filename "w"]
  regsub -all "HSD" $spdb "HIS" spdb
  regsub -all "HSE" $spdb "HIS" spdb
  regsub -all "URA" $spdb "  U" spdb
  regsub -all "ADE" $spdb "  A" spdb
  regsub -all "CYT" $spdb "  C" spdb
  regsub -all "GUA" $spdb "  G" spdb
  regsub -all "THY" $spdb "  T" spdb
  regsub -all "CYN" $spdb "CYS" spdb
  regsub -all -line {^.*CRYST.*$} $spdb " " spdb
  puts $fwpdb $spdb
  close $fwpdb
}
