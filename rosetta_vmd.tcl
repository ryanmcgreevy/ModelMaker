####################################################
# VMD Rosetta Scripts Generator
# Maximilian Scheurer, March 2016
####################################################


# JOB CONTROLLER
#### DO NOT CHANGE ANYTTHING HERE
package require mdff

package require RosettaInputGenerator
package require RosettaScoring
package require CCColor
package require FindSelection
package require RosettaUtilities
package require SSAnalysis
package require topotools 

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
	namespace export write_rosettapdb

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
		set find_sel [find_selection $::MODELMAKER::workdir/setup-$jobname/$mol $selections $find_cfg 5]
	} else {
		set find_sel [find_selection $::MODELMAKER::workdir/setup-$jobname/$mol $selections $find_cfg 0]
	}
	set spans [lindex $find_sel 0]
	set exclude [lindex $find_sel 1]
	set chains [lindex $find_sel 2]
	set constraints [lindex $find_sel 3]

  set anchor [get_anchor_residue $::MODELMAKER::workdir/setup-$jobname/$mol $anchor]

	set ros_config [list $chains $spans $exclude [list $anchor $constraints]]

	######################
	# Cleanup input file
	write_rosettapdb $::MODELMAKER::workdir/setup-$jobname/$mol.pdb
  ######################
	puts "Rosetta started"
	refine_with_rosetta $jobname $mol.pdb $mapname $mapresolution $score_dens $nstruct $cluster $nPerTask $ros_config $cartesian
	file attributes $jobname.sh -permissions +x
	file rename $jobname.sh $::MODELMAKER::workdir/run-$jobname/
	file copy -force $mapname.mrc $::MODELMAKER::workdir/setup-$jobname/

	if {!$scoreOnly} {
		file mkdir $::MODELMAKER::workdir/run-$jobname/sc_out
		file mkdir $::MODELMAKER::workdir/run-$jobname/pdb_out
		file mkdir $::MODELMAKER::workdir/run-$jobname/OUTPUT_FILES
    set logfile [open $::MODELMAKER::workdir/run-$jobname/rosetta_log_$jobname.log a]
    set script [open "| $::MODELMAKER::workdir/run-$jobname/$jobname.sh $jobname $mol.pdb 2>@stderr &" r]
    set numfiles 0
    while {[gets $script line] >= 0} {
      puts $logfile $line
      set current [llength [glob -nocomplain $::MODELMAKER::workdir/run-$jobname/pdb_out/*.pdb ] ]
      if { $current > $numfiles} {
        puts "Current number of files: $current - [expr double($current)/($nstruct) * 100.0] % complete"
        set numfiles $current
      }
    }

    close $logfile
    puts "Rosetta finished."
	
  }

	# Scoring
  puts "Scoring normal run."
	score_refinement ${jobname}_$mol $bestN $jobname
  
  foreach pdb [glob -nocomplain $::MODELMAKER::workdir/run-$jobname/pdb_out/*.pdb ] {
    set number [string trimleft [file tail $pdb] "${jobname}_${mol}"]
    set newname "${jobname}_${number}"
    file rename -force $pdb  [file join $::MODELMAKER::workdir/run-$jobname/pdb_out/ $newname]
  }


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
	set find_sel [find_selection $::MODELMAKER::workdir/setup-$jobname/$mol $selections $find_cfg 0]
	set spans [lindex $find_sel 0]
	set exclude [lindex $find_sel 1]
	set chains [lindex $find_sel 2]
	set constraints [lindex $find_sel 3]

	set anchor [get_anchor_residue $::MODELMAKER::workdir/setup-$jobname/$mol $anchor]
	#puts "anchor residue converted: $anchor"

	set ros_config [list $chains $spans $exclude [list $anchor $constraints]]
	######################
	# Cleanup input file
	write_rosettapdb $::MODELMAKER::workdir/setup-$jobname/$mol.pdb
	######################
  puts "Rosetta sidechain started"
	refine_with_rosetta $jobname $mol.pdb $mapname $mapresolution $score_dens $nstruct $cluster $nPerTask $ros_config
	file attributes $jobname.sh -permissions +x
  file rename $jobname.sh $::MODELMAKER::workdir/run-$jobname/
	file copy -force $mapname.mrc $::MODELMAKER::workdir/setup-$jobname/

	if {!$scoreOnly} {
		file mkdir $::MODELMAKER::workdir/run-$jobname/sc_out
		file mkdir $::MODELMAKER::workdir/run-$jobname/pdb_out
		file mkdir $::MODELMAKER::workdir/run-$jobname/OUTPUT_FILES
    set logfile [open $::MODELMAKER::workdir/run-$jobname/rosetta_log_$jobname.log a]
    set script [open "| $::MODELMAKER::workdir/run-$jobname/$jobname.sh $jobname $mol.pdb 2>@stderr &" r]
    set numfiles 0
    while {[gets $script line] >= 0} {
      puts $logfile $line
      set current [llength [glob -nocomplain $::MODELMAKER::workdir/run-$jobname/pdb_out/*.pdb ] ]
      if { $current > $numfiles} {
        puts "Current number of files: $current - [expr double($current)/($nstruct) * 100.0] % complete"
        set numfiles $current
      }
    }

    close $logfile
    puts "Rosetta finished."
	}

  puts "Scoring normal run."
  score_refinement ${jobname}_$mol $bestN $jobname
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

	set anchor [get_anchor_residue $::MODELMAKER::workdir/setup-$jobname/$mol $anchor]

	set ros_config [list $chains $spans $exclude [list $anchor $constraints]]

	######################
	# Cleanup input file
	write_rosettapdb full_length_model/$mol.pdb
	######################
	puts "Rosetta basic refinement (without density) started"
	rosetta_basic_refinement $jobname $mol.pdb $nstruct $cluster $nPerTask $ros_config

	exec chmod +x $jobname.sh
	exec mv $jobname.sh rosetta_output_$jobname/
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
		}
		puts $output
		puts "Rosetta finished"
	}

	# Scoring
  puts "Scoring normal run."
  score_refinement ${jobname}_$mol $bestN
	cd ..
}



proc start_rosetta_abinitio {jobname mol selections anchor fragfiles nstruct {cluster 0} {nPerTask 5} {testrun 0} args} \
{
	set selection_length [llength $selections]
	set find_cfg []
	for {set i 0} {$i < $selection_length} {incr i} {
		# allow backbone and sidechain movement for selections
		lappend find_cfg {1 1}
	}

	# Cleanup input file
	write_rosettapdb $::MODELMAKER::workdir/setup-$jobname/$mol.pdb


	# MOL selections config {offset 4}
	set find_sel [find_selection $::MODELMAKER::workdir/setup-$jobname/$mol $selections $find_cfg 0]
	set spans [lindex $find_sel 0]
	set exclude [lindex $find_sel 1]
	set chains [lindex $find_sel 2]
	set constraints [lindex $find_sel 3]

	set anchor [get_anchor_residue $::MODELMAKER::workdir/setup-$jobname/$mol $anchor]

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

	file mkdir $::MODELMAKER::workdir/run-$jobname/sc_out
	file mkdir $::MODELMAKER::workdir/run-$jobname/pdb_out
	file mkdir $::MODELMAKER::workdir/run-$jobname/OUTPUT_FILES

  if {[glob -nocomplain $::MODELMAKER::workdir/run-$jobname/pdb_out/*.pdb] != ""} {
    set prefix "tmp"
  } else { 
    set prefix ""  
  }
	
  set username $::MODELMAKER::settings(username)
	puts "Rosetta abinitio started."
	rosetta_abinitio $jobname $mol.pdb $fragfiles $nstruct $cluster $nPerTask $testrun $ros_config $chain_idents $prefix
	file attributes $::MODELMAKER::workdir/run-$jobname/$jobname.sh -permissions +x


	set logfile [open $::MODELMAKER::workdir/run-$jobname/rosetta_log_$jobname.log a]
  set script [open "| $::MODELMAKER::workdir/run-$jobname/$jobname.sh $jobname $mol.pdb 2>@stderr &" r]
	set numfiles 0
  while {[gets $script line] >= 0} {
		puts $logfile $line
		set current [llength [glob -nocomplain $::MODELMAKER::workdir/run-$jobname/pdb_out/${prefix}${jobname}*.pdb ] ]
    if { $current > $numfiles} {
      puts "Current number of files: $current - [expr double($current)/($nstruct) * 100.0] % complete"
      set numfiles $current
    }
	}
	close $logfile
  
#  #having problems forcing rosetta to keep sidechains fixed, so we will just replace the supposedly fixed
#  #parts with the input template structure
  set inpmol [mol new $mol.pdb]
  set fixsel [atomselect $inpmol "not ($selections)"]
  foreach struct [glob -nocomplain $::MODELMAKER::workdir/run-$jobname/pdb_out/${prefix}${jobname}*.pdb ] {
    set mymol [mol new $struct]
    set fitsel [atomselect $mymol "not ($selections)"]
    set mat [measure fit $fixsel $fitsel]
    $fixsel move $mat

    set mysel [atomselect $mymol $selections]
    set mol3 [::TopoTools::selections2mol "$fixsel $mysel"] 
    animate write pdb $struct $mol3 
  }

  if {$prefix != ""} {
    set nstruct [llength [glob -nocomplain $::MODELMAKER::workdir/run-$jobname/pdb_out/${jobname}*.pdb ]]
    if {$nstruct < 10000} {
      set offset 4
    } else {
      set offset [expr int(floor(log10($end))) + 1] 
    }

    set last_pdb [lindex [lsort -decreasing [glob -nocomplain $::MODELMAKER::workdir/run-$jobname/pdb_out/${jobname}*.pdb ]] 0 ]
    set last_num [lindex [split [file rootname $last_pdb] "_"] end]
    
    set i 1
    foreach tmppdb [lsort -increasing [glob -nocomplain $::MODELMAKER::workdir/run-$jobname/pdb_out/${prefix}${jobname}*.pdb  ]] {
     set newname $::MODELMAKER::workdir/run-$jobname/pdb_out/${jobname}_${mol}_[format %0${offset}i  [expr $last_num + $i]].pdb
     file rename -force $tmppdb $newname     
     foreach scorefile [glob -nocomplain $::MODELMAKER::workdir/run-$jobname/sc_out/*] {
       set frpdb [open $scorefile "r"]
       set spdb [read $frpdb]
       close $frpdb
       set fwpdb [open $scorefile "w"]
       regsub -all [file rootname [file tail $tmppdb]] $spdb [file rootname [file tail $newname]] spdb
       puts $fwpdb $spdb
       close $fwpdb
     }
     incr i
   }
   set orig_score [open $::MODELMAKER::workdir/run-$jobname/sc_out/${jobname}_score.sc "a"]
   set temp_score [open $::MODELMAKER::workdir/run-$jobname/sc_out/${prefix}${jobname}_score.sc "r"]
   set lines [split [read $temp_score] \n]
   for {set j 2} {$j < [llength $lines]} {incr j} {
     puts $orig_score [lindex $lines $j]
   }
   close $orig_score
   close $temp_score
   file delete -force $::MODELMAKER::workdir/run-$jobname/sc_out/${prefix}${jobname}_score.sc
  }
  puts "Rosetta abinitio finished."
}


proc start_rosetta_insertion {jobname mol selections fragfiles fasta nstruct {cluster 0} {nPerTask 5} args} \
{
	# prepare configuration
	set selection_length [llength $selections]
	set find_cfg []
	for {set i 0} {$i < $selection_length} {incr i} {
		# allow backbone and sidechain movement for selections
		lappend find_cfg {1 1}
	}

	# Cleanup input file
	write_rosettapdb $::MODELMAKER::workdir/setup-$jobname/$mol.pdb


	# MOL selections config {offset 4}
	set find_sel [find_selection $::MODELMAKER::workdir/setup-$jobname/$mol $selections $find_cfg 0]
	set spans [lindex $find_sel 0]
	set exclude [lindex $find_sel 1]
	set chains [lindex $find_sel 2]

	set ros_config [list $chains $spans $exclude]

	puts "Rosetta insertion folding started."
	rosetta_insertion $jobname $mol $fragfiles $fasta $nstruct $cluster $nPerTask $ros_config
	file attributes $::MODELMAKER::workdir/run-$jobname/$jobname.sh -permissions +x

	file mkdir $::MODELMAKER::workdir/run-$jobname/sc_out
	file mkdir $::MODELMAKER::workdir/run-$jobname/pdb_out
	file mkdir $::MODELMAKER::workdir/run-$jobname/OUTPUT_FILES
  file mkdir $::MODELMAKER::workdir/run-$jobname/intermediates
  set logfile [open $::MODELMAKER::workdir/run-$jobname/rosetta_log_$jobname.log a]
  set script [open "| $::MODELMAKER::workdir/run-$jobname/$jobname.sh 2>@stderr &" r]
  set numfiles 0
  while {[gets $script line] >= 0} {
    puts $logfile $line
    set current [llength [glob -nocomplain $::MODELMAKER::workdir/run-$jobname/pdb_out/$jobname*.pdb ] ]
    if { $current > $numfiles} {
      puts "Current number of files: $current - [expr double($current)/($nstruct) * 100.0] % complete"
      set numfiles $current
    }
  }
	file rename {*}[glob $::MODELMAKER::workdir/run-$jobname/pdb_out/loops_closed*.pdb] $::MODELMAKER::workdir/run-$jobname/intermediates/

	puts "Rosetta insertion folding finished."
}

proc analyze_abinitio {jobname mol template bestN nstruct cluster align_template align_rosetta analysis_components kmin kmax {insertion 0} args} \
{
	global vmdexe
	global packagePath
	cd $::MODELMAKER::workdir/run-$jobname

	puts "Ab-initio analysis started."

	# ALIGNMENT
	puts "current directory: [pwd]"
	file mkdir pdb_out_aligned

#originally this disctinction was made when cluster used to indicate Juan's code. Not sure why this was needed...
#	if {!$cluster} {
		if {$insertion != 0} {
			align_rosetta_local 1 $nstruct ${jobname}_$mol $template $insertion $align_template $align_rosetta
		} else {
 		  align_rosetta_local 1 $nstruct ${jobname}_$mol $template $mol $align_template $align_rosetta
		}
#	} else {
#		if {$insertion != 0} {
#			align_rosetta_cluster 1 $nstruct ${jobname}_$mol $template $insertion $align_template $align_rosetta
#		} else {
#			align_rosetta_cluster 1 $nstruct ${jobname}_$mol $template $mol $align_template $align_rosetta
#		}
#	}

	#SCORING MOL max_structures cluster?
	# extra variable to run scoring on cluster folded insertion runs
	set extra 0
	if {$insertion != 0} {
		set extra $mol
	}

	score_abinitio $jobname ${jobname}_$mol $bestN $cluster $extra

	# ANALYSIS
	file mkdir analysis
	file copy -force ${jobname}_${mol}_rosetta_scoring_min_$bestN.dcd analysis
	file copy -force ${jobname}_${mol}_rosetta_scoring_min_$bestN.pdb analysis
	cd analysis

	set prefix ${jobname}_${mol}

  foreach ana $analysis_components {
		set type [lindex $ana 0]
		switch -exact -- $type {
			cluster {
        puts "running cluster analysis"
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
					run_clustering $prefix $start $end $bestN $kmin $kmax
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

				puts "running SS analysis"
		    ss_analysis -molfile ${MOL}_rosetta_scoring_min_${bestN}.pdb -sel "resid $resid_start to $resid_stop and chain $chain_id" -avg 1 -output "ss_${bestN}_${resid_start}_${resid_stop}"
        #evaluate_ss_analysis $bestN $start $end
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

proc smooth_density {smoothr mapname} \
{
	volutil -smooth $smoothr $mapname.dx -o ${mapname}_smooth_$smoothr.dx
}


proc write_rosettapdb {filename} {
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
  #I believe this line is fine for Rosetta
  #regsub -all -line {^.*CRYST.*$} $spdb " " spdb
  puts $fwpdb $spdb
  close $fwpdb
}
