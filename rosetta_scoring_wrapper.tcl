### Wrapper for rosetta scoring
### Maximilian Scheurer, April 2016
### mail: mscheurer@ks.uiuc.edu

namespace eval ::RosettaScoring {
  namespace export score_refinement
  namespace export score_refinement_cluster
  namespace export score_abinitio
  # Set up Variable
  set version 0.1
  set packageDescription "RosettaVMD Plugin"
}
package provide RosettaScoring $RosettaScoring::version

proc score_refinement {args} \
{
  return [eval ::RosettaScoring::score_refinement $args]
}

proc score_refinement_cluster {args} \
{
  return [eval ::RosettaScoring::score_refinement_cluster $args]
}

proc score_abinitio {args} \
{
  return [eval ::RosettaScoring::score_abinitio $args]
}

proc ::RosettaScoring::score_refinement {MOL max_structures jobname} \
{
  #set MOL loopmodel_73-78_4OCM_B.pdb_full_length
  #set max_structures 25
  puts "scoring started"
  set pdb [::RosettaScoring::rosetta_scoring $max_structures $jobname]
  puts $pdb
  puts "scoring finished"
  mol delete all
  #mol new ./pdb_out/${MOL}_[format %04i [lindex $pdb 0]].pdb waitfor all
  for {set i 0} {$i < [llength $pdb]} {incr i} {
   if { [llength $pdb] < 10000} {
    mol new $::MODELMAKER::workdir/run-$jobname/pdb_out/${MOL}_[format %04i [lindex $pdb $i]].pdb
   } else {
    set offset [expr int(floor(log10([llength $pdb]))) + 1] 
    mol new $::MODELMAKER::workdir/run-$jobname/pdb_out/${MOL}_[format %0${offset}i [lindex $pdb $i]].pdb
   }
   [atomselect top all] writepdb $::MODELMAKER::workdir/run-$jobname/${jobname}_best[expr $i + 1].pdb
 }

 #animate write pdb ${MOL}_rosetta_scoring_min_$max_structures.pdb beg 0 end [expr $max_structures-1] skip 1 0

}


proc ::RosettaScoring::score_refinement_cluster {run MOL max_structures} \
{
  set pdb [::RosettaScoring::rosetta_scoring_cluster_consistent $max_structures $run]

  puts $pdb
  for {set i 0} {$i < [llength $pdb]} {incr i} {
    mol new ./pdb_out/${MOL}[lindex $pdb $i]_0001.pdb
    [atomselect top all] writepdb ${MOL}_best[expr $i + 1].pdb
  }
}

proc ::RosettaScoring::score_abinitio {run MOL max_structures cluster {extra 0} args} \
{
  puts $MOL
#not sure why there is a need now for this difference with new cluster script...testing using the same
#  if {$cluster} {
#    set pdb [::RosettaScoring::rosetta_scoring_cluster_consistent $max_structures $run $extra]
#    } else {
      set pdb [::RosettaScoring::rosetta_scoring $max_structures $run]
#    }
    

    #set pdblist [lsort -dictionary [glob ../rosetta_output/pdb_out_aligned/*.pdb]]
    #set file "../rosetta_output/pdb_out/${MOL}[lindex $pdb 0]_0001.pdb"
    puts $pdb
    #set file "../rosetta_output_bbfix/pdb_out/${MOL}[lindex $pdb 0]_0001.pdb"
    set file "$::MODELMAKER::workdir/run-$run/pdb_out_aligned/${MOL}_[lindex $pdb 0]_0001_aligned.pdb"
    #set file "../rosetta_output_bbfree/pdb_out_aligned/${MOL}[lindex $pdb 0]_0001_aligned.pdb"
    mol new $file waitfor all
    set sel [atomselect top "all and noh"]
    puts "[file root $file]_noh.pdb"
    $sel writepdb "[file root $file]_noh.pdb"
    $sel delete
    mol delete [molinfo top]
    set top [mol new "[file root $file]_noh.pdb"]
    file delete "[file root $file]_noh.pdb"
    for {set i 1} {$i < [llength $pdb]} {incr i} {
      #set file "../rosetta_output/pdb_out/${MOL}[lindex $pdb $i]_0001.pdb"
      #set file "../rosetta_output_bbfix/pdb_out/${MOL}[lindex $pdb $i]_0001.pdb"
      set file "$::MODELMAKER::workdir/run-$run/pdb_out_aligned/${MOL}_[lindex $pdb $i]_0001_aligned.pdb"
      #set file "../rosetta_output_bbfree/pdb_out_aligned/${MOL}[lindex $pdb $i]_0001_aligned.pdb"
      mol new $file waitfor all
      set sel [atomselect top "all and noh"]
      $sel writepdb "[file root $file]_noh.pdb"
      $sel delete
      mol delete [molinfo top]
      mol addfile "[file root $file]_noh.pdb" waitfor all $top
      file delete "[file root $file]_noh.pdb"
      update idletasks
    }

    animate write pdb ${MOL}_rosetta_scoring_min_$max_structures.pdb beg 0 end [expr $max_structures-1] skip 1 $top
    update idletasks
    update
    animate write dcd ${MOL}_rosetta_scoring_min_$max_structures.dcd beg 0 end [expr $max_structures-1] skip 1 $top
    update idletasks
    update
}



proc ::RosettaScoring::rosetta_scoring {length_tot jobname} {
  set file_list ""
  set val_list ""

  set file_list [lsort -dictionary [glob $::MODELMAKER::workdir/run-$jobname/sc_out/*.sc]]
  set length 0
  for {set i 0} {$i < [llength $file_list]} {incr i} {
    set file [open [lindex $file_list $i] r]
    
    set log [open "$::MODELMAKER::workdir/run-$jobname/rosetta_scoring.log" w+]
    while {[eof $file] != 1} {
       set line [gets $file]
       if {[lindex $line 0] == "SCORE:" && [regsub -all {[^0-9]} [lindex $line 1] ""] != ""} {
          lappend val_list [lindex $line 1]
          incr length
       }
    }
    close $file
  }
  # for {set i 0} {$i < $length} {incr i} {
  #   set file [open [lindex $file_list $i] r]
  #   set rfile [read $file]
  #   set line [lindex [split $rfile "\n"] 2]

    
  #   if {$line == ""} {
  #     set line [lindex [split $rfile "\n"] 0]

  #   }
  #   set val [lindex $line 1]
    
    
  #   #puts "$i FILES LIST [lindex $file_list $i]"
  #   #puts "$i VALUES LIST [lindex $val_list $i]"
  #   close $file
  #}

  set index_sorted [lsort -real -indices $val_list]

  #puts "\n$index_sorted"
  set number ""
  for {set i 0} {$i < [llength $index_sorted]} {incr i} {
    lappend number [expr [lindex $index_sorted $i] +1 ]
    puts $log "[expr [lindex $index_sorted $i] +1 ]\t[lindex $file_list [lindex $index_sorted $i]]\t[lindex $val_list [lindex $index_sorted $i]] "
    
  }
  puts "LENGTH $length"
  if {$length >= $length_tot} {
     return [lrange $number 0 [expr $length_tot -1]]
  } else {
    puts " Not enough file for chosen length_tot total file = $length: length_tot= $length_tot "
  }
 
} 


proc ::RosettaScoring::rosetta_scoring_cluster {length_tot} {

    set file_list [lsort -dictionary [glob ./sc_out/*.sc]]
    #set file_list [lsort -dictionary [glob ../rosetta_output_bbfix/sc_out/*.sc]]
    #set file_list [lsort -dictionary [glob ../rosetta_output_bbfree/sc_out/*.sc]]

    set length [llength $file_list]
    set log [open "rosetta_scoring.log" w+]
    for {set i 0} {$i < $length} {incr i} {
      set file [open [lindex $file_list $i] r]
      set rfile [read $file]
      set line [lindex [split $rfile "\n"] 2]


      if {$line == ""} {
        set line [lindex [split $rfile "\n"] 0]

      }
      set val [lindex $line 1]

      lappend val_list $val
      #puts "$i FILES LIST [lindex $file_list $i]"
      #puts "$i VALUES LIST [lindex $val_list $i]"
      close $file
    }

    set index_sorted [lsort -real -indices $val_list]

    #puts "\n$index_sorted"
    set number ""
    for {set i 0} {$i < [llength $index_sorted]} {incr i} {
      lappend number [expr [lindex $index_sorted $i] +1 ]
      puts $log "[expr [lindex $index_sorted $i] +1 ]\t[lindex $file_list [lindex $index_sorted $i]]\t[lindex $val_list [lindex $index_sorted $i]] "

    }

    if {$length >= $length_tot} {
     return [lrange $number 0 [expr $length_tot -1]]
     } else {
      puts " Not enough file for chosen length_tot total file = $length: length_tot= $length_tot "
    }

  }


#package require Tcl 8.5
proc ::RosettaScoring::rosetta_scoring_cluster_consistent {length_tot runname {extra 0} args} {
  set file_list ""
  set val_list ""
  set id_list ""
  # set ic_val_list ""
  if {$extra == 0} {
    set score_name "${runname}_score" 
  } else {
    set score_name "${runname}_${extra}_score"
  }
  puts $score_name
  set file_list [lsort -dictionary [glob ./sc_out/*.sc]]
  #set file_list [lsort -dictionary [glob ../rosetta_output_bbfix/sc_out/*.sc]]
  #set file_list [lsort -dictionary [glob ../rosetta_output_bbfree/sc_out/*.sc]]

  set file_id_list [] ;# for total score
  # set file_id_if_list [] ;# for interface score

  foreach file $file_list {
    puts "file: $file"
    #regexp -nocase [subst -nocommands -nobackslashes {${score_name}[0-9]*}] $file a 
    set a [glob ./pdb_out/*]
    puts "a: $a"
    # set ident [string trim $a $score_name]
    set sc_length [string length $score_name]
    set ident [string range $a $sc_length end]
    puts "ident: $ident"

    set f [open "$file" r]
    set rfile [read $f]
    set line [lindex [split $rfile "\n"] 2]
    if {$line == ""} {
      set line [lindex [split $rfile "\n"] 0]

    }
    set val [lindex $line 1]
    # set icval [lindex $line 4]
    close $f
    lappend file_id_list [list $ident $file $val]
    # lappend file_id_if_list [list $ident $file $icval]
  }

  set sort_scores [lsort -real -index 2 $file_id_list]
  # set sort_if_scores [lsort -real -index 2 $file_id_if_list]

  puts "finished reading and sorting."

  set length [llength $file_list]
  set log [open "rosetta_scoring.log" w+]
  # set ic_log [open "rosetta_if_sc.log" w+]

  set counter 0
  foreach item $sort_scores {
    if {$counter < $length_tot} {
      lappend number [lindex $item 0]
    }
    puts $log "[lindex $item 0]\t[lindex $item 1]\t[lindex $item 2]"
    incr counter
  }

  set counter 0
  # foreach item $sort_if_scores {
  #   if {$counter < $length_tot} {
  #     lappend ic_number [lindex $item 0]
  #   }
  #   puts $ic_log "[lindex $item 0]\t[lindex $item 1]\t[lindex $item 2]"
  #   incr counter
  # }


  if {$length > $length_tot} {
    # return [list [lrange $number 0 [expr $length_tot -1]] [lrange $ic_number 0 [expr $length_tot -1]]]
    return [lrange $number 0 [expr $length_tot -1]]
    } else {
      puts " Not enough file for chosen length_tot total file = $length: length_tot= $length_tot "
    }

  } 
