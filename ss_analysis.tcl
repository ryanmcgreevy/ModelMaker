# SSAnalysis
#
# Author:
#   João V. Ribeiro
#   Beckman Institute for Advanced Science and Technology
#   University of Illinois, Urbana-Champaign
#   jribeiro@ks.uiuc.edu
#   http://www.ks.uiuc.edu/~jribeiro
#
#  
# Intro:
#   Analyze the secondary structure from one or more frames, returning 
#   
#  
source /Projects/jribeiro/vmd/plugins/multiplot/multiplot.tcl
package require qwikmd
package require multiplot
 
namespace eval ::SSAnalysis {
  set version 0.1
  set description "ss analysis"
}
package provide SSAnalysis $SSAnalysis::version

proc ss_analysis { args } \
{
  return [eval ::SSAnalysis::ss_analysis $args]
}

## Calling to the SSAnalysis package to analyze
## secondary structures of a sequence of frames
## Return: the persistence of each secondary structure 
##      * per residue (in percentage) and the representative frame of 
##        the most prevalent secondary structure (average ss)
##  
##      * dcd (and psf) file containing frames with the average ss
##
##      * dcd (and psf) file containing frames with the ss seq given 
##        by the user

### ssanalysis: secondary structure analysis"
### Usage: ssanalysis \[args...\]"
### Args:"
###   -workdir    <working/project directory for job> (default: current folder)>"
###   -mol        (optional if -molfile provided) VMD molecule containing the frames to be analyzed"
###   -molfile    (optional if -mol provided) file containing the frames to be analyzed"
###   -strtcfile  (mandatory file if -molfile provided) structure file containing structure
###                information (e.g. *.psf or *.pdb)"
###   -sel        atom selection text defining region to be analyzed (default:all )"
###   -output     prefix to be used in the naming of the output files (default: output)"
###   -seq    (optional)search ss patter e.g. HHH - three consecutive residues presenting
###                alpha helix as secondary structure. the sequence's length has to match the 
###                number of residues defined in the atom selection"


proc ::SSAnalysis::ss_analysis { args } {

  set workdir [pwd]
  set mol ""
  set molfile ""
  set strtcfile ""
  set seltext "all"
  set output "output"
  set seq ""
  foreach {indarg val} $args {
    switch -exact -- $indarg {
         -workdir {
            set workdir $val
         }
         -mol {
            set mol $val
         }
         -molfile {
            set molfile $val
         }
         -strtcfile {
            set strtcfile $val
         }
         -outputname {
            set outputname $val
         }
         -sel {
            set seltext $val
         }
         -output {
            set output $val
         }
         -seq {
            set seq $val
         }
         default {
            puts "$indarg argument not defined"
         }
    } 
  }

  if {$molfile != ""} {
    if {$strtcfile == "" && [file extension $molfile] != ".pdb"} {
      puts "SSanalysis: please provide the structure file correspondent to the frame sequence $molfile"
      return
    }
    if {$mol != ""} {
      puts "SSanalysis: please use either mol or the pair \"molfile + strtcfile\" to defined the molecule\
      to be analyzed. Not both"
      return
    }
  }
  
  

  
  ### moving to the work directory and store the current location
  set pwd [pwd]
  cd ${workdir}

  ### open file to write the output containing the ss prevalence per residue (in percentage)
  set f [open "${output}.txt" "w+"]

  puts $f "[string repeat "#" 80]"
  puts $f "##"
  puts $f "## Secondary Structure analysis version $SSAnalysis::version"
  puts $f "##"
  puts $f "[string repeat "#" 80]"
  puts $f "Command executed: ssanalysis $args\n\n"

  

  ### loading the predicted structure to VMD and creating the defined atom selections
  if {$mol == ""} {
    if {[file extension $molfile] == ".pdb"} {
      set mol [mol new $molfile waitfor all]
    } else {
      set mol [mol new $strtcfile waitfor all]
      mol addfile $molfile waitfor all  
    }
    
  } elseif {[string match "*original_mol*" [file root [molinfo $mol get name]]] == 0} {
    set sel [atomselect $mol "all"]
    set oriname [file root [molinfo $mol get name]]
    $sel writepsf ${oriname}_original_mol.psf
    animate write dcd ${oriname}_original_mol.dcd beg 0 end [expr [molinfo $mol numframes] -1] skip 1 $mol waitfor all
    puts $f "Original molecule saved as psf and dcd files with the name: ${oriname}_original_mol\n\n"
    $sel delete
  }
  
  ### Splitting the secondary structure pattern provided by the user into individual list arguments.
  set numframes [molinfo $mol get numframes]
  set sel [atomselect $mol "($seltext) and name CA"]
  if {$seq != ""} {
    if {[string length $seq] != [$sel num]} {
      puts "ERROR: The number of residues defined in -sel doesn't match the length of the -pattern"
      puts "residues number in -sel: [$sel num]"
      puts "seq length: [string length $seq]"
      $sel delete
      return
    } 
  }
  if {[$sel num] == 0} {
    puts "Error!!! No atom selected."
    $sel delete
    return
  }

  ### definition list of secondary structure motives
  set list_H "" ;# Alpha Helix
  set list_T "" ;# Turn
  set list_C "" ;# Coil
  set list_G "" ;# 3-10 Helix
  set list_E "" ;# Extended Conformation
  set list_B "" ;# Bridge
  set list_I "" ;# Pi Helix
  set ssseq ""
  set residues [lsort -unique -real [$sel get resid]]

  ### populate the list with 0's for each residue
  set nres [llength $residues]
  set list_H [string repeat "0 " $nres]
  set list_T [string repeat "0 " $nres]
  set list_C [string repeat "0 " $nres]
  set list_G [string repeat "0 " $nres]
  set list_E [string repeat "0 " $nres]
  set list_B [string repeat "0 " $nres]
  set list_I [string repeat "0 " $nres]


  ### evaluating the occurrence of each secondary structure element of each frame
  for {set i 0} {$i < $numframes} {incr i} {
    animate goto $i
    display update ui
    $sel frame $i
    mol ssrecalc $mol 
    set j 0
    lappend ssseq [$sel get structure]
    foreach val [lindex $ssseq $i] {
      switch $val {   
        "H" {lset list_H $j [expr [lindex $list_H $j] + 1]}
        "T" {lset list_T $j [expr [lindex $list_T $j] + 1]}
        "C" {lset list_C $j [expr [lindex $list_C $j] + 1]}
        "G" {lset list_G $j [expr [lindex $list_G $j] + 1]}
        "E" {lset list_E $j [expr [lindex $list_E $j] + 1]}
        "B" {lset list_B $j [expr [lindex $list_B $j] + 1]}
        "I" {lset list_I $j [expr [lindex $list_I $j] + 1]}
      }
      incr j
    }

  }

  puts $f "Table with secondary structure prevalence (percentage) of each residue in all frames "
  puts $f "Number of structures analyzed: $numframes.\n" 
  puts $f [format {%10s%10s%10s%10s%10s%10s%10s%10s} "Resid" "H  " "T  " "C  " "G  " "E  " "B  " "I  "]
  puts $f [format {%10s%10s%10s%10s%10s%10s%10s%10s} "-----" "-----" "-----" "-----" "-----" "-----" "-----" "-----"]
  set formatNum {%10s%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f}
  ### Iteration over all input structures to compute the occurrence of each secondary
  ### structure element (in percentage) for each residue
  set avg_seq ""
  set plotdata [list]
  for {set i 0} {$i < [llength $residues]} {incr i} {
    set percent_H [expr [lindex $list_H $i]*1.0/$numframes*100.0]
    set percent_T [expr [lindex $list_T $i]*1.0/$numframes*100.0]
    set percent_C [expr [lindex $list_C $i]*1.0/$numframes*100.0]
    set percent_G [expr [lindex $list_G $i]*1.0/$numframes*100.0]
    set percent_E [expr [lindex $list_E $i]*1.0/$numframes*100.0]
    set percent_B [expr [lindex $list_B $i]*1.0/$numframes*100.0]
    set percent_I [expr [lindex $list_I $i]*1.0/$numframes*100.0]
    set val [lsort -indices -decreasing -real [list $percent_H $percent_T $percent_C $percent_G $percent_E $percent_B $percent_I]]
    switch [lindex $val 0] {
      0 {lappend avg_seq H}
      1 {lappend avg_seq T}
      2 {lappend avg_seq C}
      3 {lappend avg_seq G} 
      4 {lappend avg_seq E}
      5 {lappend avg_seq B}
      6 {lappend avg_seq I}

    }
    puts $f [format $formatNum [lindex $residues $i] $percent_H $percent_T $percent_C\
     $percent_G $percent_E $percent_B $percent_I]
    
    ## store data to be plotted
    lappend plotdata [list [lindex $residues $i] $percent_H $percent_T $percent_C\
     $percent_G $percent_E $percent_B $percent_I]
  }
  
  puts $f "Average sequence (residue basis):  $avg_seq\n\n"
  

  ### Plot the secondary structure prevalence per residue 
  set plot [multiplot -xsize 600 -ysize 400 -title "Secondary Structure" -xlabel ResidID -ylabel "Percentage (%)" -nolines -linewidth 2 -bkgcolor \#d1efff \
  -xmin [expr [lindex $residues 0] -0.5] -xmax [expr [lindex $residues end] + 0.5] -ymin 0 -ymax 100]


  set colors [list]
  set ssList [list H T C G E B T]
  set legendList [list "Alpha Helix" "Turn" "Coil" "3-10 Helix" "Extended\nConformation" "Bridge" "Pi Helix"]
  set index 0
  foreach ssname $ssList {
    set hexcols [QWIKMD::chooseColor $ssname]
            
    set hexred [lindex $hexcols 0]
    set hexgreen [lindex $hexcols 1]
    set hexblue [lindex $hexcols 2]
    set color "\#${hexred}${hexgreen}${hexblue}"
    lappend colors $color
    $plot add [expr [lindex $residues 0] - 0.5] 0 -legend [lindex $legendList $index] -linecolor $color
    incr index
  }

  $plot replot
  foreach res $plotdata {
    set residue [lindex $res 0]
    set data [lrange $res 1 end]
    set xmin [expr $residue - 0.5]
    set xmax [expr $residue + 0.5]
    
    set min 0.0
    for {set i 0} {$i < [llength $data]} {incr i} {
      set val [expr [lindex $data $i] + $min]
      $plot draw rectangle $xmin $min $xmax $val -fill [lindex $colors $i] -tag "${residue}ss$i"
      set min $val
    }
  }
  $plot replot

  variable [$plot namespace]::c
  [$plot getpath].f.cf move legend -150 0
  update
  update idletasks
  $c postscript -file ${output}_highest_score.ps 
  
  update
  update idletasks
  

  ### convert image from postscript to gif
  ### uses ImageMagick, only available on unix machines
  ### Windows users can use their own image editor to convert the ps file  
  global tcl_platform env

  if {$::tcl_platform(os) == "Darwin" || $::tcl_platform(os) == "Linux"} {
    eval ::ExecTool::exec "convert ${output}_highest_score.ps gif:./${output}_highest_score.gif"
    update
    update idletasks
  }
  
  # $plot quit

  ### Identify the secondary structure sequence with most prevalence secondary structure
  set score_list ""
  for {set j 0} {$j < $numframes} {incr j} {
    set score_counter 0
    set frm_seq [lindex $ssseq $j]
    for {set k 0} {$k < [llength $frm_seq]} {incr k} {
      if {[lindex $frm_seq $k] == [lindex $avg_seq $k] } {incr score_counter} 
    }
    lappend score_list $score_counter
  }

  ### Identifying the representative structures: meaning the ones representing 
  ### the secondary structure sequence closest to the average secondary structure avg secondary structure above. 
  set final_results_index [lsort -real -indices -decreasing $score_list]
  set final_results [lsort -decreasing -real $score_list]
  set frames [lsearch -all $final_results [lindex $final_results 0]]
  puts $f "Score Table"
  puts $f "Score: Percentage of finding the most prevalent secondary\n structure per residue in that particular frame."
  set formatNum {%10.0f%10.2f}
  puts $f [format "%10s%10s" "Frame" "Score"]
  puts $f [format "%10s%10s" "----" "----"]
  for {set i 0} {$i < [llength $final_results_index]} {incr i} {
    set val [expr [lindex $final_results $i]*1.0/$nres*100.0]
    puts $f [format $formatNum [lindex $final_results_index $i] $val]
  }
  if {[llength $frames] >= 1 && $frames != ""} {
    
    set lframes ""
    for {set i 0} {$i < [llength $frames]} {incr i} {
      set lframes [append lframes "[lindex $final_results_index $i] "]
    }
    puts $f "Closest frame(s) to the average ss sequence: $lframes (N=[llength $frames])"
  } else {
    puts $f "It was not possible to find any frame in the secondary structure [lindex $final_results 0]."
    close $f
    return
  }
  puts $f "\n"
  flush $f


  ### get the frames matching the avg ss sequence and save it
  set lastindexes [expr $numframes - 1]

  ### Save the frames sorted by closiness of the secondary structure to the average secondary structure
  ### on the residue basis
  animate write dcd tmp_sorted_ssseq.dcd beg 0 end [expr $numframes -1] skip 1 $mol
  set sel [atomselect $mol "all" frame [lindex $final_results_index 0]]
  $sel writepdb tmp_sorted_ssseq.pdb
  set newmol [mol new tmp_sorted_ssseq.pdb]
  for {set j 1} {$j < [llength $final_results_index]} {incr j} {
    mol addfile tmp_sorted_ssseq.dcd first [lindex $final_results_index $j] last [lindex $final_results_index $j] waitfor all
  }
  mol bondsrecalc $newmol
  mol ssrecalc $newmol
  
  animate write dcd ${output}_sorted_ssseq_traj.dcd beg 1 end [expr [llength $final_results] -1] skip 1 $newmol
  file rename -force tmp_sorted_ssseq.pdb ${output}_sorted_ssseq_highest_score.pdb
  file delete -force tmp_sorted_ssseq.dcd

  puts $f "Closest structure file: ${output}_sorted_ssseq_highest_score.pdb"
  puts $f "Remaining structures sorted by score: ${output}_sorted_ssseq_traj.dcd"
  ### Search for consecutive alpha-helix seq in the atom selection defined in -sel"
  ### Not sure what to do with this for now
  if {$seq != ""} {
    puts $f "Searching for the secondary structure sequence $seq"
    set searchindex ""
    set frames ""
    for {set i 0} {$i< [llength $ssseq]} {incr i} {
      for {set j 0} {$j < [$sel num]} {incr j} {
        set index [lindex $residrange $j]
        set ini [expr [lindex $index 0] - $resid_start] 
        set end [expr [lindex $index 1] - $resid_start] 
        set pattern [string repeat "H" [expr  $end - $ini +1] ] 
        set strng [join [lrange [lindex $ssseq $i] $ini $end]  ""]
        
        if {$strng == $pattern} {
          append searchindex "frame :[expr $i +1] resid range: $ini - $end\n"
          lappend frames $i 
        }
      }
      
    }
    # puts $f "Pattern search result : $searchindex"
    
    # #Writes all structures containing the consecutive alpha-helix pattern in the range defined above by the variable "residrange" to a pdb file
    # mol new ${MOL}_rosetta_scoring_min_${max_structures}.pdb first [lindex [lindex $frames] 0] last [lindex [lindex $frames] 0] -waitfor all
    
    # for {set j 1} {$j < [llength [lindex $frames]]} {incr j} {
    #   mol addfile ${MOL}_rosetta_scoring_min_${max_structures}.pdb first [lindex [lindex $frames $j]] last [lindex [lindex $frames $j]] -waitfor all
    # }
    # animate write pdb ss_pattern_${max_structures}.pdb beg 0 end [expr [llength $frames] -1] skip 1 top
  }
  close $f
  return
}
