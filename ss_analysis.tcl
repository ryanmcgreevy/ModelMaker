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
package require multiplot
 
namespace eval ::SSAnalysis {
  set version 0.1
  set description "ss analysis"
}

package provide SSAnalysis 0.1

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
##      * dcd (and pdb) file containing frames with the average ss
##
##      * dcd (and pdb) file containing frames with the ss pattern searched 
##        by the user
##

### ssanalysis: secondary structure analysis"
### Usage: ssanalysis \[args...\]"
### Args:"
###   -workdir    <working/project directory for job> (default: current folder)>"
###   -mol        (optional if -molfile provided) VMD molecule containing the frames to be analyzed"
###   -molfile    (optional if -mol provided) file containing the frames to be analyzed"
###   -strtcfile  (mandatory file if -molfile provided as dcd) structure file containing structure
###                information (e.g. *.psf or *.pdb)"
###   -avg        flag to calculate avg secondary structure <0 or 1>"
###   -sel        (mandatory if -avg 1)atom selection text defining region
###               to be analyzed (default:all )"
###   -showplot   show plot with average secondary structure analysis <0 or 1>
###  -seqfind     flag to search for secondary sequence pattern <0 or 1>"
###  -seq         (mandatory if -seqfind 1) nested list of searching ss element and atom selection
###               pairs. Search for consecutive ss elements in the atom selection. 
###               - e.g. {{"H" "resid 1 to 10","C" "resid 15 to 22"}} "
###  -output     prefix to be used in the naming of the output files (default: output)"


proc ::SSAnalysis::ss_analysis { args } {

  set workdir [pwd]
  set mol ""
  set molfile ""
  set strtcfile ""
  set seltext "all"
  set output "output"
  set seq ""
  set avg 0
  set seqfind 0
  set showplot 1
  set plot multiplot
  if {[catch {package require Tk}] == 1} {
    set plot "gnuplot"
  }
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
         -avg {
            set avg $val
         }
         -sel {
            set seltext $val
         }
         -plot {
            set plot $val
         }
         -showplot {
            set showplot $val
         }
         -seqfind {
            set seqfind $val
         }
         -seq {
            set seq $val
         }
         -output {
            set output $val
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

  if {$avg == 0 && $seqfind == 0} {
    puts "SSanalysis: nothing to do here."
    return
  } elseif {$seqfind == 1 && $seq == ""} {
    puts "SSanalysis: please provide the \"ss \[atom selection\]\" pairs for the secondary structure search."
    return
  }
  
  ### moving to the work directory and store the current location
  set pwd [pwd]
  cd ${workdir}

  ### open file to write the output containing the ss prevalence per residue (in percentage)
  set f [open "${output}.txt" "a+"]

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
    
  } elseif {[string match "*original_mol*" [file root [molinfo $mol get name]]] == 0 &&\
  [file exists [file root [molinfo $mol get name]]_original_mol.pdb] == 0} {
    set sel [atomselect $mol "all" frame 0]
    set oriname [file root [molinfo $mol get name]]
    $sel writepdb ${oriname}_original_mol.pdb
    animate write dcd ${oriname}_original_mol.dcd beg 1 end [expr [molinfo $mol get numframes] -1] skip 1 waitfor all $mol
    puts $f "Original molecule saved as pdb and dcd files with the name: ${oriname}_original_mol\n\n"
    $sel delete
    mol top $mol
  }
  set numframes [molinfo $mol get numframes]
  if {$avg == 1} {
    ### Splitting the secondary structure pattern provided by the user into individual list arguments.
    set mainsel [atomselect $mol "($seltext) and name CA"]
    if {[$mainsel num] == 0} {
      puts "Error!!! No atom selected."
      $mainsel delete
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
    set onechain 1; # sign if the selection involves more than one chain

    set residues [lsort -unique -real [$mainsel get residue]]
    set residChain [list]
    foreach resid [$mainsel get resid] chain [$mainsel get chain] {
      lappend residChain "${resid}_${chain}"
    }
    if {[llength [lsort -unique [$mainsel get chain]]] > 1} {
      set onechain 0
    }
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
      $mainsel frame $i
      mol ssrecalc $mol 
      set j 0
      lappend ssseq [$mainsel get structure]
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
    puts $f [format {%15s%15s%10s%10s%10s%10s%10s%10s%10s} "Resid_Chain" "Res_Number" "H  " "T  " "C  " "G  " "E  " "B  " "I  "]
    puts $f [format {%15s%15s%10s%10s%10s%10s%10s%10s%10s} "-----" "-----" "-----" "-----" "-----" "-----" "-----" "-----" "-----"]
    set formatNum {%15s%15s%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f}
    ### Iteration over all input structures to compute the occurrence of each secondary
    ### structure element (in percentage) for each residue
    set avg_seq ""
    set xplotdata [list]
    set yplotdata [list]
    set reslistplot [$mainsel get resid]
    if {$onechain == 0} {
      set reslistplot $residues
    }
    # set resnumber $residues
    for {set i 0} {$i < [llength $residues]} {incr i} {
      set percent_H [format %.1f [expr [lindex $list_H $i]*1.0/$numframes*100.0]]
      set percent_T [format %.1f [expr [lindex $list_T $i]*1.0/$numframes*100.0]]
      set percent_C [format %.1f [expr [lindex $list_C $i]*1.0/$numframes*100.0]]
      set percent_G [format %.1f [expr [lindex $list_G $i]*1.0/$numframes*100.0]]
      set percent_E [format %.1f [expr [lindex $list_E $i]*1.0/$numframes*100.0]]
      set percent_B [format %.1f [expr [lindex $list_B $i]*1.0/$numframes*100.0]]
      set percent_I [format %.1f [expr [lindex $list_I $i]*1.0/$numframes*100.0]]
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

      puts $f [format $formatNum [lindex $residChain $i] [lindex $residues $i] $percent_H $percent_T $percent_C\
       $percent_G $percent_E $percent_B $percent_I]
      
      ## store data to be plotted
      lappend xplotdata [lindex $reslistplot $i] 
      lappend yplotdata [list $percent_H $percent_T $percent_C\
       $percent_G $percent_E $percent_B $percent_I]
    }
    
    puts $f "Average sequence (residue basis):  $avg_seq\n\n"
    

    ### Plot the secondary structure prevalence per residue 
    set xlabel "Residue ID"
    if {$onechain == 0} {
      set xlabel "Residue Number"
    }
    if {[catch {package require Tk}] == 0 && $plot == "multiplot"} {
      set plot [multiplot -xsize 600 -ysize 400 -title "Secondary Structure" -xlabel $xlabel -ylabel "Percentage (%)" -nolines -linewidth 2 -bkgcolor \#d1efff \
      -xmin [expr [lindex $reslistplot 0] - 0.5] -xmax [expr [lindex $reslistplot end] + 0.5] -ymin 0 -ymax 100]


      set ssList [list H T C G E B I]
      set colors [list "#eb82eb" "#469696" "#ffffff" "#1414ff" "#ffff64" "#b4b400" "#e11414"]
      set legendList [list "Alpha Helix" "Turn" "Coil" "3-10 Helix" "Beta Sheet" "Bridge" "Pi Helix"]
      set index 0
      foreach ssname $ssList color $colors {
        $plot add [expr [lindex $reslistplot 0] - 0.5] 0 -legend [lindex $legendList $index] -linecolor $color
        incr index
      }

      $plot replot
      for {set i 0} {$i < [llength $xplotdata]} {incr i} {
        set residue [lindex $xplotdata $i]
        set data [lindex $yplotdata $i]
        set xmin [expr [lindex $xplotdata $i] - 0.5]
        set xmax [expr [lindex $xplotdata $i] + 0.5]
        
        set min 0.0
        for {set j 0} {$j < [llength $data]} {incr j} {
          set val [expr [lindex $data $j] + $min]
          $plot draw rectangle $xmin $min $xmax $val -fill [lindex $colors $j] -tag "${residue}ss$j"
          set min $val
        }
      }
      $plot replot

      variable [$plot namespace]::c
      [$plot getpath].f.cf move legend -150 0
      update
      update idletasks
      $c postscript -file plot_${output}.ps 
      
      update
      update idletasks
    
      ### convert image from postscript to gif
      ### uses ImageMagick, only available on unix machines
      ### Windows users can use their own image editor to convert the ps file  
      global tcl_platform env

      if {$::tcl_platform(os) == "Darwin" || $::tcl_platform(os) == "Linux"} {
        eval ::ExecTool::exec "convert plot_${output}.ps gif:./plot_${output}.gif"
        update
        update idletasks
      }
      
      ### close the plot in case not intended 
      if {$showplot == 0} {
        $plot quit
      }
    } else {
      eval ::ExecTool::exec "gnuplot [::SSAnalysis::make_histogram_gnuplot ${xlabel} $xplotdata $yplotdata $output]"

      update idletasks
    }

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
    ::SSAnalysis::savetrajframes $mol $final_results_index ${output}_sorted_ssseq_traj pdbdcd

    puts $f "Closest structure file: ${output}_sorted_ssseq_highest_score.pdb"
    puts $f "Remaining structures sorted by score: ${output}_sorted_ssseq_traj.dcd"
    $mainsel delete
  }
  ## Search ss elements in atom selection pair. Search for
  ## consecutive ss elements in the atom selection. - e.g. {{"H" "resid 1 to 10","C" "resid 15 to 22"}} "
  if {$seqfind == 1} {
    set ssseqindex 0
    foreach pattern $seq {
      set sseq [split $pattern ","]
      set selList [list]
      set sselement [list]
      set sslistpat [list]
      foreach pair $sseq {
        lappend sselement [lindex $pair 0]
        lappend selList [atomselect $mol "[lindex $pair 1] and name CA"]
      }
      
      ### Generate a list containing all secondary structure for all selections,
      ### for all frames
      for {set i 0} {$i < $numframes} {incr i} {
          animate goto $i
          display update ui
          mol ssrecalc $mol 
          foreach sel $selList {
            $sel frame $i
            lappend sslistpat [lsort -unique [$sel get structure]]
          }

      }

      puts $f "\n\n### Pattern search\n"
      puts $f "Searching for the secondary structure element(s) in atom selection pairs(s): \"$sseq\" in all frames(N=$numframes)."
      
      ### Store the frames found to have all the patterns specified by the user 
      ### to save later in one dcd file
      array unset foundFrame
      array set foundFrame ""
      if {[info exists foundFrame($ssseqindex)] == 0} {
        set foundFrame($ssseqindex) [list]
      }
      for {set frame 0} {$frame < $numframes} {incr frame} {
        set found 0
        
        set incrframe [expr $frame * [llength sselement]]
        for {set i 0} {$i < [llength $sselement]} {incr i} {
          ### The index for each sselement on each frame in the sslistpat is
          ### $i + [expr $frame * [llength sselement] ]

          if {[lindex $sselement $i] == [lindex $sslistpat [expr $i + $incrframe] ]} {
            incr found
          }
          
        }
        ### Expecting the secondary element is one character element
        if {$found == [llength $sselement]} {
          lappend foundFrame($ssseqindex) $frame
        }
      }
      
      if {[llength $foundFrame($ssseqindex)] > 0} {
        puts $f "The secondary structure element(s) in atom selection pairs(s): $pattern was found in the frame(s) [join $foundFrame($ssseqindex)]."
        
        ### save trajectory with the searching pattern
        set filename [string map {"," "_"} [string map {\" ""} [string map {" " ""} pattern_index_[string trim ${pattern}]_frames]]]
        puts $f "Frames containing secondary structure element(s) in atom selection pairs(s): $pattern were saved in the pair pdb&dcd with the file name $filename"
        ::SSAnalysis::savetrajframes $mol $foundFrame($ssseqindex) $filename pdbdcd
      }
      incr ssseqindex
    }
    puts $f "\n\n"
    flush $f
  }
  close $f
  puts "Secondary Analysis finished."
  return
}

### Proc to save the frames of a trajectory in any order defined by the user
### Type variable receives two values pdbdcd - if the trajectory is be saved as pdb and dcd files
### or psfdcd - if the trajectory is be saved as pdb and dcd files
proc ::SSAnalysis::savetrajframes {mol frameslist output type} {
    set numframes [molinfo $mol get numframes]
    if {[llength $frameslist] > 1} {
      animate write dcd tmp_traj.dcd beg 0 end [expr $numframes -1] skip 1 $mol
    }
    set sel [atomselect $mol "all" frame [lindex $frameslist 0]]
    switch -exact -- $type {
      pdbdcd {
        $sel writepdb tmp_initialpdb.pdb
      }
      psfdcd {
        $sel writepdb tmp_initialpdb.psf
      }
      default {puts "Error saving trajectory"}
    }
    
    $sel delete
    if {[llength $frameslist] > 1} {
      set newmol [mol new tmp_initialpdb.pdb]
      set i 1
      if {$type == "psfdcd"} {
        set i 0
      }
      for {set j $i} {$j < [llength $frameslist]} {incr j} {
        mol addfile tmp_traj.dcd first [lindex $frameslist $j] last [lindex $frameslist $j] waitfor all
      }
      mol bondsrecalc $newmol
      mol ssrecalc $newmol
      animate write dcd ${output}.dcd beg $i end [expr [llength $frameslist] -1] skip 1 $newmol
      file delete -force tmp_traj.dcd
      
    }
    set ext "pdb"
    if {$type == "psfdcd"} {
      set ext "psf"
    }
    file rename -force tmp_initialpdb.$ext ${output}.$ext
    
    if {[info exists newmol]} {mol delete $newmol}
}

### Generate the postscript file from gnuplot based Maximilian Scheurer
### initial script
### xlabel - Residue ID or Residue Number 
### xlist - list of point in the x axis Residue ID
### ylist - list of the values of percentage for each residue in xlist
###         in the specific order: "{{H T C G E B I} {H T C G E B I} }"
proc ::SSAnalysis::make_histogram_gnuplot { xlabel xlist ylist output } {

  set file [open plot_${output}.gp w+]
  
  puts $file "set terminal postscript enhanced color solid"
  puts $file "set key top left outside horizontal autotitle columnhead"
  puts $file "set output \"plot_${output}.ps\""
  puts $file "set ytics rotate out font \",9\""
  set min [expr [lindex $xlist 0] -1]
  set max [expr [lindex $xlist end] +1]
  
  puts $file "set yrange \[ 0 : 100 \]"
  puts $file "set xlabel \"$xlabel\""
  puts $file "set title \"Secondary Structure\""
  puts $file "set ylabel \"Percentage(%)\""
  puts $file "set style data histograms"
  puts $file "set style histogram rowstacked"
  puts $file "set boxwidth 1"
  puts $file "set style fill solid 1.0 border -1"
  puts $file "set samples 1000000"

  puts $file "set autoscale xy"
  puts $file "\n\nplot \'-\' using 2:xtic(1) linecolor rgb \"magenta\", \'-\' using 2 linecolor rgb \"sea-green\", \'-\' using 2 linecolor rgb \"gray90\", \'-\' using 2 linecolor rgb \"blue\", \'-\' using 2 linecolor rgb \"yellow\", \'-\' using 2 linecolor rgb \"dark-yellow\", \'-\' using 2 linecolor rgb \"red\""

  set ssList [list H T C G E B I]
  
  for {set ssind 0 } {$ssind < [llength $ssList]} {incr ssind} {
    puts $file "\"${xlabel}\"\t[lindex $ssList $ssind]"
  
    for {set i 0} {$i < [llength $xlist]} {incr i} {
      puts  $file "[lindex $xlist $i]\t[lindex [lindex $ylist $i] $ssind]"
    }
    puts $file "e"
  }
  set maxx [expr round([expr [llength $xlist]*0.50])]
  set incrint [expr int([expr [llength $xlist]/$maxx])]
  set incrint [format %.0f [expr double(round(1*$incrint))/1]]
  set xtics ""
  set index 0
  for {set i [expr $min +1]} {$i < $max} {incr i $incrint} {
    append xtics "\"$i\" $index, "
    incr index $incrint
  }
  append xtics ""
  set xtics [string trimright $xtics ", "]
  puts $file "unset xtics"
  puts $file "set terminal postscript enhanced color solid"
  puts $file "set key top left outside horizontal autotitle columnhead"
  puts $file "set output \"plot_${output}.ps\""
  puts $file "set xtics font \",9\""
  puts $file "set xtics rotate out"
  puts $file "set xtics \($xtics\)"
  puts $file "set xrange \[ -1 \: [expr $index +1] \]"
  puts $file "refresh"
  puts $file "reset"
  puts $file "clear"
  puts $file "quit"
  close $file
  return "plot_${output}.gp"
}

