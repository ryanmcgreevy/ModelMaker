
namespace eval ::SSAnalysis {
  set version 0.1
  set description "ss analysis"
}
package provide SSAnalysis $SSAnalysis::version

proc ss_analysis {args} \
{
  return [eval ::SSAnalysis::ss_analysis $args]
}

proc ::SSAnalysis::ss_analysis {MOL max_structures resid_start resid_stop chain_id} \
{
  mol delete all
 set input_structure {}
 set ss_list {}
 set residrange {1 5}
 # open a file to write the output to
 set f [open "ss_${max_structures}_${resid_start}_${resid_stop}.txt" "w+"]

 # Splitting the secondary structure pattern provided by the user into individual list arguments.
 set avg_new [split $input_structure {}]

 # loading the predicted structure to VMD and creating the defined atom selections
 mol new ${MOL}_rosetta_scoring_min_$max_structures.pdb waitfor all
 set input_num [molinfo top get numframes]
 set sel [atomselect top "chain $chain_id and resid $resid_start to $resid_stop and name CA"]
 set sel2 [atomselect top "all"]

 # puts "number of selection: [$sel num] - input num: $input_num - all: [$sel2 num]"
 # definition of the list of secondary structure motives
 set list_H "" ;# Alpha Helix
 set list_T "" ;# Turn
 set list_C "" ;# Coil
 set list_E "" ;# Extended Comformation
 set list_G "" ;# 3-10 Helix
 set list_B "" ;# Bridge
 set list_I "" ;# Pi Helix
 set ssseq ""
 set residues [lsort -unique -real [$sel get resid]]

 for {set i 0} {$i < [llength $residues]} {incr i} {
  lappend list_H 0
  lappend list_T 0
  lappend list_C 0
  lappend list_E 0
  lappend list_G 0
  lappend list_B 0
  lappend list_I 0
  }

  # evaluating the occurancy of each secondary structure element of each frame
  for {set i 0} {$i < $input_num} {incr i} {
    # puts "frame: $i"
    animate goto $i
    display update ui
    $sel frame $i
    mol ssrecalc 0

    set j 0
    lappend ssseq [$sel get structure]
    foreach val [lindex $ssseq $i] {
      switch $val {   
        "H" {lset list_H $j [expr [lindex $list_H $j] +1 ]}
        "T" {lset list_T $j [expr [lindex $list_T $j] +1 ]}
        "C" {lset list_C $j [expr [lindex $list_C $j] +1 ]}
        "E" {lset list_E $j [expr [lindex $list_E $j] +1 ]}
        "G" {lset list_G $j [expr [lindex $list_G $j] +1 ]}
        "B" {lset list_B $j [expr [lindex $list_B $j] +1 ]}
        "I" {lset list_I $j [expr [lindex $list_I $j] +1 ]}
      }
      incr j
    }

  }

  #Iteration over all input structures to compute the occurrence of each secondary structure element (in percentage) for each residue

  set avg_seq ""
  for {set i 0} {$i < [llength $residues]} {incr i} {
    set percent_H [expr [lindex $list_H $i]*1.0/$input_num*100.0]
    set percent_T [expr [lindex $list_T $i]*1.0/$input_num*100.0]
    set percent_C [expr [lindex $list_C $i]*1.0/$input_num*100.0]
    set percent_G [expr [lindex $list_G $i]*1.0/$input_num*100.0]
    set percent_E [expr [lindex $list_E $i]*1.0/$input_num*100.0]
    set percent_B [expr [lindex $list_B $i]*1.0/$input_num*100.0]
    set percent_I [expr [lindex $list_I $i]*1.0/$input_num*100.0]
    set val [lsort -indices -decreasing -real [list $percent_H $percent_T $percent_C $percent_E $percent_E $percent_B $percent_I]]
    switch [lindex $val 0] {
      0 {lappend avg_seq H}
      1 {lappend avg_seq T}
      2 {lappend avg_seq C}
      3 {lappend avg_seq E} 
      4 {lappend avg_seq G}
      5 {lappend avg_seq B}
      6 {lappend avg_seq I}

    }
    puts $f "Resid [lindex $residues $i]:\tH: [format %0.2f $percent_H] \tT: [format %0.2f $percent_T] \t C: [format %0.2f $percent_C] \t G: [format %0.2f $percent_G] \t E: [format %0.2f $percent_E] \t B: [format %0.2f $percent_B] \t I: [format %0.2f $percent_I]"
  }

  puts $f "Average sequence :  $avg_seq"


  if {[llength $avg_new] != 0} {
    puts $f "User input provided: $input_structure"
    set avg_seq $avg_new
  }
  set score_list ""

  for {set j 0} {$j < $input_num} {incr j} {
    set score_counter 0
    set frm_seq [lindex $ssseq $j]
    for {set k 0} {$k< [llength $frm_seq]} {incr k} {
      if {[lindex $frm_seq $k] == [lindex $avg_seq $k] } {incr score_counter} 
    }
    lappend score_list $score_counter
  }

  # Identifying the representative structures meaning the ones representing the secondary structure sequence closest to the average secondary structure sequence above. 


  set final_results_index [lsort -real -indices -decreasing $score_list]
  set final_results [lsort -decreasing -real $score_list]
  set range [lsearch -all $final_results [lindex $final_results 0]]
  set frames [lrange $final_results_index 0 [lindex $range end]]
  puts $f "final results (frames with highest scores - [lindex $final_results 0]) : $frames"
  set hit_number [llength $frames]
  puts $f "number of hits: $hit_number"
  
  mol new ${MOL}_rosetta_scoring_min_${max_structures}.pdb first [lindex [lindex $frames] 0] last [lindex [lindex $frames] 0] -waitfor all
  
  for {set j 1} {$j < [llength [lindex $frames]]} {incr j} {
    mol addfile ${MOL}_rosetta_scoring_min_${max_structures}.pdb first [lindex [lindex $frames $j]] last [lindex [lindex $frames $j]] -waitfor all
  }
  animate write pdb ss_average_${max_structures}.pdb beg 0 end [expr $hit_number-1] skip 1 top
  
  set score_check [lindex $final_results 0]
  set frame_cluster {}
  set num_cluster 0
  
  # Search for consecutive alpha-helix pattern in the range defined above by the variable "residrange"
  # set searchindex ""
  # set frames ""
  # for {set i 0} {$i< [llength $ssseq]} {incr i} {
  #   for {set j 0} {$j < [llength $residrange]} {incr j} {
  #     set index [lindex $residrange $j]
  #     set ini [expr [lindex $index 0] - $resid_start] 
  #     set end [expr [lindex $index 1] - $resid_start] 
  #     set pattern [string repeat "H" [expr  $end - $ini +1] ] 
  #     #puts [string trim [lindex $ssseq $i] " "]
  #     set strng [join [lrange [lindex $ssseq $i] $ini $end]  ""]
  #     #puts "pattern $pattern string $strng ssseq [lindex $ssseq $i]"
      
  #     if {$strng == $pattern} {
  #       append searchindex "frame :[expr $i +1] resid range: $ini - $end\n"
  #       lappend frames $i 
  #     }
  #   }
    
  # }
  # puts $f "Pattern search result : $searchindex"
  
  # #Writes all structures containing the consecutive alpha-helix pattern in the range defined above by the variable "residrange" to a pdb file
  # mol new ${MOL}_rosetta_scoring_min_${max_structures}.pdb first [lindex [lindex $frames] 0] last [lindex [lindex $frames] 0] -waitfor all
  
  # for {set j 1} {$j < [llength [lindex $frames]]} {incr j} {
  #   mol addfile ${MOL}_rosetta_scoring_min_${max_structures}.pdb first [lindex [lindex $frames $j]] last [lindex [lindex $frames $j]] -waitfor all
  # }
  # animate write pdb ss_pattern_${max_structures}.pdb beg 0 end [expr [llength $frames] -1] skip 1 top
  
  close $f
}