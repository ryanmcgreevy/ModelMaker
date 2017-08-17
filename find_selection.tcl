# find a selection residue range for rosetta in given PDB
# Author: Maximilian Scheurer, April 2016

namespace eval ::FindSelection {
	namespace export find_selection
	namespace export get_anchor_residue
	set version 0.1
	set description "Convert VMD selection to Rosetta selection"
}
package provide FindSelection $FindSelection::version

proc find_selection {args} \
{
	return [eval ::FindSelection::find_selection $args]
}

proc get_anchor_residue {args} \
{
	return [eval ::FindSelection::get_anchor_residue $args]
}

proc ::FindSelection::get_anchor_residue {MOL selection} \
{
	set searchmol [mol new $MOL.pdb]
	set searchsel [atomselect $searchmol "$selection"]
	set vmd_resids [lsort -unique -integer [$searchsel get residue]]
	set ros_resids []
	foreach vmd_resid $vmd_resids {
		lappend ros_resids [expr $vmd_resid +1]
	}
	mol delete $searchmol
	if {[llength $ros_resids] > 1} {
		puts "WARNING: Please specify only ONE anchor residue! Will use the first residue in the list only!"
	}
	if {[llength $ros_resids] == 0} {
		puts "ERROR: atom selection for anchor residue is empty! Exiting."
		exit
	}
	return [lindex $ros_resids 0]
}


#config should contain bb and chi settings for span for each selection, such as [list {0 1} {1 1}]
# first: chi, second: bb
proc ::FindSelection::find_selection {MOL selections config {offset 4} args} \
{
	puts "Trying to find selections."
	puts $selections
	set seltexts []
	foreach seltext $selections {
		lappend seltexts $seltext
	}
	set searchmol [mol new $MOL.pdb]

	set final_spans []
	set final_configurations []
	foreach findseltext $seltexts cf $config {
		set searchsel [atomselect $searchmol "$findseltext"]
		set vmd_resids [lsort -unique -integer [$searchsel get residue]]
		set ros_resids []
		foreach vmd_resid $vmd_resids {
			lappend ros_resids [expr $vmd_resid +1]
		}
		set counter -1
		set span_closed 1
		set current_span []
		set all_spans []
		for {set i 0} {$i < [llength $ros_resids]} {incr i} {
			if {$i == 0 && $counter == -1} {
				set counter [lindex $ros_resids $i]
				lappend current_span $counter
			} elseif {[lindex $ros_resids $i] == [expr $counter + 1]} {
				set counter [lindex $ros_resids $i]
			} else {
				lappend current_span $counter
				lappend all_spans $current_span
				set current_span []
				set counter [lindex $ros_resids $i]
				lappend current_span $counter
			}

			if {$i == [expr [llength $ros_resids] -1]} {
				lappend current_span $counter
				lappend all_spans $current_span
			}
		}
		set nSpans [llength $all_spans]
		for {set i 0} {$i < $nSpans} {incr i} {
			lappend final_configurations [list $cf]
		}
		lappend final_spans $all_spans
	}

	set out_ros [open "$MOL-spans.txt" w]
	set span_config []
	set single_spans []
	set out_spans ""
	foreach span_list $final_spans cf $final_configurations {
		foreach single_span $span_list {
			lappend single_spans $single_span
			puts $out_ros "{\"Span\" [lindex $single_span 0] [lindex $single_span 1] [lindex $cf 0] [lindex $cf 1]} " nonewline
			append out_spans "{\"Span\" [lindex $single_span 0] [lindex $single_span 1] [lindex $cf 0] [lindex $cf 1]} "
			lappend span_config "\"Span\" [lindex $single_span 0] [lindex $single_span 1] [lindex $cf 0] [lindex $cf 1]"
		}
	}
	#turn on for debugging
	puts $span_config

	set minimal_length [expr $offset * 2.0]
	#sort spans
	set sorted_spans [lsort -integer -index 0 $single_spans]
	set tot_resnum [llength [lsort -unique -integer [[atomselect $searchmol all] get residue]]]
	set tot_chainnum [llength [lsort -unique [[atomselect $searchmol all] get chain]]]
	set exclude_string ""
	set current_start 1
	set before_start 1
	set current_stop 1
	set before_end 1
	set span_open 0
	set constraint ""
	for {set i 1} {$i <= $tot_resnum} {incr i} {
		# puts $i
		foreach span $sorted_spans {
			if {$i == [lindex $span 0] && $i != 1} {
				set beg [lindex $span 0]
				set end [lindex $span 1]
				set spanlength [expr $end - $beg + 1]
				set middle [expr int(floor(double($spanlength)/2.0))]
				if {$spanlength <= $minimal_length} {
					# puts "span limit length!"
					set current_stop [expr $i - 1 + $middle] ;#[expr $i - 1 + $middle]
					append exclude_string "$current_start-$current_stop,"
					# puts "$current_start-$current_stop"
					set before_start $current_start
					set current_start [expr $end + 2 - $middle]
					set before_end $end
					set i $current_start ;# i is increased to $current_start + 1
				} else {
					set current_stop [expr $i - 1 + $offset]
					append exclude_string "$current_start-$current_stop,"
					# puts "$current_start-$current_stop"
					set before_start $current_start
					set current_start [expr $end + 1 - $offset]
					set before_end $end
					set i $current_start ;# i is increased to $current_start + 1
				}
				# puts "beginning of span found: $span"
				# puts "span length: $spanlength"
			} elseif {$i == 1 && $i == [lindex $span 0]} {
				set end [lindex $span 1]
				if {[expr $end - $offset] >= 1} {
					set end_f [expr $end - $offset]
				} else {
					set end_f $end
				}
				set current_start [expr $end_f + 1]
			}
		}
		if {$i == $tot_resnum && $i != $before_end} {
			if {$current_start > $i} {
				append exclude_string "1"
			} else {
				append exclude_string "$current_start-$i"
			}

			# puts "$current_start-$i"
		}
	}

	set current_start 1
	set before_start 1
	set current_stop 1
	set before_end 1
	set span_open 0
	set constraint ""

	for {set i 1} {$i <= $tot_resnum} {incr i} {
		# puts $i
		foreach span $sorted_spans {
			if {$i == [lindex $span 0] && $i != 1} {
				set beg [lindex $span 0]
				set end [lindex $span 1]
				set spanlength [expr $end - $beg + 1]
				set current_stop [expr $i - 1]
				append constraint "{$current_start $current_stop} "
				# puts "$current_start-$current_stop"
				set before_start $current_start
				set current_start [expr $end + 1]
				set before_end $end
				set i $current_start ;# i is increased to $current_start + 1
				# puts "beginning of span found: $span"
				# puts "span length: $spanlength"
			} elseif {$i == 1 && $i == [lindex $span 0]} {
				set end [lindex $span 1]
				set current_start [expr $end + 1]
			}
		}
		if {$i == $tot_resnum && $i != $before_end} {
			append constraint "{$current_start $i}"
			# puts "$current_start-$i"
		}
	}


	# build fix chains!!!
	set chain_config ""
	for {set i 1} {$i <= $tot_chainnum} {incr i} {
		append chain_config "{\"Chain\" $i 0 0} "
	}

	puts $out_ros $exclude_string
	mol delete $searchmol
	close $out_ros
	puts "finished finding selections."
	set lastchar [string range $exclude_string end end]
	if {$lastchar == ","} {
		return [list $out_spans [string range $exclude_string 0 [expr [string length $exclude_string]-2]] $chain_config $constraint]
	}
	return [list $out_spans $exclude_string $chain_config $constraint]
	#return [list $out_spans [string range $exclude_string 0 [expr [string length $exclude_string]-2]] $chain_config $constraint]
#  return [list $out_spans [string range $exclude_string 0 [expr $end-1]] $chain_config $constraint]
}
