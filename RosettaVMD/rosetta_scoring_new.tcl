#package require Tcl 8.5
proc rosetta_scoring {length_tot score_name} {
	set file_list ""
	set val_list ""
	set id_list ""
	set ic_val_list ""
	set file_list [lsort -dictionary [glob ./rosetta_output/sc_out/*.sc]]
	#set file_list [lsort -dictionary [glob ../rosetta_output_bbfix/sc_out/*.sc]]
	#set file_list [lsort -dictionary [glob ../rosetta_output_bbfree/sc_out/*.sc]]

	set file_id_list [] ;# for total score
	set file_id_if_list [] ;# for interface score

	foreach file $file_list {
		#puts $file
		regexp -nocase [subst -nocommands -nobackslashes {${score_name}[0-9]*}] $file a 
		# puts $a
		set ident [string trim $a $score_name]
		# puts $ident

		set f [open "$file" r]
		set rfile [read $f]
		set line [lindex [split $rfile "\n"] 2]
		if {$line == ""} {
			set line [lindex [split $rfile "\n"] 0]

		}
		set val [lindex $line 1]
		set icval [lindex $line 4]
		close $f
		lappend file_id_list [list $ident $file $val]
		lappend file_id_if_list [list $ident $file $icval]
	}

	set sort_scores [lsort -real -index 2 $file_id_list]
	set sort_if_scores [lsort -real -index 2 $file_id_if_list]

	puts "finished reading and sorting."

	set length [llength $file_list]
	set log [open "rosetta_scoring.log" w+]
	set ic_log [open "rosetta_if_sc.log" w+]

	set counter 0
	foreach item $sort_scores {
		if {$counter < $length_tot} {
			lappend number [lindex $item 0]
		}
		puts $log "[lindex $item 0]\t[lindex $item 1]\t[lindex $item 2]"
		incr counter
	}

	set counter 0
	foreach item $sort_if_scores {
		if {$counter < $length_tot} {
			lappend ic_number [lindex $item 0]
		}
		puts $ic_log "[lindex $item 0]\t[lindex $item 1]\t[lindex $item 2]"
		incr counter
	}


	if {$length > $length_tot} {
		return [list [lrange $number 0 [expr $length_tot -1]] [lrange $ic_number 0 [expr $length_tot -1]]]
		} else {
			puts " Not enough file for chosen length_tot total file = $length: length_tot= $length_tot "
		}

	} 
