### Utilities for controller_linux.tcl
### Maximilian Scheurer, April 2016

namespace eval ::RosettaUtilities {
	set description "Utilities for RosettaVMD"
	set version 0.1
}
package provide RosettaUtilities $RosettaUtilities::version

proc align_rosetta_local {args} \
{
	return [eval ::RosettaUtilities::align_rosetta_local $args]
}

proc align_rosetta_cluster {args} \
{
	return [eval ::RosettaUtilities::align_rosetta_cluster $args]
}

proc make_cluster_input {args} \
{
	return [eval ::RosettaUtilities::make_cluster_input $args]
}

proc run_clustering {args} \
{
	return [eval ::RosettaUtilities::run_clustering $args]
}

proc run_vmd_clustering {args} \
{
	return [eval ::RosettaUtilities::run_vmd_clustering $args]
}

proc run_rosetta_clustering {args} \
{
	return [eval ::RosettaUtilities::run_rosetta_clustering $args]
}

proc read_cluster_file  {args} \
{
	return [eval ::RosettaUtilities::read_cluster_file $args]
}

proc evaluate_ss_analysis {args} \
{
	return [eval ::RosettaUtilities::evaluate_ss_analysis $args]
}

proc get_unass_dens {args} \
{
	return [eval ::RosettaUtilities::get_unass_dens $args]
}

proc getSS {args} \
{
	return [eval ::RosettaUtilities::getSS $args]
}

proc assignSS {args} \
{
	return [eval ::RosettaUtilities::assignSS $args]
}

proc ::RosettaUtilities::align_rosetta_local {start end MOL template tempMol sel_temp sel_rosetta} \
{
	
  set ml [mol new $template]

	set temp [atomselect $ml "$sel_temp"]

	for {set x $start} {$x<=$end} {incr x} {
    if {$end < 10000} {
      set offset 4
    } else {
      set offset [expr int(floor(log10($end))) + 1] 
    }
    if { [file exists ./pdb_out_aligned/${MOL}_${x}_0001_aligned.pdb] != 1 } {  
      set index [mol new ./pdb_out/${MOL}_[format %0${offset}i  $x].pdb]
      set sel [atomselect $index "$sel_rosetta"]
      set transformation_matrix [measure fit $sel $temp]
      set move_sel [atomselect $index "all"]
      $move_sel move $transformation_matrix
      $move_sel writepdb ./pdb_out_aligned/${MOL}_${x}_0001_aligned.pdb
      $sel delete
      $move_sel delete
      mol delete $index
    }
	}
}

proc ::RosettaUtilities::align_rosetta_cluster {start end MOL template tempMol sel_temp sel_rosetta} \
{
	set tempmol [mol new $template]
	set temp [atomselect $tempmol "($sel_temp) and backbone and noh"]

	for {set x $start} {$x<=$end} {incr x} {
		if { [file exists "./pdb_out/${MOL}${x}_0001.pdb"] == 1} {
		set index [mol new ./pdb_out/${MOL}${x}_0001.pdb]
		set sel [atomselect $index "($sel_rosetta) and backbone and noh"]
		set transformation_matrix [measure fit $sel $temp]
		set move_sel [atomselect $index "all"]
		$move_sel move $transformation_matrix
		$move_sel writepdb ./pdb_out_aligned/${MOL}_${x}_0001_aligned.pdb
		$sel delete
		$move_sel delete
		mol delete $index
		}
	}
}

proc ::RosettaUtilities::make_cluster_input {MOL start end max_structures} \
{
	#set start 267
	#set end 286
	#set MOL ubp6_5a5b_8_g4S
	#set max_structures 2500

	mol new ${MOL}_rosetta_scoring_min_${max_structures}.pdb waitfor all

	set all [atomselect top "resid $start to $end" frame 0]
	#Move to center
	# for {set i 0} {$i < [molinfo top get numframes]} {incr i} {
	# 	set all [atomselect top "$selection" frame $i]
		#set sel [atomselect top "backbone and resid $start to $end" frame $i]
		#$all moveby [vecinvert [measure center $sel]]
	# }

	#align
	#for {set i 0} {$i < [molinfo top get numframes]} {incr i} {
		#        set all [atomselect top "resid $start to $end" frame $i]
		#		set sel [atomselect top "backbone and resid 162 to 172 or resid 139 to 146" frame $i]
		#        set sel1 [atomselect top "backbone and resid 162 to 172 or resid 139 to 146" frame [expr $i-1]]
		#        set transformation_matrix [measure fit $sel $sel1]
		#        $all move $transformation_matrix
		#}

		animate write dcd cluster_input_${start}_${end}_${max_structures}.dcd beg 0 end [expr $max_structures-1] skip 1 sel $all top
		animate write pdb cluster_input_${start}_${end}_${max_structures}.pdb beg 0 end [expr $max_structures-1] skip 1 sel $all top
		animate write pdb cluster_input_start_${start}_${end}_${max_structures}.pdb beg 0 end 0 skip 1 sel $all top


}

proc ::RosettaUtilities::run_clustering {mol start end bestN kmin kmax} \
{
  puts "Running python clustering..."
  set name cluster_input_${start}_${end}_${bestN}

  set mol [mol new ${name}.pdb]
  set sel [atomselect $mol "all"]
  #$sel set occupancy 1
  #$sel set beta 1
  #$sel writepdb ${name}.pdb
  $sel writepsf ${name}.psf

#  set mol [mol new ${name}.pdb]
#  set sel [atomselect $mol "all"]
#  set chains [$sel get chain]
#  set segnames [$sel get segname]
##get unique list
#  foreach chain $chains {dict set tmp $chain 1}
#  set chains [dict keys $tmp]
#  foreach segname $segnames {dict set tmp2 $segname 1}
#  set segnames [dict keys $tmp2]

#  foreach chain $chains segname $segnames {
#    lappend chseg  [list ${name}.pdb $chain $segname]
#  }
  
  global env 
  #set topfiles [list [file join $env(CHARMMTOPDIR) top_all36_prot.rtf] \
    [file join $env(CHARMMTOPDIR) top_all36_lipid.rtf] \
    [file join $env(CHARMMTOPDIR) top_all36_na.rtf] \
    [file join $env(CHARMMTOPDIR) top_all36_carb.rtf] \
    [file join $env(CHARMMTOPDIR) top_all36_cgenff.rtf] \
    [file join $env(CHARMMTOPDIR) toppar_all36_carb_glycopeptide.str] \
    [file join $env(CHARMMTOPDIR) toppar_water_ions_namd.str] ]

  #auto_makepsf ${name}.pdb $topfiles "" ""

  exec python $env(RosettaVMDDIR)/clusterK.py -psf ${name}.psf -dcd ${name}.dcd -kmin $kmin -kmax $kmax > clustering.log
  
  #create dcds for each cluster
  set fp [open "clu_labels.txt" r]
  set file_data [read $fp]
  close $fp
  set data [split $file_data "\n"]
  foreach line $data {
    if {$line != ""} {
      lappend ar([lindex $line 0]) [lindex $line 1]
    }
  }
  foreach cluster [array names ar] {
    mol new ${name}.psf
    foreach framenum $ar($cluster) {
      animate read dcd ${name}.dcd beg $framenum end $framenum
    }
    animate write dcd "cluster-${cluster}.dcd"
  }

  #create a pdb for each medoid
  set fp [open "clu_medoids.txt" r]
  set file_data [read $fp]
  close $fp
  set data [split $file_data "\n"]
  foreach line $data {
    if {$line != ""} {
      lappend arm([lindex $line 0]) [lindex $line 1]
    }
  }
  foreach cluster [array names arm] {
    mol new ${name}.psf
    set framenum $arm($cluster)
    animate read dcd ${name}.dcd beg $framenum end $framenum
    set sel [atomselect top all]
    $sel writepdb "cluster-${cluster}-medoid.pdb"
    $sel delete
  }

}

proc ::RosettaUtilities::run_vmd_clustering {mol start end bestN cutoff cluster_number} \
{
	puts "VMD clustering"
	#VMD rmsd based clustering
	package require math::statistics
	set clust_num $cluster_number ;# number of clusters
	# set cutoff 0.25 ;# RMSD within each cluster

	# open a file to write the output text to
	set log [open cluster_out.txt w+]

	# loading the predicted structures sorted by the rosetta energy scoring script to VMD and creating the defined atom selections
	# mol new ${MOL}_rosetta_scoring_min_${max_structures}.pdb waitfor all
	set cl_mol [mol new cluster_input_${start}_${end}_${bestN}.pdb waitfor all]
	set clustList [measure cluster [atomselect $cl_mol "protein and backbone"] num $clust_num distfunc rmsd cutoff $cutoff first 0 last -1 step 1]
	puts $log $clustList
	# performing a RMSD based cluster analysis
	for {set j 0} {$j < [expr [llength $clustList] -1]} {incr j} {
		set a [lindex [lindex $clustList $j] 0]
		set a1 [expr $a+1]
		set index [mol new cluster_input_${start}_${end}_${bestN}.pdb first $a last $a waitfor all]
		set l  [llength [lindex $clustList $j] ]
		puts $log "------------------------------------------"
		puts $log "Cluster $j $l"
		set l1 [expr $l-1]
		set clustList_aux [lindex $clustList $j]
		for {set i 1} {$i < $l} {incr i} {
			set x [lindex $clustList_aux $i]
			set x1 [expr $x+1]
			puts "x $x1"
			puts "l $l"
			mol addfile cluster_input_${start}_${end}_${bestN}.pdb first $x last $x waitfor all
		}
  		# write the structures of each cluster to a seperate pdb file
  		animate write pdb cluster_$j.pdb beg 0 end $l1 skip 1 $index
  		animate dup frame $l1 $index
  		# calculate the representative structure of each cluster and write it to a pdb file
  		set sel [atomselect top "all"]
  		set av_struct [atomselect top all frame $l]
  		set av [measure avpos $sel]
  		set a [$av_struct set {x y z} $av]
  		set numframe [molinfo top get numframes]
  		set rmsd ""
  		for {set frame 0} {$frame < [expr $numframe -1]} {incr frame} {
  			$sel frame $frame
  			lappend rmsd [format %.2f [measure rmsd $sel $av_struct]]
  		}
		$av_struct delete

		set min [math::statistics::min $rmsd]
		set minframe [lsearch $rmsd $min]
		mol delete all
		mol new cluster_input_${start}_${end}_${bestN}.pdb waitfor all
		animate write pdb cluster_rep_${j}.pdb beg [lindex $clustList_aux $minframe] end [lindex $clustList_aux $minframe]

		puts $log "Cluster $j {$clustList_aux}"
		puts $log "Cluster $j representative [lindex $clustList_aux $minframe]"
		puts $log "Cluster $j {$min}"
		puts $log "------------------------------------------"
		puts $log "animate write pdb cluster_${j}_rep.pdb beg [lindex $clustList_aux $minframe] end [lindex $clustList_aux $minframe]"
	}
	close $log
}

proc run_rosetta_clustering {mol start end bestN radius cluster_max} \
{
	puts "running Rosetta clustering"
	global rosettapath
	global rosettaDBpath
	global platform
	set cl_mol [mol new cluster_input_${start}_${end}_${bestN}.pdb waitfor all]
	set f [molinfo $cl_mol get numframes]
	set pdblist [open "pdblist" "w"]
	set blacklist {}
	for {set i 0} {$i < $f} {incr i} {
		animate write pdb "cluster_single_${start}_${end}_${bestN}_${i}.pdb" beg $i end $i
		puts $pdblist "cluster_single_${start}_${end}_${bestN}_${i}.pdb"
		lappend blacklist "cluster_single_${start}_${end}_${bestN}_${i}.pdb"
	}
	close $pdblist
	exec [glob $rosettapath/cluster.$platform] -database $rosettaDBpath -in:file:l pdblist -in:file:fullatom -cluster:radius $radius -cluster:limit_clusters $cluster_max > cluster.log
	foreach item $blacklist {
		exec rm $item
	}
	exec rm pdblist
}

proc ::RosettaUtilities::read_cluster_file {MOL max_structures max_cluster} \
{
	set file [open "cluster_list_$max_structures.txt" r]
	set outfile [open "cluster_output.txt" w+]
	set info [split [read -nonewline $file] "\n"]

	set clusters ""
	for {set i 1} {$i <= $max_cluster} {incr i} {
		lappend clusters [lsearch -all -exact $info $i]
		puts "cluster $i [lindex $clusters [expr $i -1]]"
	}


	puts $outfile "Cluster\tNum_Frames\tFrames_Index"
	for {set i 0} {$i < $max_cluster} {incr i} {
		# mol new ${MOL}_rosetta_scoring_min_${max_structures}.pdb first [lindex [lindex $clusters $i] 0] last [lindex [lindex $clusters $i] 0] -waitfor all
		# for {set j 1} {$j < [llength [lindex $clusters $i] ]} {incr j} {
		# 	mol addfile ${MOL}_rosetta_scoring_min_${max_structures}.pdb first [lindex [lindex $clusters $i] $j] last [lindex [lindex $clusters $i] $j] -waitfor all
		# }
		# animate write pdb cluster_[expr $i +1].pdb beg 0 end [expr [llength [lindex $clusters $i] ] -1] skip 1 top

		puts $outfile "[expr $i +1]\t[llength [lindex $clusters $i] ]\t[lindex $clusters $i]"
		mol delete top
	}

	close $file


	set file [open "cluster_${max_structures}_out2.log" r]

	set info [split [read -nonewline $file] "\n"]

	set line [lsearch -all -exact $info {[1] "Most representative structure from each cluster  (trajectory index) :"}]

	set representative [lindex $info  [expr $line + 1]]
	set represent_cluster ""

	#mol new ${MOL}_rosetta_scoring_min_${max_structures}.pdb
	for {set i 1} {$i <= $max_cluster} {incr i} {
		puts $outfile "Representative struct from cluster $i [lindex $representative $i]"
		mol new ${MOL}_rosetta_scoring_min_${max_structures}.pdb first [lindex $representative $i] last [lindex $representative $i] -waitfor all
		animate write pdb cluster_rep_$i.pdb beg 0 end 0 skip 1 top
		mol delete top

	}

	close $outfile
}

proc ::RosettaUtilities::evaluate_ss_analysis {max_structures resid_start resid_stop} \
{
	global gnuplotexe
	set out [open "histogram_${max_structures}_${resid_start}_${resid_stop}.gp" w]
	set diff [expr double($resid_stop - $resid_start)]
	puts $out [::RosettaUtilities::make_histogram_gnuplot $max_structures $resid_start $resid_stop $diff]
	close $out
	exec $gnuplotexe "histogram_${max_structures}_${resid_start}_${resid_stop}.gp"
}

proc ::RosettaUtilities::make_histogram_gnuplot {max_structures resid_start resid_stop diff} \
{
	return "
	reset
	set term postscript enhanced color
	set output \"plot_${max_structures}_${resid_start}_${resid_stop}.ps\"
	set xtics font \",9\"
	set xtics rotate out
	set xrange \[ -1 \: $diff \]
	set yrange \[ 0 : 100 \]

	set xlabel \"Residue ID\"
	set ylabel \"Secondary Structure\"

	set noytics
	# set noxtics

	#set grid y
	set style data histograms
	set style histogram rowstacked
	set boxwidth 1
	set style fill solid 1.0 border -1


	plot   \'ss_${max_structures}_${resid_start}_${resid_stop}.txt\' using 4:xticlabels(2) t \"H\" linecolor rgb \"magenta\", \'\' using 10 t \"G\" linecolor rgb \"blue\", \'\' using 6 t \"T\" linecolor rgb \"sea-green\", \'\' using 8 t \"C\" linecolor rgb \"gray90\", \'\' using 12 t \"E\" linecolor rgb \"yellow\", \'\' using 14 t \"B\" linecolor rgb \"dark-yellow\", \'\' using 16 t \"I\" linecolor rgb \"red\"

	quit
	"
}

# Get unassigned density
# Till Rudack, March 2016
#

#volutile only works with dx densities do griddx twice on the ccp4 to get a dx density
# mdff griddx -i ccp4 -o potential.dx
# mdff griddx -i potential.dx -o density.dx

proc ::RosettaUtilities::get_unass_dens {MOL selection mapname cutoff outname} \
{
	set comp [mol new $MOL.pdb]

	set sel [atomselect $comp "$selection"]

	volmap mask $sel -o mask.dx -cutoff $cutoff
	mdff griddx -i mask.dx -o mask_invert.dx
	volutil -mult ${mapname}.dx mask_invert.dx -o ${outname}_${cutoff}_${mapname}_unassigned_density.dx
	exec rm -f mask.dx
	exec rm -f mask_invert.dx
#mdff griddx -i ${MOL}_${cutoff}_${DENS_NAME}_density.dx -o ${MOL}_${cutoff}_${DENS_NAME}_potential.dx

#volutil ${MOL}_${cutoff}_${DENS_NAME}_density.dx -o ${MOL}_${cutoff}_${DENS_NAME}_density.situs
#volutil ${MOL}_${cutoff}_${DENS_NAME}_potential.dx -o ${MOL}_${cutoff}_${DENS_NAME}_potential.situs

#situs
#map2map ${MOL}_${cutoff}_density.situs ${MOL}_${cutoff}_density.mrc
# option 1
}


#set fileout "ssfile.ss"
#set filein "ssfile.ss"
#set molid 0
#set breakby "chain"
#set tmpdir "./temp"


# Just run:
#    getSS $fileout $molid $breakby $tmpdir
#    assignSS $filein $molid $breakby $tmpdir

# Don't edit below this line
proc ::RosettaUtilities::getSS {fileout molid breakby tmpdir} {

    set fh [open $fileout "w"]

    set breaks [lsort -unique [[atomselect $molid "protein"] get $breakby]]
    set N [llength $breaks]
    puts "Detected $N breaks"

    foreach brk $breaks {
	set sel [atomselect $molid "protein and ${breakby} $brk"]
	set tmpfile "${tmpdir}/ss${breakby}${brk}.pdb"
	$sel writepdb $tmpfile
	set cmol [mol new $tmpfile]
	set ss [[atomselect $cmol "protein"] get structure]
	puts $fh $ss
	$sel delete
	mol delete $cmol
    }
    close $fh
}


proc ::RosettaUtilities::assignSS {filein molid breakby} {
    set breaks [lsort -unique [[atomselect $molid "protein"] get $breakby]]
    set fh [open $filein "r"]
    set file_data [read $fh]
    close $fh
    set data [split $file_data "\n"]
    set cdata [lrange $data 0 [expr [llength $data] - 2 ]]
    foreach line $cdata brk $breaks {
	# do some line processing here
	set sel [atomselect $molid "protein and ${breakby} ${brk}"]
	$sel set structure $line
	$sel delete
    }
}
