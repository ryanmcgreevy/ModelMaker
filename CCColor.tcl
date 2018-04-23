# Color pdb via CCC of Secondary Structure and residues
# Maximilian Scheurer, March 2016
#

namespace eval ::CCColor {
	namespace export cccolor
	set version 0.1
	set description "CCColoring tool as VMD plugin"
}
package provide CCColor $CCColor::version

proc cccolor {args} \
{
	return [eval ::CCColor::cccolor $args]
}
#####################
# input files
####################
proc ::CCColor::cccolor {MOL mapname resolution threshold spacing ss_on res_on gen_sel res_sel {tempfolder "/usr/tmp"} {bb_only 0} args} {
	#################
	# Program settings
	#################
	####color ss
	#set ss_on 1
	####color residues
	#set res_on 0
	# if "all" is specified as range, all residues will be analyzed, otherwise specify "resid x to y"
	set range $res_sel

	#######################
	## MAIN PROGRAM
	########################
	set pdb [mol new $MOL]

  set MOL [file rootname [file tail $MOL]]
	
  ### set beta for all
	set defaultBeta 0
	[atomselect $pdb all] set beta $defaultBeta

	set atomSel [atomselect $pdb "$gen_sel"]

	set dens_mol [mol new $mapname]

	#source getSS.tcl
	
	getSS $MOL.ss $pdb segname $tempfolder
	assignSS $MOL.ss $pdb segname

	set chains [lindex [list [lsort -unique [$atomSel get fragment]]] 0]
	set log [open "${MOL}_cclog.txt" w]
	puts $log $chains
	# do ccc-coloring for SS
	if {$ss_on} {
		puts $log "Running SS coloring"
		foreach chain $chains {
			puts $log "current chain: $chain"
			set chainCAS [atomselect $pdb "fragment $chain and name CA"]
			puts $chainCAS
			set ssList [$chainCAS get structure]
			puts $log $ssList
			puts $log [llength $ssList]
			set resList [lsort -unique -integer [$chainCAS get resid]]
			set resListLength [llength $resList]
			puts $log "reslist $resListLength"	
			set id_begin [lindex $resList 0]
			set id_end [expr $id_begin + $resListLength - 1]
			set id_0 $id_begin
			for {set id_i $id_begin} {$id_i <= $id_end} {incr id_i} {
				set ss_0 [lindex $ssList [expr $id_i - $id_begin]]
				set ss_1 [lindex $ssList [expr $id_i - $id_begin + 1]]
				#	puts $log "$id_i"
				if {$ss_0 != $ss_1} {
          set ssSel [atomselect $pdb "fragment $chain and (resid $id_0 to $id_i)"]
          if {[$ssSel num] == 0} {
						puts "no atoms in selection"
					}
					puts $log [$ssSel num]
					#    volmap mask $ssSel -o mask.dx -cutoff $cutoff
					#   volutil -mult ${emdbName}_density.dx mask.dx -o compare.dx
					#set CCSS [mdff ccc $ssSel -i $mapname.dx -res $resolution -spacing $spacing -thresholddensity $threshold]
					if {$spacing != -1} {
           set CCSS [mdffi cc $ssSel -mol $dens_mol -res $resolution -spacing $spacing -thresholddensity $threshold]
					} else {
           set CCSS [mdffi cc $ssSel -mol $dens_mol -res $resolution -thresholddensity $threshold]
          }
          #set CCSS 0.0
					set str [string range "$CCSS" 1 3]
					if {[string equal $str "nan"] || [string equal $str "NaN"]} {
						$ssSel set beta -1.0
						} else {
							$ssSel set beta $CCSS
						}
						#rm mask.dx compare.dx
						puts $log "$id_0 $id_i $CCSS"
						set id_0 [expr $id_i +1]
					}
        }
			}
      [atomselect $pdb all] writepdb $MOL-ss.pdb
		}

		if {$res_on} {
			puts $log "Running resid coloring"
			[atomselect $pdb all] set beta $defaultBeta
			foreach chain $chains {
				if {$range=="all"} {
					set resids [lsort -unique -integer [[atomselect $pdb "fragment $chain"] get resid]]
					} else {
						set resids [lsort -unique -integer [[atomselect $pdb "fragment $chain and (${range})"] get resid]]
					}
					puts $log "$resids"
					foreach res $resids {
						if {$bb_only} {
							set resSel [atomselect $pdb "fragment $chain and resid $res and backbone"]
						} else {
							set resSel [atomselect $pdb "fragment $chain and resid $res"]
						}
						puts $log "[$resSel num] --- $chain"

						#volmap mask $resSel -o mask.dx -cutoff $cutoff
						#volutil -mult ${emdbName}_density.dx mask.dx -o compare.dx
						#set CCres [mdff ccc $resSel -i $mapname.dx -res $resolution -spacing $spacing -thresholddensity $threshold]
						# set CCres [mdffi cc $resSel -mol $dens_mol -res $resolution -spacing $spacing -thresholddensity $threshold]
            if {$spacing != -1} {
             set CCres [mdffi cc $resSel -mol $dens_mol -res $resolution -spacing $spacing -thresholddensity $threshold]
            } else {
             set CCres [mdffi cc $resSel -mol $dens_mol -res $resolution -thresholddensity $threshold]
            }
						set str [string range "$CCres" 1 3]
						if {[string equal $str "nan"] || [string equal $str "NaN"]} {
							$resSel set beta -1.0
							} else {
								$resSel set beta $CCres
							}
							#rm mask.dx compare.dx
							puts $log "$res $CCres"
							unset resSel	
						}
					}

					if {$bb_only} {
						[atomselect $pdb all] writepdb $MOL-res-bb.pdb
					} else {
						[atomselect $pdb all] writepdb $MOL-res.pdb
					}

				}
				close $log

			}
