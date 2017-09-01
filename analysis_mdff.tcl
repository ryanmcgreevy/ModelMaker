namespace eval ::AnalysisMDFF {

	set packageDescription "AnalysisMDFF Plugin for RosettaVMD protocol"
	set version 0.1
}
package provide AnalysisMDFF $AnalysisMDFF::version

proc analyze_mdff {args} \
{
	return [eval ::AnalysisMDFF::analyze_mdff $args]
}

proc ::AnalysisMDFF::launch_gnuplot {} {
	global gnuplot
	set gnuplot [open "gnuplot" "w+"]
	fconfigure $gnuplot -buffering none -blocking 0
}


proc ::AnalysisMDFF::gnuplot {a} {
    global gnuplot
    fileevent $gnuplot writable [puts $gnuplot $a]
}

proc ::AnalysisMDFF::setGnuplotOutput {output} {
	::AnalysisMDFF::gnuplot "set term png size 2048,1536 font \"/usr/share/fonts/dejavu/DejaVuSans.ttf\" 14"
	::AnalysisMDFF::gnuplot "set output '$output.png'"
}


proc ::AnalysisMDFF::analyze_mdff {mol traj mapname res ch_seg} \
{
	global gnuplotexe
mol delete all	
set pdb [mol new ${mol}.psf]
mol addfile $traj.dcd type dcd first 0 last -1 step 1 waitfor all

set dName $mapname
set outName [lindex [molinfo $pdb get filename] 0 1]
#set res $res

set dens_mol [mol new $dName]

#mdff check -ccc -map $dName.dx -res $res -cccseltext "protein and noh" -cccfile ${outName}_ccc.txt
set cc_file [open ${outName}_ccc.txt w]
set numframes [molinfo $pdb get numframes]
puts "frames: $numframes"
for {set i 0} {$i < $numframes} {incr i} {
	puts $i
	# TODO: Ryan mdffi
	# set current_cc [mdffi cc [atomselect $pdb "protein and noh" frame $i] -mol $dens_mol -res $res -thresholddensity 1.0]
	# set current_cc 0
  set current_cc [mdffi cc [atomselect $pdb "protein and noh" frame $i] -mol $dens_mol -res $res -thresholddensity 0.1]
  #[mdff ccc [atomselect $pdb "protein and noh" frame $i] -i $dName -res $res]
	puts $cc_file "$i $current_cc"
}
close $cc_file

mdff check -mol $pdb -rmsd -rmsdseltext "backbone" -rmsdfile ${outName}_rmsd.txt

# delete density
mol delete $dens_mol

animate write pdb $outName-last.pdb beg [expr $numframes -1] end [expr $numframes -1] sel [atomselect $pdb all] 
set last_ccc [mdff ccc [atomselect $pdb "protein and noh" frame [expr $numframes -1]] -i $dName -res $res]
set f [open "last_ccc.txt" "w"]
puts $f $last_ccc
close $f

mol delete all
mol new $outName-last.pdb
foreach var $ch_seg {
	puts "name:[lindex $var 0] chain:[lindex $var 1] segname:[lindex $var 2]" 
	set chain [lindex $var 1]
    	set sel [atomselect [molinfo top] "segname [lindex $var 2]"]
    	$sel set chain $chain
}

# set tot [atomselect [molinfo top] "all" frame [expr $numframes -1]]
# $tot writepdb $outName-last.pdb

::AnalysisMDFF::launch_gnuplot
::AnalysisMDFF::gnuplot "reset"
::AnalysisMDFF::setGnuplotOutput "$outName-ccc"
::AnalysisMDFF::gnuplot "plot \"${outName}_ccc.txt\" using 1:2 w lines t \"CCC: $outName\""
exec gnuplot "gnuplot"

::AnalysisMDFF::launch_gnuplot
::AnalysisMDFF::gnuplot "reset"
::AnalysisMDFF::setGnuplotOutput "$outName-rmsd"
::AnalysisMDFF::gnuplot "plot \"${outName}_rmsd.txt\" using 1:2 w lines t \"RMSD: $outName\""
exec gnuplot "gnuplot"

}
