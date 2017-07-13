# vmd -dispdev text -e
#
#	HELIX-TURN
#	
#	in combination with python script regression.py
#	
#	- calculates linear regression of given selection 
#	and rotates a selection of interest around the calculated axis
#
#	by Maximilian Scheurer, Feb 2016

#Define Molecule
# set MOL Rpn8_mdff
# #Define De
# set Density Rpn8_mdff_4_6479_density 
# set Res 3.5
# set Threshold 1 
proc helixRotateMove {MOL seltext angleStep xLowerLimit xUpperLimit xStep Density Res Threshold inputfolder} {
	mol delete all
	# load molecule to for 3d linear regression
	mol new $inputfolder/$MOL.pdb


	# set selection for regression and write output files containing coordinates
	set sel [atomselect top "$seltext and backbone"]
	set fx [open "${inputfolder}/x.dat" w]
	set fy [open "${inputfolder}/y.dat" w]
	set fz [open "${inputfolder}/z.dat" w]
	foreach coord [$sel get {x y z}] {
		puts $fx [lindex $coord 0]
		puts $fy [lindex $coord 1]
		puts $fz [lindex $coord 2]
	}
	close $fx
	close $fy
	close $fz

	#cd ${inputfolder}
	#execute python script that writes file reg.txt with start and end points of averaging line
	#exec "/usr/bin/python2.6" "--version"
	#exec "setenv PYTHONPATH /home/mscheurer/local/lib/python2.6/site-packages"
	
	#does not work
	#exec python2 regression.py
	#wait forall
	#exec "cd -"
	source ${inputfolder}/reg.txt
	#draw the calculated line in vmd
	draw color red
	set startCoord [list $x1 $y1 $z1]
	set endCoord [list $x2 $y2 $z2]
	#draw cylinder $startCoord $endCoord radius 0.5

	#for {set i 0} {$i <= 180} {incr i 36} {
	#	#set selection for rotation around axis from startCoord to endCoord
	#	set rotSel [atomselect top "(resid 241 to 285)"]
	#	#set angle for rotation in DEGREES
	#	set angle $i
	#	$rotSel move [trans bond $startCoord $endCoord $angle deg]#
	#	set tot [atomselect top all]
	#	$tot writepdb $MOL-$angle-deg.pdb
	#}

	set outputfile [open ${MOL}-ccc.dat w+]
	set outfiles [list]
	for {set i $xLowerLimit} {$i <= $xUpperLimit} {incr i $xStep} {
		#set selection for move around and along regression axis from startCoord to endCoord
		set moveSel [atomselect top "$seltext"]
		# move along regression axis
		set vec [vecsub $endCoord $startCoord]
		set vNorm [vecnorm $vec]
		# set length to move along vec in Angstrom
		set dist $i
		set moveCoords [list [expr [lindex $vNorm 0]*$dist] [expr [lindex $vNorm 1]*$dist] [expr [lindex $vNorm 2]*$dist]]
		puts $moveCoords 
		$moveSel moveby $moveCoords
		set tot_mov [atomselect top "$seltext"]
	#	$tot_mov writepdb ${MOL}-${dist}.pdb

		for {set j 0} {$j < 360} {incr j $angleStep} {
			#set selection for rotation around axis from startCoord to endCoord
			set rotSel [atomselect top "$seltext"]
			#set angle for rotation in DEGREES
			set angle $j
			$rotSel move [trans bond $startCoord $endCoord $angle deg]
			set sel_ccc [atomselect top "$seltext"]		
			set tot [atomselect top "$seltext"]
			set filename ${MOL}-${dist}-${angle}-deg
			
			if {[file exists $filename] == 0} {
				file mkdir $filename
			}
			$tot writepdb $filename/${filename}.pdb
			lappend outfiles $filename
			#set ccc [mdff ccc $sel_ccc -i $inputfolder/$Density.dx -res $Res -threshold $Threshold]
			
			#puts $outputfile "$i $j $ccc"
			
		}

		#set tot [atomselect top all]
		#$tot writepdb ${MOL}-${dist}.pdb
	}

	close $outputfile

	#write new file to pdb (if desired)
	if {0} {
		set tot [atomselect top all]
		$tot writepdb $MOL-$angle-deg.pdb
	}

	#cleanup files (if desired)
	if {0} {
		rm x.dat y.dat z.dat reg.txt
	}
	return $outfiles
}
#quit
