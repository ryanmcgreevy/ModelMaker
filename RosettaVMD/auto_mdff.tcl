####################################################
# Customized Auto MDFF
# Maximilian Scheurer, March 2016
####################################################
namespace eval ::AutoMDFF {
# Set up Variable
  namespace export auto_mdff_init
  set version 0.1
  set packageDescription "AutoMDFF Plugin for RosettaVMD protocol"
 
    # Variable for the path of the script
    variable home [file join [pwd] [file dirname [info script]]]
}
package provide AutoMDFF $AutoMDFF::version

proc auto_mdff_init {args} \
{
	return [eval ::AutoMDFF::auto_mdff_init $args]
}

proc ::AutoMDFF::auto_mdff_init {jobname MOL mapname fixedselection gscale minSteps num ch_seg mutations topdir topfiles parfiles} \
{
	package require MakePsf
	set short $jobname

	#set EMDB_ID $::EMDB_ID

	#set cutoff $::cutoff

	#protein restraints
	set seltext $fixedselection


	#IMD
	#set ::IMD on
	#set ::IMDPort 3000
	#set ::IMDfreq 1

	set gridscale $gscale
	#set num $::num
	#set minSteps $::minSteps

	set outname mdff_$short
	#IMD
	set IMD 0
	#set IMDPort $::IMDPort
	#set IMDfreq $::IMDfreq

	package require ssrestraints
	package require cispeptide
	package require chirality
	package require mdff
	package require autopsf

	#get density
	#mdff griddx -i $emdbName.ccp4 -o ${emdbName}_potential.dx
	#mdff griddx -i ${emdbName}_potential.dx -o ${emdbName}_density.dx

	#crop density
	#mol new $MOL.pdb 

	# set sel [atomselect top all] 
	# volmap mask $sel -o mask.dx -cutoff $cutoff 
	# volutil -mult ${emdbName}_density.dx mask.dx -o ${MOL}_${cutoff}_${EMDB_ID}_density.dx

	#make psf
	#mol delete all
	#mol new $MOL.pdb

	#make _potential
	mdff griddx -i $mapname.dx -o ${mapname}_potential.dx

	auto_makepsf $MOL $topdir $topfiles $ch_seg $mutations
	#autopsf -mol top -top ../top_all27_prot_lipid_na.inp


	#make grid
	mdff gridpdb -psf ${MOL}.psf -pdb ${MOL}-psfout.pdb -o ${MOL}-psfout-grid.pdb

	#make ssrestraints
	ssrestraints -psf ${MOL}.psf -pdb ${MOL}-psfout.pdb -o ${MOL}-extrabonds.txt -hbonds

	#make cispeptide
	mol delete all
	mol new ${MOL}.psf
	mol addfile ${MOL}-psfout.pdb
	cispeptide restrain -o ${MOL}-cispeptide.txt

	#chirality
	chirality restrain -o ${MOL}-chirality.txt

	#fix pdb
	set all [atomselect top all]
	$all set occupancy 0.0
	set restraint [atomselect top "$seltext"]
	$restraint set occupancy 1.0
	$all writepdb fixed.pdb

	set pfiles []
	foreach par $parfiles {
		lappend pfiles $topdir/$par
	}

	puts [pwd]
	#setup mdff
	cd [pwd]
	if {$IMD} {
		#mdff setup -o [pwd]/$::outname -psf ${MOL}_autopsf.psf -pdb ${MOL}_autopsf.pdb -griddx ${MOL}_${cutoff}_${EMDB_ID}_density.dx -gridpdb ${MOL}_autopsf-grid.pdb -extrab [list ${MOL}-extrabonds.txt ${MOL}-cispeptide.txt ${MOL}-chirality.txt] -gscale $gridscale -minsteps $minSteps -numsteps $num --imd -imdport $IMDPort -imdfreq $IMDfreq -fixpdb fixed.pdb
		} else {
			mdff setup -o $outname -psf ${MOL}.psf -pdb ${MOL}-psfout.pdb -griddx ${mapname}_potential.dx -gridpdb ${MOL}-psfout-grid.pdb -extrab [list ${MOL}-extrabonds.txt ${MOL}-cispeptide.txt ${MOL}-chirality.txt] -gscale $gridscale -minsteps $minSteps -numsteps $num -fixpdb fixed.pdb -dir [pwd] -parfiles [list $pfiles]
	}
}

# config item: gscale smoothr temperature timesteps
# config: [list {} {} {} {}]
proc auto_mdff_casc_init {jobname MOL mapname fixedselection config minSteps dcdFreq ch_seg mutations topdir topfiles parfiles inputfolder gridconfig} \
{
	package require MakePsf
	set short $jobname

	#number of items in config determines number of cascade runs
	set cascades [llength $config]


	#set EMDB_ID $::EMDB_ID

	#set cutoff $::cutoff

	#protein restraints
	set seltext $fixedselection


	#IMD
	#set ::IMD on
	#set ::IMDPort 3000
	#set ::IMDfreq 1

	#set gridscale $gscale
	#set num $::num
	#set minSteps $::minSteps

	set IMD 0
	#set IMDPort $::IMDPort
	#set IMDfreq $::IMDfreq

	#source ../makepsf.tcl

	package require ssrestraints
	package require cispeptide
	package require chirality
	package require mdff
	package require autopsf

	#get density
	#mdff griddx -i $emdbName.ccp4 -o ${emdbName}_potential.dx
	#mdff griddx -i ${emdbName}_potential.dx -o ${emdbName}_density.dx

	#crop density
	#mol new $MOL.pdb 

	# set sel [atomselect top all] 
	# volmap mask $sel -o mask.dx -cutoff $cutoff 
	# volutil -mult ${emdbName}_density.dx mask.dx -o ${MOL}_${cutoff}_${EMDB_ID}_density.dx

	#make psf
	#mol delete all
	#mol new $MOL.pdb

	#make _potential
	
	#mdff griddx -i $mapname.dx -o ${mapname}_potential.dx

	auto_makepsf $MOL $topdir $topfiles $ch_seg mutations
	#autopsf -mol top -top ../top_all27_prot_lipid_na.inp


	#make grid
	mdff gridpdb -psf ${MOL}.psf -pdb ${MOL}-psfout.pdb -o ${MOL}-psfout-grid.pdb

	#modify gridpdb, resmap is supposed to be in [pwd]/..
	if {$gridconfig != 0} {
		set resmapmap [lindex $gridconfig 0]
		set min [lindex $gridconfig 1]
		set max [lindex $gridconfig 2]
		set stepwidth [lindex $gridconfig 3]
		set factor_min [lindex $gridconfig 4]
		set factor_max [lindex $gridconfig 5]
		modify_gridpdb ${MOL}-psfout-grid ../$resmapmap $min $max $stepwidth $factor_min $factor_max
	}

	set gridpdb ${MOL}-psfout-grid

	#make ssrestraints
	ssrestraints -psf ${MOL}.psf -pdb ${MOL}-psfout.pdb -o ${MOL}-extrabonds.txt -hbonds

	#make cispeptide
	mol delete all
	mol new ${MOL}.psf
	mol addfile ${MOL}-psfout.pdb
	cispeptide restrain -o ${MOL}-cispeptide.txt

	#chirality
	chirality restrain -o ${MOL}-chirality.txt

	#fix pdb
	set all [atomselect top all]
	$all set occupancy 0.0
	set restraint [atomselect top "$seltext"]
	$restraint set occupancy 1.0
	$all writepdb fixed.pdb
	set fixpdb [pwd]/fixed

	puts [pwd]
	#setup mdff
	# cd [pwd]
	# if {$IMD} {
	# 	#mdff setup -o [pwd]/$::outname -psf ${MOL}_autopsf.psf -pdb ${MOL}_autopsf.pdb -griddx ${MOL}_${cutoff}_${EMDB_ID}_density.dx -gridpdb ${MOL}_autopsf-grid.pdb -extrab [list ${MOL}-extrabonds.txt ${MOL}-cispeptide.txt ${MOL}-chirality.txt] -gscale $gridscale -minsteps $minSteps -numsteps $num --imd -imdport $IMDPort -imdfreq $IMDfreq -fixpdb fixed.pdb
	# 	} else {
	# 		mdff setup -o [pwd]/$outname -psf ${MOL}.psf -pdb ${MOL}-psfout.pdb -griddx ${mapname}_potential.dx -gridpdb ${MOL}-psfout-grid.pdb -extrab [list ${MOL}-extrabonds.txt ${MOL}-cispeptide.txt ${MOL}-chirality.txt] -gscale $gridscale -minsteps $minSteps -numsteps $num -fixpdb fixed.pdb -dir [pwd] -parfiles $topdir/$parfile
	# }
	puts "Configuring MDFF."

	set step 1

	set itemp 300
	# config item: gscale smoothr temperature timesteps
	foreach item $config {
		set gscale [lindex $item 0]
		set ftemp [lindex $item 2]
		set timesteps [lindex $item 3]
		set file [open "mdff_${MOL}_step$step.namd" w+]
		if {$step > 1} {
			set minSteps 0
		}
		generateMDFF $MOL $dcdFreq $minSteps $timesteps $gscale $itemp $ftemp $fixpdb $gridpdb $inputfolder $topdir $parfiles $file $step
		close $file
		incr step
	}
}

proc generateMDFF {MOL dcdFreq minsteps runsteps gscale itemp ftemp fixpdb gridpdb inputfolder parfolder parfiles file step} {
	puts $file "###  Cascade MDFF Customized to RosettaVMD protocol"
	puts $file "set inputfolder $inputfolder"
	puts $file "set dcdFreq $dcdFreq"

	puts $file "set PSFFILE $MOL.psf"
	puts $file "set PDBFILE ${MOL}-psfout.pdb"
	puts $file "set GRIDPDB $gridpdb.pdb"
	puts $file "set GBISON 0"
	puts $file "set DIEL 80"
	puts $file "set SCALING_1_4 1.0"
	puts $file "set ITEMP $itemp"
	puts $file "set FTEMP $ftemp"
	puts $file "set GRIDFILE $inputfolder/density_grids/griddx$step.dx"
	puts $file "set GSCALE $gscale"
	puts $file "set EXTRAB \{${MOL}-extrabonds.txt ${MOL}-cispeptide.txt ${MOL}-chirality.txt\}"
	puts $file "set CONSPDB 0"
	puts $file "set FIXPDB $fixpdb.pdb"
	puts $file "set GRIDON 1"
	if {$step > 1} {
		set before [expr $step - 1]
		set INPUTNAME mdff_${MOL}_step$before
		puts $file "BinVelocities $INPUTNAME.restart.vel"
		puts $file "BinCoordinates $INPUTNAME.restart.coor"
	} else {
		puts $file "temperature \$ITEMP"
	}
	puts $file ""
 
	puts $file "set OUTPUTNAME mdff_${MOL}_step$step"
	 
	puts $file "set TS $runsteps"
	puts $file "set MS $minsteps"
 
	puts $file "set MARGIN 0"
	
	puts $file "####################################"

	puts $file "structure \$PSFFILE"
	puts $file "coordinates \$PDBFILE"
	
	puts $file "paraTypeCharmm on"
	foreach par $parfiles {
		puts $file "parameters $parfolder/$par"
	}

	puts $file "source $inputfolder/mdff_template_casc.namd"
}
