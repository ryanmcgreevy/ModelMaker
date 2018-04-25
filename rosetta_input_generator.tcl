####################################################
# VMD Rosetta Input Generator
# Maximilian Scheurer, March 2016
####################################################


#################################
# create rosetta input files
# Maximilian Scheurer, March 2016
#################################

#### DEFAULT FRAGMENT PATH:
#set fragpath [pwd]

namespace eval ::RosettaInputGenerator {
    namespace export refine_with_rosetta
    # namespace export refine_sidechains_rosetta
    namespace export rosetta_abinitio

	# Set up Variable
	set version 0.1
	set packageDescription "RosettaInputGenerator Plugin"

    # Variable for the path of the script
    variable home [file join [pwd] [file dirname [info script]]]
}
package provide RosettaInputGenerator $RosettaInputGenerator::version

proc refine_with_rosetta {args} \
{
	return [eval ::RosettaInputGenerator::refine_with_rosetta $args]
}

proc rosetta_abinitio {args} \
{
	return [eval ::RosettaInputGenerator::rosetta_abinitio $args]
}

proc rosetta_insertion {args} \
{
	return [eval ::RosettaInputGenerator::rosetta_insertion $args]
}

proc rosetta_basic_refinement {args} \
{
	return [eval ::RosettaInputGenerator::rosetta_basic_refinement $args]
}

proc ::RosettaInputGenerator::refine_with_rosetta {jobname MOL mapname res score_dens nstruct cluster nPerTask configuration {cartesianSampler 0}} \
{
	###############################
	#	CONFIGURATION
	###############################
	set chains [lindex $configuration 0]
	set spans [lindex $configuration 1]
	set exclude [lindex $configuration 2]
	set constraints [lindex $configuration 3]
	set anchor_residue [lindex $constraints 0]
	set anchor_spans [lindex $constraints 1]

	set converted []
	set conv_anchor_spans []
	#foreach ch $chains {
	#	lappend converted $ch
	#}
  puts $anchor_spans
	foreach sp $spans {
		lappend converted $sp
	}
	foreach anch $anchor_spans {
		lappend conv_anchor_spans $anch
    lappend converted [list "Span" [lindex $anch 0] [lindex $anch 1] 0 0]
	}
	puts $converted

	###############################
	#	ROSETTA XML SCRIPT
	###############################
	set allMovers {}

	set tot ""
	#START SCRIPT
	append tot "<ROSETTASCRIPTS>\n"

	#LOAD SCRFXNS
	#puts [make_default_scfxn]
	append tot [::RosettaInputGenerator::make_default_scfxn]

	#START MOVERS
	append tot "<MOVERS>\n"

	set denssetup [::RosettaInputGenerator::make_density_setup]
	append tot $denssetup

	append tot [::RosettaInputGenerator::make_centroid_mover "centroid"]

	set min [::RosettaInputGenerator::make_cenmin cenmin 0 0 $converted]
	append tot $min

# {"Span" 75 96 1 1}
	set re_cst [::RosettaInputGenerator::make_relaxcart relaxcart 1 $converted]
	append tot $re_cst

	if {$cartesianSampler} {
		set cart [::RosettaInputGenerator::make_cart_sample [list cen5_50 $score_dens cen cen dens_soft auto density %%rms%% 200 0 0 25 4 7 25 "$exclude"]]
		append tot $cart
	}

	# append tot [::RosettaInputGenerator::make_full_atom_mover "fullatom"]

	# set fr [::RosettaInputGenerator::make_fastrelax relax 1 [list {"Chain" 1 0 0}]]
	# append tot $fr

	set cst [::RosettaInputGenerator::make_coord_cst [list coordcst 20 CA $anchor_residue 0] $conv_anchor_spans]
	append tot $cst

# "centroid" , "coordcst" "cenmin"
	lappend allMovers "setupdens" "loaddens" "coordcst" "relaxcart";#
	if {$cartesianSampler} {
		lappend allMovers "coordcst" "cen5_50" "coordcst" "relaxcart";
	}

	append tot "</MOVERS>\n"
	#END MOVERS

	#START PROTOCOLS

	set prot [::RosettaInputGenerator::make_protocol $allMovers]
	append tot $prot

	#END PROTOCOLS

	#append tot "<OUTPUT scorefxn=dens/> \n"
	append tot "</ROSETTASCRIPTS>\n"

	set f [open "$jobname.xml" w]
	file rename $jobname.xml $::MODELMAKER::workdir/setup-$jobname/
	puts $f $tot
	close $f

	###############################
	#	BATCH SUBMISSION SCRIPT
	###############################
	set bashscript [::RosettaInputGenerator::make_bash_script 25 1.5 $res $mapname.mrc $jobname $nstruct]

	set script [open "$jobname.sh" w]
	puts $script $bashscript
	close $script
}

# basic refinement without density
proc ::RosettaInputGenerator::rosetta_basic_refinement {jobname MOL nstruct cluster nPerTask configuration} \
{
	###############################
	#	CONFIGURATION
	###############################
	set chains [lindex $configuration 0]
	set spans [lindex $configuration 1]
	set exclude [lindex $configuration 2]
	set constraints [lindex $configuration 3]
	set anchor_residue [lindex $constraints 0]
	set anchor_spans [lindex $constraints 1]

	set converted []
	set conv_anchor_spans []
	foreach ch $chains {
		lappend converted $ch
	}
	foreach sp $spans {
		lappend converted $sp
	}
	foreach anch $anchor_spans {
		lappend conv_anchor_spans $anch
	}
	puts $converted

	exec mkdir -p "rosetta_input_$jobname"
	exec mkdir -p "rosetta_output_$jobname"
	###############################
	#	ROSETTA XML SCRIPT
	###############################
	set allMovers {}

	set tot ""
	#START SCRIPT
	append tot "<ROSETTASCRIPTS>\n"

	#LOAD SCRFXNS
	#puts [make_default_scfxn]
	# append tot [::RosettaInputGenerator::make_default_scfxn]

	#START MOVERS
	append tot "<MOVERS>\n"

	# set denssetup [::RosettaInputGenerator::make_density_setup]
	# append tot $denssetup

	# append tot [::RosettaInputGenerator::make_centroid_mover "centroid"]

	# set min [::RosettaInputGenerator::make_cenmin cenmin 0 0 $converted]
	# append tot $min

# {"Span" 75 96 1 1}
	# set re_cst [::RosettaInputGenerator::make_relaxcart relaxcart 1 $converted]
	# append tot $re_cst

	# if {$cartesianSampler} {
		# set cart [::RosettaInputGenerator::make_cart_sample [list cen5_50 $score_dens cen cen dens_soft auto density %%rms%% 200 0 0 25 4 7 25 "$exclude"]]
		# append tot $cart
	# }

	# append tot [::RosettaInputGenerator::make_full_atom_mover "fullatom"]

	# set fr [::RosettaInputGenerator::make_fastrelax relax 1 [list {"Chain" 1 0 0}]]
	# append tot $fr

	set cst [::RosettaInputGenerator::make_coord_cst [list coordcst 20 CA $anchor_residue 0] $conv_anchor_spans]
	append tot $cst

	set fr [make_fastrelax relax 1 $converted]
	append tot $fr

# "centroid" , "coordcst" "cenmin"
	lappend allMovers "coordcst" "relax";#

	append tot "</MOVERS>\n"
	#END MOVERS

	#START PROTOCOLS

	set prot [::RosettaInputGenerator::make_protocol $allMovers]
	append tot $prot

	#END PROTOCOLS

	#append tot "<OUTPUT scorefxn=dens/> \n"
	append tot "</ROSETTASCRIPTS>\n"

	set f [open "$jobname.xml" w]
	exec mv $jobname.xml rosetta_input_$jobname
	puts $f $tot
	close $f

	###############################
	#	BATCH SUBMISSION SCRIPT
	###############################
		# jobname nstruct
  set bashscript [::RosettaInputGenerator::make_abinitio_local_script $jobname $nstruct]
	
  set script [open "$jobname.sh" w]
	puts $script $bashscript
	close $script
}

proc ::RosettaInputGenerator::rosetta_abinitio {jobname MOL fragfiles nstruct cluster nPerTask test configuration chain_idents {prefix ""}} \
{
	###############################
	#	CONFIGURATION
	###############################
	set chains [lindex $configuration 0]
	set spans [lindex $configuration 1]
	set exclude [lindex $configuration 2]
	set constraints [lindex $configuration 3]
	set anchor_residue [lindex $constraints 0]
	set anchor_spans [lindex $constraints 1]

	set converted []
	set conv_anchor_spans []
	foreach ch $chains {
		lappend converted $ch
	}
	foreach sp $spans {
		lappend converted $sp
	}
	foreach anch $anchor_spans {
		lappend conv_anchor_spans $anch
	}
	puts $converted
	##############################

  ###############################
	#	ROSETTA XML SCRIPT
	###############################
	set allMovers {}

	set tot ""
	#START SCRIPT
	append tot "<ROSETTASCRIPTS>\n"

	#LOAD SCRFXNS
	#puts [make_default_scfxn]

	#append tot [make_default_scfxn]

	#RESIDUE SELECTORS
	# TODO: multiple chains!
	set resid_selector_list [list [list "Index" "fix" "$exclude"] ]
	foreach chain $chain_idents {
		puts $chain
		lappend resid_selector_list [list "Chain" "Chain$chain" "$chain"]
	}
	append tot [make_resid_selectors $resid_selector_list]

	#START MOVERS
	append tot "<MOVERS>\n"

	# set denssetup [make_density_setup]
	# append tot $denssetup

	append tot [make_centroid_mover "centroid"]

	append tot [make_rigid_chunk chunk fix $MOL fix centroid $jobname]

	# TODO: multiple chains!
	# check if right fragfiles are assigned to the right chain!
	set frag {}
	set counter 0
	foreach chain $chain_idents {
		lappend frag [list "$::MODELMAKER::workdir/setup-$jobname/[lindex $fragfiles $counter 0]" "$::MODELMAKER::workdir/setup-$jobname/[lindex $fragfiles $counter 1]" "Chain$chain"]
		incr counter
	}

	append tot [make_abscript_mover abinitio 2 $frag]

	append tot [make_full_atom_mover "fullatom"]

	append tot [make_environment env true {chunk} {abinitio}]

	set fr [make_fastrelax relax 1 $converted]
	append tot $fr

	# set cst [make_coord_cst {coordcst 20 CA 95 0} {{90 389}}]
	# append tot $cst

	lappend allMovers "centroid" "env" "fullatom" "relax"

	#lappend allMovers "coordcst"


	#lappend allMovers "relax"


	#lappend allMovers "cen5_50"


	#lappend allMovers "relaxcart"


	#lappend allMovers "cenmin"

	append tot "</MOVERS>\n"
	#END MOVERS

	#START PROTOCOLS

	set prot [make_protocol $allMovers]
	append tot $prot

	#END PROTOCOLS

	#append tot "<OUTPUT scorefxn=dens/> \n"
	append tot "</ROSETTASCRIPTS>\n"

	set f [open "$::MODELMAKER::workdir/run-$jobname/$jobname.xml" w]
	puts $f $tot
	close $f

	###############################
	#	BATCH SUBMISSION SCRIPT
	###############################
	if {$test} {
		set bashscript [make_abinitio_test_script $jobname $nstruct]
	} elseif {!$test && !$cluster} {
		set bashscript [make_abinitio_local_script $jobname $nstruct $prefix]
	}


	set script [open "$::MODELMAKER::workdir/run-$jobname/$jobname.sh" w]
	puts $script $bashscript
	close $script
}


proc ::RosettaInputGenerator::rosetta_insertion {jobname MOL fragfiles fasta nstruct cluster nPerTask configuration} \
{
	###############################
	#	CONFIGURATION
	###############################
	set chains [lindex $configuration 0]
	set spans [lindex $configuration 1]
	set exclude [lindex $configuration 2]
	# set constraints [lindex $configuration 3]
	# set anchor_residue [lindex $constraints 0]
	# set anchor_spans [lindex $constraints 1]

	# foreach ch $chains {
	# 	lappend converted $ch
	# }
	# foreach sp $spans {
	# 	lappend converted $sp
	# }
	# foreach anch $anchor_spans {
	# 	lappend conv_anchor_spans $anch
	# }
	# puts $converted
	##############################

	puts $exclude
	set exclusions [split $exclude ","]
	set f [open "$::MODELMAKER::workdir/setup-$jobname/protein.rigid" "w"]
	foreach ex $exclusions {
		set xsplit [split $ex "-"]
		puts $f "RIGID [lindex $xsplit 0] [lindex $xsplit 1] 0 0 0"
	}
	close $f

	set f [open "$::MODELMAKER::workdir/setup-$jobname/input.tpb" "w"]
	puts $f "CLAIMER RigidChunkClaimer"
	puts $f "REGION_FILE $::MODELMAKER::workdir/setup-$jobname/protein.rigid"
	puts $f "PDB_FILE $::MODELMAKER::workdir/setup-$jobname/$MOL.pdb"
	puts $f "END_CLAIMER"
	close $f

	set frag3 [lindex $fragfiles 1]
	set frag9 [lindex $fragfiles 0]
  set bashscript [::RosettaInputGenerator::make_insertion_local_script $jobname $MOL $nstruct $frag3 $frag9 $fasta]
	set script [open "$::MODELMAKER::workdir/run-$jobname/$jobname.sh" w]
	puts $script $bashscript
	close $script

}

#<CartesianSampler name=cen5_50 automode_scorecut=-0.5 scorefxn=cen mcscorefxn=cen
#fascorefxn=dens_soft strategy="auto" fragbias="density" rms=%%rms%% ncycles=200 fullatom=0
#bbmove=0 nminsteps=25 temp=4 fraglens=7 nfrags=25 residues_to_exclude="50-116"/>
# usage


#<Chain number=1 chi=0 bb=0 />
#<Span begin=1 end=61 chi=1 bb=1 />
# <Jump number=1 setting=0/>
# each component has to be like:
#
# name <args>, e.g.
# make_movemap [list {"Chain" 1 0 1} {"Span" 1 61 0 1} {"Jump" 1 0}]
#
proc ::RosettaInputGenerator::make_movemap {components} \
{
	set out ""
	append out "<MoveMap> \n"
	foreach comp $components {
		set name [lindex $comp 0]
		switch -exact -- $name {
			"Chain" {
				if {[llength $comp] != 4} {
					wrong_input "chain length"
				}
				append out "<Chain number=\"[lindex $comp 1]\" chi=\"[lindex $comp 2]\" bb=\"[lindex $comp 3]\" /> \n"
			}
			"Span" {
				if {[llength $comp] != 5} {
					wrong_input "span length"
				}
				append out "<Span begin=\"[lindex $comp 1]\" end=\"[lindex $comp 2]\" chi=\"[lindex $comp 3]\" bb=\"[lindex $comp 4]\" /> \n"
			}
			"Jump" {
				if {[llength $comp] != 3} {
					wrong_input
				}
				append out "<Jump number=\"[lindex $comp 1]\" setting=\"[lindex $comp 2]\" /> \n"
			}
			default {
				::RosettaInputGenerator::wrong_input
			}
		}
	}
	append out "</MoveMap> \n"
	return $out
}

#<CoordinateCst name=coordcst stddev=20 atom=CA anchor_res=90 jump=0  >
#<Span begin=66 end=104 />
#</CoordinateCst>
# usage:
# make_coord_cst {coordcst 20 CA 90 0} {{1 61 0 1}}
#
proc ::RosettaInputGenerator::make_coord_cst {config spans} \
{
	set out ""
	if {[llength $config] != 5} {
		::RosettaInputGenerator::wrong_input
	}
	append out "<CoordinateCst name=\"[lindex $config 0]\" stddev=\"[lindex $config 1]\" atom=\"[lindex $config 2]\" anchor_res=\"[lindex $config 3]\" jump=\"[lindex $config 4]\"> \n"
	foreach span $spans {
		if {[llength $span] != 2} {
			::RosettaInputGenerator::wrong_input
		}
		append out "<span begin=\"[lindex $span 0]\" end=\"[lindex $span 1]\" /> \n"
	}
	append out "</CoordinateCst> \n"
	return $out
}

#<FastRelax name="relax" repeats="1">
#			<MoveMap>
#				<Chain number=1 chi=1 bb=0 /> Defines degrees of freedom for chain A
#			</MoveMap>
#		</FastRelax>

proc ::RosettaInputGenerator::make_fastrelax {name repeats movemap_comps} \
{
	set out ""
	append out "<FastRelax name=\"$name\" repeats=\"$repeats\">\n"
	append out [::RosettaInputGenerator::make_movemap $movemap_comps]
	append out "</FastRelax>\n"
	return $out
}


# <AbscriptMover name="abinitio" cycles="1" >
#         <Fragments large_frags="../rosetta_input/rpt4_human_1-114_frag9" small_frags="../rosetta_input/rpt4_human_1-114_frag3" selector="ChainA" />
#         <Fragments large_frags="../rosetta_input/rpt5_human_1-168_frag9" small_frags="../rosetta_input/rpt5_human_1-168_frag3" selector="ChainA" />
#       </AbscriptMover>
#	USAGE:
# name cycles [list {frag9 frag3 ChainA}]
proc ::RosettaInputGenerator::make_abscript_mover {name cycles fraglist} \
{
	set out ""
	append out "<AbscriptMover name=\"$name\" cycles=\"$cycles\" >\n"
	foreach frag $fraglist {
		set name [lindex $frag 0]
		append out "<Fragments large_frags=\"[lindex $frag 0]\" small_frags=\"[lindex $frag 1]\" selector=\"[lindex $frag 2]\" /> \n"
	}
	append out "</AbscriptMover> \n"
	return $out
}


# <Environment name="env" auto_cut="true">
#         <Register mover="chunk" />
#         <Apply mover="abinitio" />
#       </Environment>
proc ::RosettaInputGenerator::make_environment {name autocut registerList applyList} \
{
	set out ""
	append out "<Environment name=\"$name\" auto_cut=\"$autocut\"> \n"
	foreach reg $registerList {
		append out "<Register mover=\"$reg\" />\n"
	}
	foreach app $applyList {
		append out "<Apply mover=\"$app\" />\n"
	}
	append out "</Environment>\n"
	return $out
}


#<RigidChunkCM name="chunk" region_selector="fix" template="../full_length_model/rpt4_5_human_complete.pdb" selector="fix" apply_to_template="centroid" />
proc ::RosettaInputGenerator::make_rigid_chunk {name regselector template selector apply_to_template jobname} \
{
	set out "<RigidChunkCM name=\"$name\" region_selector=\"$regselector\" template=\"$::MODELMAKER::workdir/setup-$jobname/$template\" selector=\"$selector\" apply_to_template=\"$apply_to_template\" />\n"
	return $out
}


#<CartesianSampler name=cen5_50 automode_scorecut=-0.5 scorefxn=cen mcscorefxn=cen
#fascorefxn=dens_soft strategy="auto" fragbias="density" rms=%%rms%% ncycles=200 fullatom=0
#bbmove=0 nminsteps=25 temp=4 fraglens=7 nfrags=25 residues_to_exclude="50-116"/>
# usage

proc ::RosettaInputGenerator::make_cart_sample {arg} \
{
	set out ""
	if {[llength $arg] < 15} {
		::RosettaInputGenerator::wrong_input
	}
	if {[llength $arg] == 16} {
		append out "<CartesianSampler name=\"[lindex $arg 0]\" automode_scorecut=\"[lindex $arg 1]\" scorefxn=\"[lindex $arg 2]\" mcscorefxn=\"[lindex $arg 3]\" fascorefxn=\"[lindex $arg 4]\" strategy=\"[lindex $arg 5]\" fragbias=\"[lindex $arg 6]\" rms=\"[lindex $arg 7]\" ncycles=\"[lindex $arg 8]\" fullatom=\"[lindex $arg 9]\" bbmove=\"[lindex $arg 10]\" nminsteps=\"[lindex $arg 11]\" temp=\"[lindex $arg 12]\" fraglens=\"[lindex $arg 13]\" nfrags=\"[lindex $arg 14]\" residues_to_exclude=\"[lindex $arg 15]\" />"
	} else {
		append out "<CartesianSampler name=\"[lindex $arg 0]\" automode_scorecut=\"[lindex $arg 1]\" scorefxn=\"[lindex $arg 2]\" mcscorefxn=\"[lindex $arg 3]\" fascorefxn=\"[lindex $arg 4]\" strategy=\"[lindex $arg 5]\" fragbias=\"[lindex $arg 6]\" rms=\"[lindex $arg 7]\" ncycles=\"[lindex $arg 8]\" fullatom=\"[lindex $arg 9]\" bbmove=\"[lindex $arg 10]\" nminsteps=\"[lindex $arg 11]\" temp=\"[lindex $arg 12]\" fraglens=\"[lindex $arg 13]\" nfrags=\"[lindex $arg 14]\" />"
	}
	append out "\n"
	return $out
}


#<FastRelax name=relaxcart scorefxn=dens repeats=1 cartesian=1>
#			<MoveMap>
#				<Chain number=1 chi=1 bb=0 /> Defines degrees of freedom for chain A
#				<!-- <Jump number=1 setting=0/> Turns off minimization between chain A and B -->
#				<Span begin=1 end=61 chi=1 bb=1/> Defines degrees of freedom for resid x to y Rosetta internal numbering starting from 1
#			</MoveMap>
#		</FastRelax>
proc ::RosettaInputGenerator::make_relaxcart {name repeats movemap_comps} \
{
	set out ""
	append out "<FastRelax name=\"$name\" scorefxn=\"dens\" repeats=\"$repeats\" cartesian=\"1\"> \n"
	append out [::RosettaInputGenerator::make_movemap $movemap_comps]
	append out "</FastRelax> \n"
	return $out
}


#<MinMover name=cenmin scorefxn=cen type=lbfgs_armijo_nonmonotone max_iter=200 tolerance=0.00001 chi=0 bb=0>
#			<MoveMap>
#				<Chain number=1 chi=0 bb=0 />
#				<Span begin=1 end=61 chi=1 bb=1 />
#			</MoveMap>
#		</MinMover>

proc ::RosettaInputGenerator::make_cenmin {name chi bb movemap_comps} \
{
	set out ""
	append out "<MinMover name=\"$name\" scorefxn=\"cen\" type=\"lbfgs_armijo_nonmonotone\" max_iter=\"200\" tolerance=\"0.00001\" chi=\"$chi\" bb=\"$bb\"> \n"
	append out [::RosettaInputGenerator::make_movemap $movemap_comps]
	append out "</MinMover>\n"
	return $out
}

#<Index name="fix" resnums="21-114,167-211,235-282" /> Defines the selector fix for all atoms in the range resnum_a-resnum_b,resnum_c-resnum_d
#      <Chain name="ChainA" chains="A" /> Defines the selector ChainA for all atoms with chain ID A in the input file
#      <Chain name="ChainB" chains="B" />
proc ::RosettaInputGenerator::make_resid_selectors {components} \
{
	set out ""
	append out "<RESIDUE_SELECTORS>\n"
	foreach comp $components {
		set name [lindex $comp 0]
		switch -exact -- $name {
			"Index" {
				append out "<Index name=\"[lindex $comp 1]\" resnums=\"[lindex $comp 2]\" />\n"
			}
			"Chain" {
				append out "<Chain name=\"[lindex $comp 1]\" chains=\"[lindex $comp 2]\" /> \n"
			}
			default {
				::RosettaInputGenerator::wrong_input
			}
		}
	}
	append out "</RESIDUE_SELECTORS>\n"
	return $out
}

proc ::RosettaInputGenerator::make_protocol {movernames} \
{
	set out ""
	append out "<PROTOCOLS> \n"
	foreach mover $movernames {
		append out "<Add mover=\"$mover\"/> \n"
	}
	append out "</PROTOCOLS> \n"
	return $out
}


#DEFAULT SCOREFXNS

proc ::RosettaInputGenerator::make_default_scfxn {} \
{
	return "<SCOREFXNS>
		<ScoreFunction name=\"cen\" weights=\"score4_smooth_cart\">
			<Reweight scoretype=\"elec_dens_fast\" weight=\"20\"/>
		</ScoreFunction>
		<ScoreFunction name=\"dens_soft\" weights=\"soft_rep\">
			<Reweight scoretype=\"cart_bonded\" weight=\"0.5\"/>
			<Reweight scoretype=\"pro_close\" weight=\"0.0\"/>
			<Reweight scoretype=\"elec_dens_fast\" weight=\"%%denswt%%\"/>
		</ScoreFunction>
		<ScoreFunction name=\"dens\" weights=\"talaris2014_cart\">
			<Reweight scoretype=\"elec_dens_fast\" weight=\"%%denswt%%\"/>
			<Set scale_sc_dens_byres=\"R:0.76,K:0.76,E:0.76,D:0.76,M:0.76,C:0.81,Q:0.81,H:0.81,N:0.81,T:0.81,S:0.81,Y:0.88,W:0.88,A:0.88,F:0.88,P:0.88,I:0.88,L:0.88,V:0.88\"/>
		</ScoreFunction>
	</SCOREFXNS> \n"
}

# density setup

proc ::RosettaInputGenerator::make_density_setup {} \
{
	return "<SetupForDensityScoring name=\"setupdens\"/> \n <LoadDensityMap name=\"loaddens\" mapfile=\"%%map%%\"/>\n"
}


proc ::RosettaInputGenerator::make_bash_script {denswt rms res map jobname nstruct} \
{
	global rosettapath
	global rosettaDBpath
	global platform
  
  if { [string match "*mpi*" $platform] } {
    set mpi_args "mpiexec -np $::MODELMAKER::MPINP"      
  } else {
    set mpi_args ""
  }
	return "#!/bin/bash
	JOBNAME=\"\${1}\"
	MOL=\"\${2}\"
	if \[ \-z \"\$1\" \]\; then
  		echo Need job name!
  		exit
	fi
	$mpi_args $rosettapath/rosetta_scripts.$platform \\
			-database $rosettaDBpath \\
			-nstruct $nstruct \\
		-parser::script_vars denswt=$denswt rms=$rms reso=$res map=$::MODELMAKER::workdir/setup-$jobname/$map testmap=$::MODELMAKER::workdir/setup-$jobname/$map \\
		-edensity::mapreso $res \\
	        -out::prefix \${JOBNAME}_ \\
          -out:path:pdb $::MODELMAKER::workdir/run-$jobname/pdb_out/ \\
          -out:path:score $::MODELMAKER::workdir/run-$jobname/sc_out/ \\
		-s $::MODELMAKER::workdir/setup-$jobname/\${MOL} \\
	        -parser::protocol $::MODELMAKER::workdir/setup-$jobname/$jobname.xml \\
	        -parser:view \\
		-restore_talaris_behavior \\
	        -overwrite \\
	        -ignore_zero_occupancy false

	"
}

proc ::RosettaInputGenerator::make_abinitio_test_script {jobname nstruct} \
{
	global rosettapath
	global rosettaDBpath
	global platform
  
  if { [string match "*mpi*" $platform] } {
    set mpi_args "mpiexec -np $::MODELMAKER::MPINP"      
  } else {
    set mpi_args ""
  }
	return "

#!/bin/bash
JOBNAME=\"\${1}\"
MOL=\"\${2}\"

if \[ -z \"\$1\" \]\; then
  echo Need job name!
  exit
fi

$mpi_args $rosettapath/rosetta_scripts.$platform \\
    -database $rosettaDBpath \\
	-nstruct $nstruct \\
	-run:test_cycles \\
    -out::prefix \${JOBNAME}_ \\
	-out:path:pdb $::MODELMAKER::workdir/run-$jobname/pdb_out/ \\
  -out:path:score $::MODELMAKER::workdir/run-$jobname/sc_out/ \\
	-s $::MODELMAKER::workdir/setup-$jobname/\${MOL} \\
    -parser::protocol $::MODELMAKER::workdir/run-$jobname/$jobname.xml \\
    -parser:view \\
    -ignore_zero_occupancy false\\
    -overwrite

#mv \${JOBNAME}*.pdb ./pdb_out/
#mv *.sc ./sc_out/

"
}


proc ::RosettaInputGenerator::make_insertion_local_script {jobname mol nstruct frag3 frag9 fasta} \
{
	global rosettapath
	global rosettaDBpath
	global platform
  if { [string match "*mpi*" $platform] } {
    set mpi_args "mpiexec -np $::MODELMAKER::MPINP"      
  } else {
    set mpi_args ""
  }
	return "
#!/bin/bash

$mpi_args $rosettapath/minirosetta.$platform \\
	-run::shuffle \\
	-abinitio::close_loops \\
	-short_frag_cycles 2 \\
	-scored_frag_cycles 2 \\
	-non_ideal_loop_closing \\
	-alternative_closure_protocol \\
	-fast_loops:window_accept_ratio .01 \\
	-fast_loops:nr_scored_sampling_passes 4 \\
	-fast_loops:min_breakout_good_loops 5 \\
	-fast_loops:min_breakout_fast_loops 80 \\
	-fast_loops:min_fast_loops 3 \\
	-fast_loops:min_good_loops 0 \\
	-fast_loops:nr_scored_fragments 20 \\
	-fast_loops:vdw_delta 0.5 \\
	-fast_loops:give_up 1000 \\
	-random_grow_loops_by  4 \\
	-increase_cycles 1 \\
	-jumps:ramp_chainbreaks \\
	-overlap_chainbreak \\
	-relax::fast \\
	-relax::default_repeats 1 \\
	-out::shuffle_nstruct $nstruct \\
	-out::prefix ${jobname}_${mol}_ \\
	-out:path:pdb $::MODELMAKER::workdir/run-$jobname/pdb_out/ \\
  -out:path:score $::MODELMAKER::workdir/run-$jobname/sc_out/ \\
  -run::protocol broker \\
	-in:file:fasta $::MODELMAKER::workdir/setup-$jobname/$fasta.fasta \\
	-broker:setup $::MODELMAKER::workdir/setup-$jobname/input.tpb \\
	-frag3 $::MODELMAKER::workdir/setup-$jobname/$frag3 \\
	-frag9 $::MODELMAKER::workdir/setup-$jobname/$frag9 \\
	-nstruct $nstruct \\
	-overwrite

"
}

proc ::RosettaInputGenerator::make_abinitio_local_script {jobname nstruct {prefix ""}} \
{
	global rosettapath
	global rosettaDBpath
	global platform
  
  if { [string match "*mpi*" $platform] } {
    set mpi_args "mpiexec -np $::MODELMAKER::MPINP"      
  } else {
    set mpi_args ""
  }
	return "
#!/bin/bash
JOBNAME=\"\${1}\"
MOL=\"\${2}\"

if \[ -z \"\$1\" \]\; then
  echo Need job name!
  exit
fi

$mpi_args $rosettapath/rosetta_scripts.$platform \\
    -database $rosettaDBpath \\
	  -nstruct $nstruct \\
    -out::prefix ${prefix}\${JOBNAME}_ \\
	  -out:path:pdb $::MODELMAKER::workdir/run-$jobname/pdb_out/ \\
    -out:path:score $::MODELMAKER::workdir/run-$jobname/sc_out/ \\
    -s $::MODELMAKER::workdir/setup-$jobname/\${MOL} \\
    -parser::protocol $::MODELMAKER::workdir/run-$jobname/$jobname.xml \\
    -parser:view \\
    -ignore_zero_occupancy false\\
    -overwrite


"
}

# RESIDUE TYPE MOVERS

proc ::RosettaInputGenerator::make_centroid_mover {name} \
{
	return "<SwitchResidueTypeSetMover name=\"$name\" set=\"centroid\" /> \n"
}

proc ::RosettaInputGenerator::make_full_atom_mover {name} \
{
	return "<SwitchResidueTypeSetMover name=\"$name\" set=\"fa_standard\" />\n"
}

proc ::RosettaInputGenerator::wrong_input {{str ""} args} \
{
	puts "wrong input!"
	puts $str
	exit
}
