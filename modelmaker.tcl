# general note: all underlying 'exec unix calls' should be changed to
#use built in tcl commands (e.g., 'glob' to list files instead of 'exec ls -1v')
#because it should not assume the presence or behavior of those commands, especially
#since people (like me) have some of those commands aliased.
package require RosettaVMD

package provide modelmaker 0.1

namespace eval ::MODELMAKER {
  variable defaultGapfindSel "all"
  #this should get auto-set depending on OS, but
  #should maybe be changed to not rely on this because it would be
  #best to be independent of any future naming conventions of rosetta. Use wildcard instead?
  switch $tcl_platform(os) {
    "Darwin" {
      variable rosettaEXE "macosclangrelease"
    }
    "Linux" {
      variable rosettaEXE "linuxclangrelease"
    }
    default {
      variable rosettaEXE "Unrecognized"
    }
  }
  #this environment variable needs to be set by user somehow once, probably
  #in .vmdrc but also exposed in eventual GUI. If not set, whether running through command line
  #or gui, if not present, perhaps pop up a dialog with a Browse button to allow user to get the path and
  #proceed...
  if {[info exists env(ROSETTADB)]} {
    variable rosettadbpath $env(ROSETTADB)
  } else {
    variable rosettadbpath ""
  }

  #currently set in .vmdrc like ROSETTADB, !!HOWEVER!!
  #this should not be necessary to specify if we tell the user that
  #the rosetta bin directory is on their PATH. Just setting this variable
  #to PATH does not work currently because of underlying functions in the package but
  #there should not be a reason this could be possible. Alternatively, we could have one
  #single "ROSETTAPATH" variable that just points to where the user installed rosetta's
  #top level directory. Then we can just use that and append the bin/ and database/ locations as needed.
  if {[info exists env(ROSETTAPATH)]} {
    variable rosettaPath $env(ROSETTAPATH)
  } else {
    variable rosettaPath ""
  }
 # variable defaultGapfindMol "top"
  variable DefaultNStruct 5000
  variable DefaultCluster 0
  variable DefaultNPerTask 1
  variable DefaultTestRun 0
  variable DefaultResStart 1
  variable DefaultAlignTemplate "all"
  variable DefaultInsertion "no"
  variable DefaultWorkDir "[pwd]/workdir"
  variable workdir $DefaultWorkDir
  variable DefaultScore -0.3
  variable DefaultCSFlag 1
  variable DefaultSidechain "no"

  variable settings
  set ::MODELMAKER::settings(username) ""
}

proc modelmaker { args } { return [eval ::MODELMAKER::modelmaker $args] }

proc ::MODELMAKER::modelmaker_usage { } {
  puts "Usage: modelmaker <command> \[args...\]"
  puts "Commands:"
  puts "  gapfind               -- find gaps in structure"
  puts "  full_length_model     -- generate a complete template model"
  puts "  abinitio              -- run Rosetta abinitio structure prediction"
  puts "  insertion             -- run Rosetta structure prediction for insertion folding"
  puts "  analyze               -- analyze results of Rosetta abinitio structure prediction"
  puts "  refine                -- refine a structure with Rosetta using a density" 
  return

}

proc ::MODELMAKER::modelmaker { args } {

  set nargs [llength $args]
  if { $nargs == 0 } {
    modelmaker_usage
    error ""
  }

  # parse command
  set command [lindex $args 0]
  set args [lreplace $args 0 0]

  if { $command == "gapfind" } {
    return [eval ::MODELMAKER::gapfind $args]
  } elseif { $command == "full_length_model" } {
    return [eval ::MODELMAKER::full_length_model $args]
  } elseif { $command == "abinitio" } {
    return [eval ::MODELMAKER::abinitio $args]
  } elseif { $command == "analyze" } {
    return [eval ::MODELMAKER::analyze $args]
  } elseif { $command == "insertion" } {
    return [eval ::MODELMAKER::insertion $args]
  } elseif { $command == "refine" } {
    return [eval ::MODELMAKER::refine $args]
  } else {
    modelmaker_usage
    error "Unrecognized command."
  }

  return

}


proc ::MODELMAKER::insertion_usage  { } {
  variable DefaultNStruct
  variable DefaultWorkDir
  puts "Usage: modelmaker insertion -model <full length template pdb> -fragfiles <list of fragment files> \
    -sel <list of atomselection texts with selections to fold> -fasta <fasta file> \
     ?options?"
  puts "Options:"
  puts "  -jobname    <name prefix for job> (default: taken from -model)> "
  puts "  -nstruct    <number of structures to predict> (default: $DefaultNStruct)> "
  puts "  -workdir    <working/project directory for job> (default: $DefaultWorkDir)>"
}


proc ::MODELMAKER::insertion { args } {
  variable rosettaEXE
  variable rosettadbpath
  variable DefaultNStruct
  variable rosettaPath
  variable DefaultWorkDir
 #These need to be changed in the underlying package to refer to the variable instead.
#e.g., $::MODELMAKER::rosettaDBpath
#instead of using these 'global' variables which gets confusing and dangerous.
 global rosettapath
 global rosettaDBpath
 global platform

#$::env(PATH)

  set rosettapath $rosettaPath
  set rosettaDBpath $rosettadbpath
  set platform $rosettaEXE

  set nargs [llength [lindex $args 0]]
  if {$nargs == 0} {
    insertion_usage
    error ""
  }

  foreach {name val} $args {
    switch -- $name {
      -jobname { set arg(jobname) $val }
      -model { set arg(model) $val }
      -sel { set arg(sel) $val }
      -fasta { set arg(fasta) $val }
      -fragfiles { set arg(fragfiles) $val }
      -nstruct { set arg(nstruct) $val }
      -workdir { set arg(workdir) $val }
    }
  }

  if { [info exists arg(model)] } {
    #set model $arg(model)
  #NOTE: Right now, because of RosettaVMD package, this needs to be the pdb name without the .pdb
  #extension. Need to change RosettaVMD to not require this.
    set model [string range $arg(model) 0 [expr [string last ".pdb" $arg(model)] - 1 ]]
  } else {
    error "A full model pdb file must be specified!"
  }

  if { [info exists arg(sel)] } {
    set sel [list $arg(sel)]
  } else {
    error "An atomselection text must be specified!"
  }

  if { [info exists arg(fragfiles)] } {
    set fragfiles $arg(fragfiles)
  } else {
    error "At least two fragment files must be specified!"
  }

  if { [info exists arg(fasta)] } {
  #NOTE: Right now, because of RosettaVMD package, this needs to be the fasta name without the .fasta
  #extension. Need to change RosettaVMD to not require this.
    set fasta [string range $arg(fasta) 0 [expr [string last ".fasta" $arg(fasta)] - 1 ]]
  #  set fasta $arg(fasta)
  } else {
    error "A fasta file must be specified!"
  }


  if { [info exists arg(nstruct)] } {
    set nstruct $arg(nstruct)
  } else {
    set nstruct $DefaultNStruct
  }

  if { [info exists arg(jobname)] } {
    set jobname $arg(jobname)
  } else {
    set jobname $model
  }

  if { [info exists arg(workdir)] } {
    set ::MODELMAKER::workdir $arg(workdir)
  } else {
    set ::MODELMAKER::workdir $DefaultWorkDir
  }

  if { [file exists $::MODELMAKER::workdir] } {
    puts "The working directory already exists!"
    exit 1
  }


  file mkdir $::MODELMAKER::workdir

  ## preparing files ##
  # folder with all the files needed for setup/run
  file mkdir $::MODELMAKER::workdir/setup-$jobname
  file mkdir $::MODELMAKER::workdir/run-$jobname

  file copy $model.pdb $::MODELMAKER::workdir/setup-$jobname
  file copy $fasta.fasta $::MODELMAKER::workdir/setup-$jobname
  foreach fragfile $fragfiles {
    puts $fragfile
    file copy $fragfile $::MODELMAKER::workdir/setup-$jobname
  }
  #start_rosetta_abinitio $jobname $model [list "$sel"] $anchor [list $fragfiles] $nstruct $cluster $npertask $testrun
  set currentPWD [pwd]

  #start_rosetta_insertion rpn11_insertion rpn11_yeast_23-306_complete [list "resid 138 to 157"] [list "rpn11_yeast_23-306_frag9" "rpn11_yeast_23-306_frag3"] [pwd]/input rpn11_yeast_23-306 $nstruct
  #cd $::MODELMAKER::workdir
  start_rosetta_insertion $jobname $model $sel $fragfiles $fasta $nstruct
  #cd $currentPWD
  set temp_mol [mol new $arg(model)]
  set temp_sel [atomselect $temp_mol all]
  set resstart [lindex [lsort -integer [$temp_sel get resid]] 0]

  foreach pdb [glob "$::MODELMAKER::workdir/run-$jobname/pdb_out/*"] {
    set full_mol [mol new $pdb]
    set full_sel [atomselect $full_mol all]
    renumber $full_sel $resstart
    $full_sel writepdb $pdb
    mol delete $full_mol
  }
}

proc ::MODELMAKER::refine_usage { } {
  variable DefaultNStruct
  variable DefaultCluster
  variable DefaultNPerTask
  variable DefaultTestRun
  variable DefaultWorkDir
  variable DefaultScore
  variable DefaultCSFlag
  variable DefaultSidechain
  puts "Usage: modelmaker refine -model <full length template pdb> \
    -sel <list of atomselection texts with selections to fold> -anchor <anchor residue for coordinate restraints> \
    -density <density file to refine against in .mrc format> -res <resolution of the density in Angstroms> \
     ?options?"
  puts "Options:"
  puts "  -jobname    <name prefix for job> (default: taken from -model)> "
  puts "  -workdir    <working/project directory for job> (default: $DefaultWorkDir)>"
  puts "  -nstruct    <number of structures to predict> (default: $DefaultNStruct)> "
  puts "  -bestN      <number of structures with best scores to save> (default: same as -nstruct)> "
  puts "  -score      <Rosetta density score; lower values indicate lower weight> (default: $DefaultScore)> "
  puts "  -csflag     <CartesianSample flag (0 or 1) Only used for non-sidechain refinement> (default: $DefaultCSFlag)> "
  puts "  -sidechain  <refine side chains? (yes or no) (default: $DefaultSidechain)> "
}

proc ::MODELMAKER::refine { args } {

  variable rosettaEXE
  variable rosettadbpath
  variable DefaultNStruct
  variable rosettaPath
  variable DefaultWorkDir
  variable DefaultScore
  variable DefaultCSFlag
  variable DefaultSidechain
 #These need to be changed in the underlying package to refer to the variable instead.
#e.g., $::MODELMAKER::rosettaDBpath
#instead of using these 'global' variables which gets confusing and dangerous.
 global rosettapath
 global rosettaDBpath
 global platform
 global packagePath

#$::env(PATH)
  set packagePath $::env(RosettaVMDDIR)
  set rosettapath $rosettaPath
  set rosettaDBpath $rosettadbpath
  set platform $rosettaEXE

  set nargs [llength [lindex $args 0]]
  if {$nargs == 0} {
    refine_usage
    error ""
  }

  foreach {name val} $args {
    switch -- $name {
      -jobname { set arg(jobname) $val }
      -model { set arg(model) $val }
      -sel { set arg(sel) $val }
      -anchor { set arg(anchor) $val }
      -nstruct { set arg(nstruct) $val }
      -csflag { set arg(csflag) $val }
      -density { set arg(density) $val }
      -res { set arg(res) $val }
      -score { set arg(score) $val }
      -bestN { set arg(bestN) $val }
      -workdir { set arg(workdir) $val }
      -sidechain { set arg(sidechain) $val }
    }
  }

  if { [info exists arg(model)] } {
    set model [string range $arg(model) 0 [expr [string last ".pdb" $arg(model)] - 1 ]]
    #set model $arg(model)
  } else {
    error "A full model pdb file must be specified!"
  }

  if { [info exists arg(sel)] } {
    set sel $arg(sel)
  } else {
    error "An atomselection text must be specified!"
  }

  if { [info exists arg(anchor)] } {
    set anchor $arg(anchor)
  } else {
    error "An anchor residue id must be specified!"
  }
  
  if { [info exists arg(density)] } {
    #set density $arg(density)
    set density [string range $arg(density) 0 [expr [string last ".dx" $arg(density)] - 1 ]]
  } else {
    error "A density file must be specified!"
  }
  
  if { [info exists arg(res)] } {
    set res $arg(res)
  } else {
    error "The resolution of the density must be specified!"
  }
  
  if { [info exists arg(csflag)] } {
    set csflag $arg(csflag)
  } else {
    set csflag $DefaultCSFlag
  }
  
  if { [info exists arg(sidechain)] } {
    set sidechain $arg(sidechain)
  } else {
    set sidechain $DefaultSidechain
  }
  
  if { [info exists arg(score)] } {
    set score $arg(score)
  } else {
    set score $DefaultScore
  }
  
  if { [info exists arg(nstruct)] } {
    set nstruct $arg(nstruct)
  } else {
    set nstruct $DefaultNStruct
  }
  
  if { [info exists arg(bestN)] } {
    set bestN $arg(bestN)
  } else {
    set bestN $nstruct
  }

  if { [info exists arg(jobname)] } {
    set jobname $arg(jobname)
  } else {
    set jobname $model
  }


  if { [info exists arg(workdir)] } {
    set ::MODELMAKER::workdir $arg(workdir)
  } else {
    set ::MODELMAKER::workdir $DefaultWorkDir 

  }

  if { [file exists $::MODELMAKER::workdir] } {
    puts "The working directory already exists!"
    exit 1
  }


  file mkdir $::MODELMAKER::workdir

  ## preparing files ##
  # folder with all the files needed for setup/run
  file mkdir $::MODELMAKER::workdir/setup-$jobname
  file mkdir $::MODELMAKER::workdir/run-$jobname

  file copy $model.pdb $::MODELMAKER::workdir/setup-$jobname
  
  if { $sidechain == "yes" } {
    start_rosetta_refine_sidechains_density $jobname $model $sel $anchor $density $res $score $bestN $nstruct
  } else {
    start_rosetta_refine $jobname $model $sel $anchor $csflag $density $res $score $bestN $nstruct
  }
}

proc ::MODELMAKER::abinitio_usage { } {
  variable DefaultNStruct
  variable DefaultCluster
  variable DefaultNPerTask
  variable DefaultTestRun
  variable DefaultWorkDir
  puts "Usage: modelmaker abinitio -model <full length template pdb> -fragfiles <list of fragment files> \
    -sel <list of atomselection texts with selections to fold> -anchor <anchor residue for coordinate restraints> \
     ?options?"
  puts "Options:"
  puts "  -jobname    <name prefix for job> (default: taken from -model)> "
  puts "  -workdir    <working/project directory for job> (default: $DefaultWorkDir)>"
  puts "  -nstruct    <number of structures to predict> (default: $DefaultNStruct)> "
  puts "  -testrun    <test run flag (0 or 1)> (default: $DefaultTestRun)> "
#hide these until wider functionality
#  puts "  -cluster    <run on cluster flag (0 or 1)> (default: $DefaultCluster)> "
#  puts "  -npertask   <tasks per job on cluster> (default: $DefaultNPerTask)> "
}

proc ::MODELMAKER::abinitio { args } {

  variable rosettaEXE
  variable rosettadbpath
  variable DefaultNStruct
  variable DefaultCluster
  variable DefaultNPerTask
  variable DefaultTestRun
  variable rosettaPath
  variable DefaultWorkDir
 #These need to be changed in the underlying package to refer to the variable instead.
#e.g., $::MODELMAKER::rosettaDBpath
#instead of using these 'global' variables which gets confusing and dangerous.
 global rosettapath
 global rosettaDBpath
 global platform

#$::env(PATH)

  set rosettapath $rosettaPath
  set rosettaDBpath $rosettadbpath
  set platform $rosettaEXE

  set nargs [llength [lindex $args 0]]
  if {$nargs == 0} {
    abinitio_usage
    error ""
  }

  foreach {name val} $args {
    switch -- $name {
      -jobname { set arg(jobname) $val }
      -model { set arg(model) $val }
      -sel { set arg(sel) $val }
      -anchor { set arg(anchor) $val }
      -fragfiles { set arg(fragfiles) $val }
      -nstruct { set arg(nstruct) $val }
      -cluster { set arg(cluster) $val }
      -npertask { set arg(npertask) $val }
      -testrun { set arg(testrun) $val }
      -workdir { set arg(workdir) $val }
    }
  }

  if { [info exists arg(model)] } {
    set model [string range $arg(model) 0 [expr [string last ".pdb" $arg(model)] - 1 ]]
    #set model $arg(model)
  } else {
    error "A full model pdb file must be specified!"
  }

  if { [info exists arg(sel)] } {
    set sel [list $arg(sel)]
  } else {
    error "An atomselection text must be specified!"
  }

  if { [info exists arg(anchor)] } {
    set anchor $arg(anchor)
  } else {
    error "An anchor residue id must be specified!"
  }

  if { [info exists arg(fragfiles)] } {
    set fragfiles $arg(fragfiles)
  } else {
    error "Two fragment files must be specified!"
  }

  if { [info exists arg(nstruct)] } {
    set nstruct $arg(nstruct)
  } else {
    set nstruct $DefaultNStruct
  }

  if { [info exists arg(cluster)] } {
    set cluster $arg(cluster)
  } else {
    set cluster $DefaultCluster
  }

  if { [info exists arg(npertask)] } {
    set npertask $arg(npertask)
  } else {
    set npertask $DefaultNPerTask
  }

  if { [info exists arg(testrun)] } {
    set testrun $arg(testrun)
  } else {
    set testrun $DefaultTestRun
  }

  if { [info exists arg(jobname)] } {
    set jobname $arg(jobname)
  } else {
    set jobname $model
  }


  if { [info exists arg(workdir)] } {
    set ::MODELMAKER::workdir $arg(workdir)
  } else {
    set ::MODELMAKER::workdir $DefaultWorkDir 

  }

  if { [file exists $::MODELMAKER::workdir] } {
    puts "The working directory already exists!"
    exit 1
  }


  file mkdir $::MODELMAKER::workdir

  ## preparing files ##
  # folder with all the files needed for setup/run
  file mkdir $::MODELMAKER::workdir/setup-$jobname
  file mkdir $::MODELMAKER::workdir/run-$jobname

  file copy $model.pdb $::MODELMAKER::workdir/setup-$jobname
  foreach fragfile [lindex $fragfiles 0] {
    file copy $fragfile $::MODELMAKER::workdir/setup-$jobname
  }
  #start_rosetta_abinitio $jobname $model [list "$sel"] $anchor [list $fragfiles] $nstruct $cluster $npertask $testrun
  start_rosetta_abinitio $jobname $model $sel $anchor $fragfiles $nstruct $cluster $npertask $testrun
}

proc ::MODELMAKER::analyze_usage { } {
  variable DefaultCluster
  variable DefaultAlignTemplate
  variable DefaultInsertion

  puts "Usage: modelmaker analyze -model <full length template pdb>  -nstruct <number of structures to analyze> \
    -comps <list of analysis tasks> ?options?"
  puts "Options:"
  puts "  -jobname    <name prefix for job> (default: taken from -model)> "
  puts "  -bestN      <best N structures to analyze> (default: same as -nstruct)> "
  puts "  -align_template      <selection text of template PDB for alignment> (default: $DefaultAlignTemplate)> "
  puts "  -align_rosetta       <selection text of predicted models for alignment> (default: taken from -align_template)> "
  puts "  -insertion   <analyzing output from insertion? 'yes' or 'no'> (default: $DefaultInsertion)"
#hide these until wider functionality
#  puts "  -cluster    <run on cluster flag (0 or 1)> (default: $DefaultCluster)> "

}

proc ::MODELMAKER::analyze { args } {
  variable rosettaEXE
  variable rosettadbpath
  variable DefaultCluster
  variable rosettaPath
  variable DefaultAlignTemplate
  variable DefaultInsertion

  set nargs [llength [lindex $args 0]]
  if {$nargs == 0} {
    analyze_usage
    error ""
  }

  foreach {name val} $args {
    switch -- $name {
      -jobname          { set arg(jobname) $val }
      -model            { set arg(model) $val }
      -template            { set arg(template) $val }
      -nstruct          { set arg(nstruct) $val }
      -cluster          { set arg(cluster) $val }
      -bestN            { set arg(bestN) $val }
      -align_template   { set arg(align_template) $val }
      -align_rosetta    { set arg(align_rosetta) $val }
      -comps            { set arg(comps) $val }
      -insertion        { set arg(insertion) $val }
    }
  }

  puts $arg(model)
  if { [info exists arg(model)] } {
  #NOTE: Right now, because of RosettaVMD package, this needs to be the pdb name without the .pdb
  #extension. Need to change RosettaVMD to not require this.
    set model [string range $arg(model) 0 [expr [string last ".pdb" $arg(model)] - 1 ]]
  } else {
    error "A full model pdb file must be specified!"
  }

  if { [info exists arg(template)] } {
  #NOTE: Right now, because of RosettaVMD package, this needs to be the pdb name without the .pdb
  #extension. Need to change RosettaVMD to not require this.
    # get the full path of the template file
    set dir [file dirname $arg(template)]
    if { $dir == "." } {
      set dir [pwd]
    }
    set template "$dir/$arg(template)"
  } else {
    error "A template pdb file must be specified!"
  }

  #can we determine this automatically?
  if { [info exists arg(nstruct)] } {
    set nstruct $arg(nstruct)
  } else {
    error "Number of structures to analyze must be specified!"
  }

#This is pretty cumbersome for the user and fairly opaque. Will have to find a succinct way of listing the
#possibilities and summarizing functionality in the usage command...
 #Also, like fragment files, is a list of lists so the user will have to do somehing like
#{{fragfile1 fragfile2}} if they want just one list....
  if { [info exists arg(comps)] } {
    set comps $arg(comps)
  } else {
    error "Analysis tasks must be specified!"
  }

  if { [info exists arg(align_template)] } {
    set align_template $arg(align_template)
  } else {
    set align_template "$DefaultAlignTemplate"
  }

  if { [info exists arg(align_rosetta)] } {
    set align_rosetta $arg(align_rosetta)
  } else {
    set align_rosetta $align_template
  }

  if { [info exists arg(bestN)] } {
    set bestN $arg(bestN)
  } else {
    set bestN $nstruct
  }

  if { [info exists arg(jobname)] } {
    set jobname $arg(jobname)
  } else {
    set jobname $model
  }

  if { [info exists arg(insertion)] } {
    set insertion $arg(insertion)
  } else {
    set insertion $DefaultInsertion
  }

  if { [info exists arg(cluster)] } {
    set cluster $arg(cluster)
  } else {
    set cluster $DefaultCluster
  }

  if { [info exists arg(workdir)] } {
    set ::MODELMAKER::workdir $arg(workdir)
  } else {
    set ::MODELMAKER::workdir [pwd]/workdir
  }

  if { ![file exists $::MODELMAKER::workdir] } {
    puts "The working directory does not exist!"
    exit 1
  }

 #These need to be changed in the underlying package to refer to the variable instead.
#e.g., $::MODELMAKER::rosettaDBpath
#instead of using these 'global' variables which gets confusing and dangerous.
 global rosettapath
 global rosettaDBpath
 global platform

 global packagePath
 global vmdexe
#can we avoid this or make it optional?
 global gnuplotexe

#$::env(PATH)

  set rosettapath $rosettaPath
  set rosettaDBpath $rosettadbpath
  set platform $rosettaEXE

  set packagePath $::env(RosettaVMDDIR)

#assumes in PATH which should be reasonable
  set vmdexe "vmd"
  set gnuplotexe "gnuplot"

  switch $insertion {
    "yes" {
      set modelname "$model\_S"
      set insert_model $model
    }
    "no"  {
      set modelname $model
      set insert_model ""
    }
  }
  #jobname mol bestN nstruct cluster align_template align_rosetta analysis_components
  analyze_abinitio $jobname $modelname $template $bestN $nstruct $cluster $align_template \
    $align_rosetta $comps {*}$insert_model
}

proc ::MODELMAKER::renumber { sel start } {
	if { [$sel num] == 0 } {
		puts "Error in renumber: empty selection!"
		return
	}
	set oresid [ $sel get resid ]
	set delta [ expr $start - [ lindex $oresid 0] ]
	set nresid { }
	foreach r $oresid {
		lappend nresid [ expr $r + $delta ]
	}
	$sel set resid $nresid
}

proc ::MODELMAKER::full_length_model_usage { } {
  variable DefaultResStart

  puts "Usage: modelmaker full_length_model -template <template pdb> -fragfiles <list of fragment files> \
    -fasta <fasta file> ?options?"
  puts "Options:"
 # puts "  -output <output file prefix (default: template name)> "
  puts "  -resstart <id of first residue> (default: $DefaultResStart)> "

}

proc ::MODELMAKER::full_length_model { args } {
  variable rosettaEXE
  variable DefaultResStart
  variable rosettaPath

  set nargs [llength [lindex $args 0]]
  if {$nargs == 0} {
    full_length_model_usage
    error ""
  }

  foreach {name val} $args {
    switch -- $name {
      -template { set arg(template) $val }
      #-name { set arg(name) $val }
      -fragfiles { set arg(fragfiles) $val }
      -fasta { set arg(fasta) $val }
      -resstart { set arg(resstart) $val }
      #-output { set arg(output) $val }
    }
  }

  if { [info exists arg(template)] } {
    set template $arg(template)
  } else {
    error "A template pdb file must be specified!"
  }

  if { [info exists arg(fragfiles)] } {
    set fragfiles $arg(fragfiles)
  } else {
    error "A list of two fragment files must be specified!"
  }

  if { [info exists arg(fasta)] } {
    set fasta $arg(fasta)
  } else {
    error "A fasta file must be specified!"
  }

  #if { [info exists arg(name)] } {
  #  set name $arg(name)
  #} else {
  #  error "A model name must be specified!"
  #}

  if { [info exists arg(resstart)] } {
    set resstart $arg(resstart)
  } else {
    set resstart $DefaultResStart
  }

#  if { [info exists arg(output)] } {
#    set output $arg(output)
#  } else {
#    #this mimics the default rosetta behaviour
#    set output $template
#  }



  #this assume only 9 and 3 length fragment files and in a specific order. Can
  #we implement some logic to look at the provided files and determine what all we have?
  exec $rosettaPath/full_length_model.$rosettaEXE -in:file:fasta $fasta \
    -loops:frag_files [lindex $fragfiles 0] [lindex $fragfiles 1] none \
    -loops:frag_sizes 9 3 1 \
    -in:file::s $template \
    -overwrite >> $template-full_length_model.log


  set full_mol [mol new ${template}_full_length.pdb]
  set full_sel [atomselect $full_mol all]
  renumber $full_sel $resstart
  $full_sel writepdb ${template}_full_length.pdb
  mol delete $full_mol

}





proc ::MODELMAKER::gapfind_usage { } {

  variable defaultGapfindSel
 # variable defaultGapfindMol

  puts "Usage: mdodelmaker gapfind ?options?"
  puts "Options:"
  puts "  -i <input pdb> "
  puts "  -sel <atom selection> (default: $defaultGapfindSel)"
  #puts "  -mol <molid> (find gaps in already loaded molecule) (default: $defaultGapfindMol)"
  puts "  -mol <molid> (find gaps in already loaded molecule)"

}



proc ::MODELMAKER::gapfind { args } {

  variable defaultGapfindSel
 # variable defaultGapfindMol

  set nargs [llength [lindex $args 0]]
  if {$nargs == 0} {
    gapfind_usage
    error ""
  }

#  set MOL coot_corr_of_9_refine_4_fin_fixed

#  mol delete all
#  mol new $MOL.pdb

  foreach {name val} $args {
    switch -- $name {
      -i { set arg(i) $val }
      -sel { set arg(sel) $val }
      -mol { set arg(mol) $val }
    }
  }

  if { [info exists arg(i)] && [info exists arg(mol)] } {
    error "only an input pdb OR a mol id can be specified at one time"
  } elseif { ![info exists arg(i)] && ![info exists arg(mol)] } {
    error "either an input pdb OR a mol id must be specified"
  }

  if { [info exists arg(i)] } {
    set inputpdb $arg(i)
  } else {
    set inputpdb ""
  }

  if { [info exists arg(sel)] } {
    set inputsel $arg(sel)
  } else {
    set inputsel $defaultGapfindSel
  }

  if { [info exists arg(mol)] } {
    set inputmol $arg(mol)
  } else {
  #  set inputmol $defaultGapfindMol
    set inputmol ""
  }

  if { $inputpdb != ""} {
    set MOLID [mol new $inputpdb]
  } else {
    set MOLID [molinfo $inputmol get id]
  }

  set molname [molinfo $MOLID get name]

  set chains [lsort -unique [[atomselect $MOLID all] get chain]]
#set generalOut [open "$MOL-missing.html" w]
  foreach chain $chains {
    puts $chain
    set counter 0
    set sel [atomselect $MOLID "chain $chain"]
    set residues [lsort -integer -unique [$sel get resid]]
    set missing []
    #puts $residues
    #set first [lindex $residues 0]
    set first 1
    set last [lindex $residues end]
    set delta [expr $last - $first]
    #puts "$first $last $delta"
    for {set i $first} {$i <= $delta} {incr i} {
      #puts $i
      set d [expr $i - $first - $counter]
      if {$i != [lindex $residues $d]} {
        #puts "$i"
        lappend missing $i
        incr counter
      }
    }
    if {1} {
    set f [open "$molname-$chain-missing.txt" w]
    foreach var $missing {
      puts $f "$var"
    }
    close $f
    }
    #puts $generalOut "<head>"
    #puts $generalOut "<link rel=\"stylesheet\" type=\"text/css\" href=\"style.css\">"
    #puts $generalOut "</head>"
    #puts $generalOut "<h1> chain $chain </h1>"
    #puts $generalOut [html::tableFromList $missing]
  }

}

#proc ::MODELMAKER::decr { int { n 1 } } {
#    if { [ catch {
#        uplevel incr $int -$n
#    } err ] } {
#        return -code error "decr: $err"
#    }
#    return [ uplevel set $int ]
#}
