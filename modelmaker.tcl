package require RosettaVMD
package require MakePsf
package require ssrestraints
package require cispeptide
package require chirality
package require mdff
package require SSAnalysis

package provide modelmaker 0.1

namespace eval ::MODELMAKER {
  variable defaultGapfindSel "protein"
  if {[info exists env(ROSETTAEXE)]} {
    variable rosettaEXE $env(ROSETTAEXE)
  } else {
    switch $tcl_platform(os) {
      "Darwin" {
        variable rosettaEXE "static.macos*release"
      }
      "Linux" {
        variable rosettaEXE "static.linux*release"
      }
      default {
        variable rosettaEXE "Unrecognized"
      }
    }
  }
  if {[info exists env(ROSETTADB)]} {
    variable rosettadbpath $env(ROSETTADB)
  } else {
    variable rosettadbpath ""
  }

  if {[info exists env(ROSETTAPATH)]} {
    variable rosettaPath $env(ROSETTAPATH)
  } else {
    variable rosettaPath ""
  }
  variable DefaultNStruct 5000
  variable DefaultCluster 0
  variable DefaultNPerTask 1
  variable DefaultTestRun 0
  variable DefaultResStart 1
  variable DefaultAlignTemplate "all"
  variable DefaultInsertion "no"
  variable DefaultKmin 2
  variable DefaultKmax 10
  variable DefaultWorkDir "[pwd]/workdir"
  variable DefaultMDFFWorkDir "[pwd]/mdff-workdir"
  variable workdir $DefaultWorkDir
  variable DefaultScore -0.3
  variable DefaultRefineMode "backbone"
  variable DefaultPDB2SeqSel "protein" 
  variable DefaultFixed "none"
  variable DefaultBestN 1
  variable DefaultNumSteps 10000
  variable DefaultMinSteps 200
  variable DefaultGScale 0.3
  variable DefaultGridPDBSel "noh and (protein or nucleic)"
  variable DefaultDCDFreq 100
  variable DefaultTopFiles [list [file join $env(CHARMMTOPDIR) top_all36_prot.rtf] \
    [file join $env(CHARMMTOPDIR) top_all36_lipid.rtf] \
    [file join $env(CHARMMTOPDIR) top_all36_na.rtf] \
    [file join $env(CHARMMTOPDIR) top_all36_carb.rtf] \
    [file join $env(CHARMMTOPDIR) top_all36_cgenff.rtf] \
    [file join $env(CHARMMTOPDIR) toppar_all36_carb_glycopeptide.str] \
    [file join $env(CHARMMTOPDIR) toppar_water_ions_namd.str] ]

  variable DefaultParFiles [list [file join $env(CHARMMPARDIR) par_all36_prot.prm]\
    [file join $env(CHARMMPARDIR) par_all36_lipid.prm] \
  [file join $env(CHARMMPARDIR) par_all36_na.prm] [file join $env(CHARMMPARDIR) par_all36_carb.prm] \
  [file join $env(CHARMMPARDIR) par_all36_cgenff.prm] [file join $env(CHARMMPARDIR) toppar_all36_carb_glycopeptide.str]\
  [file join $env(CHARMMPARDIR) toppar_water_ions_namd.str]]
 
  variable DefaultThreshold 0.1 
  variable DefaultMaskCutoff 2 
  
  variable MPINP
  variable settings
  set ::MODELMAKER::settings(username) ""
  
  variable DefaultModelSel "all"
  variable DefaultModelN 1
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
  puts "  pdb2seq               -- get the single-letter amino acid sequence" 
  puts "  seqsub                -- extract subrange of fasta sequence" 
  puts "  makepsf               -- make a .psf file from a pdb" 
  puts "  mdff                  -- run an MDFF simulation" 
  puts "  cccolor               -- calculate local cross correlations"
  puts "  ssanalysis            -- analyze the secondary structure of the generated structures" 
  puts "  get_empty_density     -- find the unassigned empty density around a structure" 
  puts "  model                 -- generate a homology model with MODELLER" 
  puts "  fragments             -- generate fragment files" 
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
  } elseif { $command == "pdb2seq" } {
    return [eval ::MODELMAKER::pdb2seq $args]
  } elseif { $command == "seqsub" } {
    return [eval ::MODELMAKER::seqsub $args]
  } elseif { $command == "makepsf" } {
    return [eval ::MODELMAKER::makepsf $args]
  } elseif { $command == "mdff" } {
    return [eval ::MODELMAKER::quick_mdff $args]
  } elseif { $command == "cccolor" } {
    return [eval ::MODELMAKER::cccolor $args]
  } elseif { $command == "get_empty_density" } {
    return [eval ::MODELMAKER::get_empty_density $args]
  } elseif { $command == "ssanalysis" } {
    return [eval ::MODELMAKER::ssanalysis $args]
  } elseif { $command == "model" } {
    return [eval ::MODELMAKER::model $args]
  } elseif { $command == "fragments" } {
    return [eval ::MODELMAKER::fragments $args]
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
  puts "  -np         <Number of processors to use. MPI version only> "
}


proc ::MODELMAKER::insertion { args } {
  variable rosettaEXE
  variable rosettadbpath
  variable DefaultNStruct
  variable rosettaPath
  variable DefaultWorkDir
  variable MPINP
  global rosettapath
  global rosettaDBpath
  global platform

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
      -np { set arg(np) $val }
      default { puts "Unknown argument $name"; return }
    }
  }

  if { [info exists arg(model)] } {
    set model [file rootname [file tail $arg(model)]]
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
  
  if { [info exists arg(np)] } {
    set MPINP $arg(np)
  } elseif { [string match "*mpi.*" $rosettaEXE] }  {
    error "number of processors (-np) must be specified when using MPI version of Rosetta!"
  }

  if { [info exists arg(workdir)] } {
    set ::MODELMAKER::workdir $arg(workdir)
  } else {
    set ::MODELMAKER::workdir $DefaultWorkDir
  }

#  if { [file exists $::MODELMAKER::workdir] } {
#    puts "The working directory already exists!"
#    exit 1
#  }


  if {![file exists $::MODELMAKER::workdir]} { file mkdir $::MODELMAKER::workdir }

  ## preparing files ##
  # folder with all the files needed for setup/run
  file mkdir $::MODELMAKER::workdir/setup-$jobname
  file mkdir $::MODELMAKER::workdir/run-$jobname

  set tmpmol [mol new $arg(model)]
  set tmpsel [atomselect $tmpmol "noh"]
  $tmpsel writepdb "$::MODELMAKER::workdir/setup-$jobname/[file tail $arg(model)]"
  
  file copy -force $fasta.fasta $::MODELMAKER::workdir/setup-$jobname
  foreach fragfile $fragfiles {
    puts $fragfile
    file copy -force $fragfile $::MODELMAKER::workdir/setup-$jobname
  }
  set currentPWD [pwd]

  start_rosetta_insertion $jobname $model $sel $fragfiles $fasta $nstruct
  
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
  variable DefaultRefineMode
  puts "Usage: modelmaker refine -model <full length template pdb> \
    -sel <list of atomselection texts with selections to fold> -anchor <anchor residue for coordinate restraints> \
    -density <density file to refine against in .mrc format> -res <resolution of the density in Angstroms> \
     ?options?"
  puts "Options:"
  puts "  -mode       <refinement mode (backbone, sidechain, or cartesian> (default: $DefaultRefineMode)"
  puts "  -jobname    <name prefix for job> (default: taken from -model) "
  puts "  -workdir    <working/project directory for job> (default: $DefaultWorkDir)"
  puts "  -nstruct    <number of structures to predict> (default: $DefaultNStruct) "
  puts "  -bestN      <number of structures with best scores to save> (default: same as -nstruct) "
  puts "  -score      <Rosetta density score; lower values indicate lower weight> (default: $DefaultScore) "
  puts "  -np         <Number of processors to use. MPI version only> "
}

proc ::MODELMAKER::refine { args } {

  variable rosettaEXE
  variable rosettadbpath
  variable DefaultNStruct
  variable rosettaPath
  variable DefaultWorkDir
  variable DefaultScore
  variable MPINP
  variable DefaultRefineMode
  global rosettapath
  global rosettaDBpath
  global platform
  global packagePath

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
      -density { set arg(density) $val }
      -res { set arg(res) $val }
      -score { set arg(score) $val }
      -bestN { set arg(bestN) $val }
      -workdir { set arg(workdir) $val }
      -mode { set arg(mode) $val }
      -np { set arg(np) $val }
      default { puts "Unknown argument $name"; return  }
    }
  }

  if { [info exists arg(model)] } {
    set model [file rootname [file tail $arg(model)]]
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
    set density [string range $arg(density) 0 [expr [string last ".mrc" $arg(density)] - 1 ]]
  } else {
    error "A density file must be specified!"
  }
  
  if { [info exists arg(res)] } {
    set res $arg(res)
  } else {
    error "The resolution of the density must be specified!"
  }
  
  if { [info exists arg(mode)] } {
    set mode $arg(mode)
  } else {
    set mode $DefaultRefineMode
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
  
  if { [info exists arg(np)] } {
    set MPINP $arg(np)
  } elseif { [string match "*mpi.*" $rosettaEXE] }  {
    error "number of processors (-np) must be specified when using MPI version of Rosetta!"
  }


  if { [info exists arg(workdir)] } {
    set ::MODELMAKER::workdir $arg(workdir)
  } else {
    set ::MODELMAKER::workdir $DefaultWorkDir 

  }

#  if { [file exists $::MODELMAKER::workdir] } {
#    puts "The working directory already exists!"
#    exit 1
#  }


  if {![file exists $::MODELMAKER::workdir]} { file mkdir $::MODELMAKER::workdir }

  ## preparing files ##
  # folder with all the files needed for setup/run
  file mkdir $::MODELMAKER::workdir/setup-$jobname
  file mkdir $::MODELMAKER::workdir/run-$jobname

  set tmpmol [mol new $arg(model)]
  set tmpsel [atomselect $tmpmol "noh"]
  $tmpsel writepdb "$::MODELMAKER::workdir/setup-$jobname/[file tail $arg(model)]"
    
  if { $mode == "cartesian" } {
    set csflag 1
  } else {
    set csflag 0
  }

  if { $mode == "sidechain" } {
    start_rosetta_refine_sidechains_density $jobname $model $sel $anchor $density $res $score $bestN $nstruct
  } else {
    start_rosetta_refine $jobname $model $sel $anchor $csflag $density $res $score $bestN $nstruct
  }
  
  set inmol [mol new $arg(model)]
  foreach pdb [glob -nocomplain $::MODELMAKER::workdir/run-$jobname/*_best*.pdb] {
    set outmol [mol new $pdb]
    regen_segnames $inmol $outmol
    set sel [atomselect $outmol all]
    $sel writepdb $pdb
    $sel delete  
  }

}

proc ::MODELMAKER::abinitio_usage { } {
  variable DefaultNStruct
  variable DefaultCluster
  variable DefaultNPerTask
  variable DefaultTestRun
  variable DefaultWorkDir
  puts "Usage: modelmaker abinitio -model <full length template pdb> -fragfiles <list of fragment files> \
    -sel <list of atomselection texts with selections to fold> -anchor <atomselection text for anchor residue for coordinate restraints> \
     ?options?"
  puts "Options:"
  puts "  -jobname    <name prefix for job> (default: taken from -model)> "
  puts "  -workdir    <working/project directory for job> (default: $DefaultWorkDir)>"
  puts "  -nstruct    <number of structures to predict> (default: $DefaultNStruct)> "
  puts "  -testrun    <test run flag (0 or 1)> (default: $DefaultTestRun)> "
  puts "  -np         <Number of processors to use. MPI version only> "
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
  variable MPINP
  global rosettapath
  global rosettaDBpath
  global platform


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
      -np { set arg(np) $val }
      default { puts "Unknown argument $name"; return  }
    }
  }

  if { [info exists arg(model)] } {
    set model [file rootname [file tail $arg(model)]]
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
  
  if { [info exists arg(np)] } {
    set MPINP $arg(np)
  } elseif { [string match "*mpi.*" $rosettaEXE] }  {
    error "number of processors (-np) must be specified when using MPI version of Rosetta!"
  }


  if { [info exists arg(workdir)] } {
    set ::MODELMAKER::workdir $arg(workdir)
  } else {
    set ::MODELMAKER::workdir $DefaultWorkDir 

  }

#  if { [file exists $::MODELMAKER::workdir] } {
#    puts "The working directory already exists!"
#    exit 1
#  }


  if {![file exists $::MODELMAKER::workdir]} { file mkdir $::MODELMAKER::workdir }

  ## preparing files ##
  # folder with all the files needed for setup/run
  file mkdir $::MODELMAKER::workdir/setup-$jobname
  file mkdir $::MODELMAKER::workdir/run-$jobname

  set tmpmol [mol new $arg(model)]
  set tmpsel [atomselect $tmpmol "noh"]
  $tmpsel writepdb "$::MODELMAKER::workdir/setup-$jobname/[file tail $arg(model)]"
  
  foreach fragfile [lindex $fragfiles 0] {
    file copy -force $fragfile $::MODELMAKER::workdir/setup-$jobname
  }
  start_rosetta_abinitio $jobname $model $sel $anchor $fragfiles $nstruct $cluster $npertask $testrun
}

proc ::MODELMAKER::analyze_usage { } {
  variable DefaultCluster
  variable DefaultAlignTemplate
  variable DefaultInsertion
  variable DefaultKmin
  variable DefaultKmax

  puts "Usage: modelmaker analyze -model <full length template pdb>  -nstruct <number of structures to analyze> \
    -comps <list of analysis tasks> ?options?"
  puts "Options:"
  puts "  -jobname    <name prefix for job> (default: taken from -model)> "
  puts "  -bestN      <best N structures to analyze> (default: same as -nstruct)> "
  puts "  -align_template      <selection text of template PDB for alignment> (default: $DefaultAlignTemplate)> "
  puts "  -align_rosetta       <selection text of predicted models for alignment> (default: taken from -align_template)> "
  puts "  -insertion   <analyzing output from insertion? 'yes' or 'no'> (default: $DefaultInsertion)"
#change these to more user-friendly designations and explain caveats
  puts "  -cluster    <clustering mode to use. 0 for Rosetta, 1 for Python> (default: $DefaultCluster)> "
  puts "  -kmin       <Min number of clusters for silhouette calc, cluster mode 1 only> (default: $DefaultKmin)> "
  puts "  -kmax       <Max number of clusters for silhouette calc, cluster mode 1 only> (default: $DefaultKmax)> "

}

proc ::MODELMAKER::analyze { args } {
  variable rosettaEXE
  variable rosettadbpath
  variable DefaultCluster
  variable rosettaPath
  variable DefaultAlignTemplate
  variable DefaultInsertion
  variable DefaultKmin
  variable DefaultKmax

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
      -kmin             { set arg(kmin) $val }
      -kmax        { set arg(kmax) $val }
      default { puts "Unknown argument $name"; return  }
    }
  }

  puts $arg(model)
  if { [info exists arg(model)] } {
    set model [string range $arg(model) 0 [expr [string last ".pdb" $arg(model)] - 1 ]]
  } else {
    error "A full model pdb file must be specified!"
  }

  if { [info exists arg(template)] } {
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
  
  if { [info exists arg(kmin)] } {
    set kmin $arg(kmin)
  } else {
    set kmin $DefaultKmin
  }
  
  if { [info exists arg(kmax)] } {
    set kmax $arg(kmax)
  } else {
    set kmax $DefaultKmax
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

  global rosettapath
  global rosettaDBpath
  global platform

  global packagePath
  global vmdexe
#can we avoid this or make it optional?
  global gnuplotexe


  set rosettapath $rosettaPath
  set rosettaDBpath $rosettadbpath
  set platform $rosettaEXE

  set packagePath $::env(RosettaVMDDIR)

#assumes in PATH which should be reasonable requirement
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
  analyze_abinitio $jobname $modelname $template $bestN $nstruct $cluster $align_template \
    $align_rosetta $comps $kmin $kmax {*}$insert_model
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
      -fragfiles { set arg(fragfiles) $val }
      -fasta { set arg(fasta) $val }
      -resstart { set arg(resstart) $val }
      default { puts "Unknown argument $name"; return  }
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


  if { [info exists arg(resstart)] } {
    set resstart $arg(resstart)
  } else {
    set resstart $DefaultResStart
  }

  set name [file rootname [file tail $template]]
  set mol [mol new $template]
  set sel [atomselect $mol "protein and noh"] 
  format_pdb $sel "$name.formatted.pdb" 
  #this assume only 9 and 3 length fragment files and in a specific order. Can
  #we implement some logic to look at the provided files and determine what all we have?
  exec [glob $rosettaPath/full_length_model.$rosettaEXE] -in:file:fasta $fasta \
    -loops:frag_files [lindex $fragfiles 0] [lindex $fragfiles 1] none \
    -loops:frag_sizes 9 3 1 \
    -in:file::s "$name.formatted.pdb" \
    -overwrite >> $name-full_length_model.log
  
  file delete "$name.formatted.pdb"
  set full_mol [mol new ${name}.formatted.pdb_full_length.pdb]
 
  set chain [lindex [$sel get chain] 0]
  set segname [lindex [$sel get segname] 0]
  set full_sel [atomselect $full_mol all]
  if {$chain != ""} { $full_sel set chain $chain } 
  if {$segname != ""} { $full_sel set segname $segname }
  renumber $full_sel $resstart
  
  $full_sel writepdb ${name}_full_length.pdb
  mol delete $full_mol
  file delete ${name}.formatted.pdb_full_length.pdb

}





proc ::MODELMAKER::gapfind_usage { } {
  variable defaultGapfindSel

  puts "Usage: mdodelmaker gapfind ?options?"
  puts "Options:"
  puts "  -i <input pdb> "
  puts "  -sel <atom selection text> (default: $defaultGapfindSel)"
  puts "  -mol <molid> (find gaps in already loaded molecule)"

}



proc ::MODELMAKER::gapfind { args } {

  variable defaultGapfindSel

  set nargs [llength [lindex $args 0]]
  if {$nargs == 0} {
    gapfind_usage
    error ""
  }

  foreach {name val} $args {
    switch -- $name {
      -i { set arg(i) $val }
      -sel { set arg(sel) $val }
      -mol { set arg(mol) $val }
      default { puts "Unknown argument $name"; return  }
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
    set inputmol ""
  }

  if { $inputpdb != ""} {
    set MOLID [mol new $inputpdb]
  } else {
    set MOLID $inputmol ;
  }

  set molname [molinfo $MOLID get name]

  set chains [lsort -unique [[atomselect $MOLID $inputsel] get chain]]
  foreach chain $chains {
    puts $chain
    set counter 0
    set sel [atomselect $MOLID "chain $chain"]
    set residues [lsort -integer -unique [$sel get resid]]
    set missing []
    #set first [lindex $residues 0]
    set first 1
    set last [lindex $residues end]
    set delta [expr $last - $first]
    for {set i $first} {$i <= $delta} {incr i} {
      set d [expr $i - $first - $counter]
      if {$i != [lindex $residues $d]} {
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
  }

}

proc ::MODELMAKER::format_pdb {sel filename} {
  $sel set occupancy 1
  $sel set beta 1
  $sel writepdb $filename

  set frpdb [open $filename "r"]
  set spdb [read $frpdb]
  close $frpdb
  set fwpdb [open $filename "w"]
  regsub -all "HSD" $spdb "HIS" spdb
  regsub -all "HSE" $spdb "HIS" spdb
  regsub -all "URA" $spdb "  U" spdb
  regsub -all "ADE" $spdb "  A" spdb
  regsub -all "CYT" $spdb "  C" spdb
  regsub -all "GUA" $spdb "  G" spdb
  regsub -all "THY" $spdb "  T" spdb
  regsub -all "CYN" $spdb "CYS" spdb
  regsub -all -line {^.*CRYST.*$} $spdb " " spdb
  puts $fwpdb $spdb
  close $fwpdb
}

proc ::MODELMAKER::pdb2seq_usage { } {
  variable DefaultPDB2SeqSel 
  puts "Usage: mdodelmaker pdb2seq ?options?"
  puts "Options:"
  puts "  -i <input pdb> "
  puts "  -sel <atom selection text> (default: $DefaultPDB2SeqSel)"
  puts "  -mol <molid> (get sequence of an already loaded molecule)"
  puts "  -o   <filename> (output sequence to a file instead)"

}

proc ::MODELMAKER::pdb2seq {args} {
  variable DefaultPDB2SeqSel 
  
  set nargs [llength [lindex $args 0]]
  if {$nargs == 0} {
    pdb2seq_usage
    error ""
  }
  
  foreach {name val} $args {
    switch -- $name {
      -i { set arg(i) $val }
      -o { set arg(o) $val }
      -mol { set arg(mol) $val }
      -sel { set arg(sel) $val }
      default { puts "Unknown argument $name"; return  }
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
    set sel $arg(sel)
  } else {
    set sel $DefaultPDB2SeqSel
  }

  if { [info exists arg(mol)] } {
    set inputmol $arg(mol)
  } else {
    set inputmol ""
  }
  
  if { [info exists arg(o)] } {
    set output $arg(o)
  } else {
    set output ""
  }

  if { $inputpdb != ""} {
    set MOLID [mol new $inputpdb]
  } else {
    set MOLID $inputmol
  }
  

  set seqsel [atomselect $MOLID $sel]
 
  set sequence ""
  foreach resid [lsort -integer -unique [$seqsel get resid]] {
    set ressel [atomselect $MOLID "resid $resid"]
    set resname [lindex [$ressel get resname] 0]  
    switch $resname {
      GLY { set res "G" }
      ALA { set res "A" }
      LEU { set res "L" }
      MET { set res "M" }
      PHE { set res "F" }
      TRP { set res "W" }
      LYS { set res "K" }
      GLN { set res "Q" }
      GLU { set res "E" }
      SER { set res "S" }
      PRO { set res "P" }
      VAL { set res "V" }
      ILE { set res "I" }
      CYS { set res "C" }
      CYN { set res "C" }
      TYR { set res "Y" }
      HIS { set res "H" }
      HSE { set res "H" }
      HSD { set res "H" }
      ARG { set res "R" }
      ASN { set res "N" }
      ASP { set res "D" }
      THR { set res "T" }
      default { set res "X" }
    }
    
    lappend sequence $res
    $ressel delete
  }
  
  if {$output != ""} {
   set f [open $output w] 
   foreach seq $sequence {
    puts -nonewline $f $seq
   }
   close $f
  } else {
    return $sequence
  }
}

proc ::MODELMAKER::seqsub_usage { } {
  puts "Usage: mdodelmaker seqsub -i <input fasta> -o <filename> (output file name) \
    -start <residue number> (starting residue number) \
    -end <residue number> (ending residue number) "

}

proc ::MODELMAKER::seqsub { args } {
  set nargs [llength [lindex $args 0]]
  if {$nargs == 0} {
    seqsub_usage
    error ""
  }
  
  foreach {name val} $args {
    switch -- $name {
      -i { set arg(i) $val }
      -o { set arg(o) $val }
      -start { set arg(start) $val }
      -end { set arg(end) $val }
      default { puts "Unknown argument $name"; return  }
    }
  }

  if { [info exists arg(i)] } {
    set inputfasta $arg(i)
  } else {
    error "input fasta file required"
  }
  
  if { [info exists arg(o)] } {
    set output $arg(o)
  } else {
    error "output file name required"
  }
  
  if { [info exists arg(start)] } {
    set start $arg(start)
  } else {
    error "starting residue number required"
  }
  
  if { [info exists arg(end)] } {
    set end $arg(end)
  } else {
    error "ending residue number required"
  }

  set infile [open $inputfasta r]
  while { [gets $infile line] >=0 } {
    if { [string range $line 0 0] != ">"} {
     append sequence $line 
    }
  }
  close $infile
  set outfile [open $output w]
  puts $outfile [string range $sequence [expr $start - 1] [expr $end - 1] ]
  close $outfile
}

proc ::MODELMAKER::makepsf_usage {} {
  variable DefaultTopFiles
  puts "Usage: mdodelmaker makepsf -pdb <input pdb file> ?options?"
  puts "Options:"
  puts "  -topfiles <list of topology files to use>(Default: $DefaultTopFiles) "
  puts "  -chseg    <file containing lines of desired 'name' 'chain' 'segname'> (Default: from input pdb)"
  puts "  -prot     <list of resids and their protonation state, e.g, 279 HSE> (Default: none)"
}

proc ::MODELMAKER::makepsf { args } {
  variable DefaultTopFiles

  set nargs [llength [lindex $args 0]]
  if {$nargs == 0} {
    makepsf_usage
    error ""
  }
  
  foreach {name val} $args {
    switch -- $name {
      -pdb { set arg(pdb) $val }
      -topfiles { set arg(topfiles) $val }
      -chseg { set arg(chseg) $val }
      -prot { set arg(prot) $val }
      default { puts "Unknown argument $name"; return  }
    }
  }
  
  if { [info exists arg(pdb)] } {
    set pdb $arg(pdb)
  } else {
    error "input pdb file required"
  }
  
  if { [info exists arg(topfiles)] } {
    set topfiles $arg(topfiles)
  } else {
    set topfiles $DefaultTopFiles
  }
  
  if { [info exists arg(prot)] } {
    set prot $arg(prot)
  } else {
    set prot ""
  }
  
  if { [info exists arg(chseg)] } {
  
   ### input file structure: name chain segname 
    set inputfile [open $arg(chseg) r]
    set data [read $inputfile]
    close $inputfile
    set chseg [split $data "\n"]
  
  } else {
    set mol [mol new $pdb]
    set sel [atomselect $mol "all"]
    set chains [$sel get chain]
    set segnames [$sel get segname]
    #get unique list
    foreach chain $chains {dict set tmp $chain 1}
    set chains [dict keys $tmp]
    foreach segname $segnames {dict set tmp2 $segname 1}
    set segnames [dict keys $tmp2]

    foreach chain $chains segname $segnames {
      lappend chseg  [list $pdb $chain $segname]
    }
  
  }
  auto_makepsf $pdb $topfiles $chseg $prot
}

proc ::MODELMAKER::quick_mdff_usage { } {
  variable DefaultFixed
  variable DefaultBestN
  variable DefaultNumSteps
  variable DefaultMinSteps
  variable DefaultGScale
  variable DefaultTopFiles
  variable DefaultParFiles
  variable DefaultDCDFreq
  variable DefaultMDFFWorkDir
  variable DefaultGridPDBSel
  
  puts "Usage: mdodelmaker mdff -pdb <input pdb file> -density <input density file> \
    -res <resolution of density in Angstroms> ?options?"
  puts "Options:"
  puts "  -gridpdbsel <atom selections for MDFF> (default: $DefaultGridPDBSel)"
  puts "  -workdir    <working/project directory for job> (default: $DefaultMDFFWorkDir)"
  puts "  -jobname    <name prefix for job> (default: taken from -model) "
  puts "  -fixed      <atomselect text for fixed atoms> (default: $DefaultFixed) "
  puts "  -gscale     <grid scaling factor for fitting forces> (default: $DefaultGScale) "
  puts "  -minsteps   <number of minimization steps> (default: $DefaultMinSteps) "
  puts "  -numsteps   <number of simulation steps> (default: $DefaultNumSteps) "
  puts "  -chseg      <file containing lines of desired 'name' 'chain' 'segname'> (Default: from input pdb)"
  puts "  -topfiles   <list of topology files to use>(Default: $DefaultTopFiles) "
  puts "  -parfiles   <list of parameter files to use>(Default: $DefaultParFiles) "
  puts "  -dcdfreq    <frequencey of dcd output>(Default: $DefaultDCDFreq) "
  puts "  -namdargs   <arguments to pass to namd> (default: none)"
  puts "  -ss         <secondary structure restraint file> (default: auto-generated)"
  puts "  -cis        <cispeptide restraint file> (default: auto-generated)"
  puts "  -chir       <chirality restraint file> (default: auto-generated)"
  puts "  -extrab     <list of additional extrabonds files> (default: none)"

}

proc ::MODELMAKER::quick_mdff { args } {
  variable DefaultFixed
  variable DefaultBestN
  variable DefaultNumSteps
  variable DefaultMinSteps
  variable DefaultGScale
  variable DefaultTopFiles
  variable DefaultParFiles
  variable DefaultDCDFreq
  variable DefaultMDFFWorkDir
  variable DefaultGridPDBSel
  
  set nargs [llength [lindex $args 0]]
  if {$nargs == 0} {
    quick_mdff_usage
    error ""
  }
  
  foreach {name val} $args {
    switch -- $name {
      -workdir { set arg(workdir) $val }
      -jobname { set arg(jobname) $val }
      -pdb { set arg(pdb) $val }
      -density { set arg(density) $val }
      -gridpdbsel { set arg(gridpdbsel) $val }
      -fixed { set arg(fixed) $val }
      -gscale { set arg(gscale) $val }
      -density { set arg(density) $val }
      -minsteps { set arg(minsteps) $val }
      -numsteps { set arg(numsteps) $val }
      -res { set arg(res) $val }
      -chseg { set arg(chseg) $val }
      -topfiles { set arg(topfiles) $val }
      -parfiles { set arg(parfiles) $val }
      -dcdfreq { set arg(dcdfreq) $val }
      -namdargs { set arg(namdargs) $val }
      -ss { set arg(ss) $val }
      -cis { set arg(cis) $val }
      -chir { set arg(chir) $val }
      -extrab { set arg(extrab) $val }
      default { puts "Unknown argument $name"; return  }
    }
  }
  
  if { [info exists arg(pdb)] } {
    set pdb $arg(pdb)
  } else {
    error "input pdb file required!"
  }
  
  
  if { [info exists arg(density)] } {    
    set density $arg(density)
  } else {
    error "A density file must be specified!"
  }
  
  if { [info exists arg(res)] } {
    set res $arg(res)
  } else {
    error "resolution of density must be specified!"
  }
  
  if { [info exists arg(gridpdbsel)] } {
    #if more than one map, we already should have a list
    if {[llength $density] > 1} {
      set gridpdbsel $arg(gridpdbsel)
    #we have one map, so need to make this a list so the entire first entry is the atom selection text  
    } else {
      lappend gridpdbsel $arg(gridpdbsel)
    }
  } else {
    foreach den $density {
      lappend gridpdbsel $DefaultGridPDBSel
    }
  }
  
  if { [info exists arg(workdir)] } {
    set workdir $arg(workdir)
  } else {
    set workdir $DefaultMDFFWorkDir
  }
  
  if { [info exists arg(jobname)] } {
    set jobname $arg(jobname)
  } else {
    set jobname [file rootname [file tail $pdb]]
  }
 
  if { [info exists arg(fixed)] } {
    set fixed $arg(fixed)
  } else {
    set fixed $DefaultFixed
  }
 
  if { [info exists arg(gscale)] } {
    set gscale $arg(gscale)
  } else {
    set gscale $DefaultGScale
  }
 
  if { [info exists arg(minsteps)] } {
    set minsteps $arg(minsteps)
  } else {
    set minsteps $DefaultMinSteps
  }
 
  if { [info exists arg(numsteps)] } {
    set numsteps $arg(numsteps)
  } else {
    set numsteps $DefaultNumSteps
  }
 
  if { [info exists arg(bestN)] } {
    set bestN $arg(bestN)
  } else {
    set bestN $DefaultBestN
  }
  
  if { [info exists arg(chseg)] } {
   ### input file structure: name chain segname 
    set inputfile [open $arg(chseg) r]
    set data [read $inputfile]
    close $inputfile
    set chseg [split $data "\n"]
  } else {
    set mol [mol new $pdb]
    set sel [atomselect $mol "all"]
    set chains [$sel get chain]
    set segnames [$sel get segname]
    #get unique list
    foreach chain $chains {dict set tmp $chain 1}
    set chains [dict keys $tmp]
    foreach segname $segnames {dict set tmp2 $segname 1}
    set segnames [dict keys $tmp2]

    foreach chain $chains segname $segnames {
      lappend chseg  [list $pdb $chain $segname]
    }
  
  }
  
  if { [info exists arg(topfiles)] } {
    set topfiles $arg(topfiles)
  } else {
    set topfiles $DefaultTopFiles
  }
  
  if { [info exists arg(parfiles)] } {
    set parfiles $arg(parfiles)
  } else {
    set parfiles $DefaultParFiles
  }
  
  if { [info exists arg(dcdfreq)] } {
    set dcdfreq $arg(dcdfreq)
  } else {
    set dcdfreq $DefaultDCDFreq
  }
 
  if { [info exists arg(namdargs)] } {
    set namdargs $arg(namdargs)
  } else {
    set namdargs ""
  }
  
  if { [info exists arg(ss)] } {
    set ss $arg(ss)
  } else {
    set ss ""
  }


  if { [info exists arg(cis)] } {
    set cis $arg(cis)
  } else {
    set cis ""
  }

  if { [info exists arg(chir)] } {
    set chir $arg(chir)
  } else {
    set chir ""
  }

  if { [info exists arg(extrab)] } {
    set extrab  ""
    foreach extrab_file $arg(extrab) {
      lappend extrab [file tail $extrab_file]
    }
  } else {
    set extrab ""
  }
  
  set MOL [file rootname [file tail $pdb]]
  set mapname $density
  
  if { [file exists $workdir] } {
    puts "The working directory already exists!"
    exit 1
  }


  file mkdir $workdir
	
   if { [info exists arg(extrab)] } {
    foreach extrab_file $arg(extrab) {
      file copy -force $extrab_file $workdir
    }
  }
  

  set mutations ""
	auto_makepsf $pdb $topfiles $chseg $mutations
  file rename -force ${MOL}.psf $workdir
  file rename -force ${MOL}-psfout.pdb $workdir

  for {set i 0} {$i < [llength $mapname]} {incr i} {
    #make _potential
    set name [lindex $mapname $i]
    set sel [lindex $gridpdbsel $i]
    mdff griddx -i $name -o [file join $workdir ${name}_potential.dx]
    lappend potentials "${name}_potential.dx"
    #make grid
    mdff gridpdb -psf [file join $workdir ${MOL}.psf] -pdb [file join $workdir ${MOL}-psfout.pdb] -o [file join $workdir ${MOL}-${name}-psfout-grid.pdb] -seltext $sel
    lappend gridpdbs "${MOL}-${name}-psfout-grid.pdb"
  }

	#make ssrestraints
	if { $ss == ""} {
    ssrestraints -psf [file join $workdir ${MOL}.psf] -pdb [file join $workdir ${MOL}-psfout.pdb] -o [file join $workdir ${MOL}-ssrestraints.txt] -hbonds
    lappend extrab ${MOL}-ssrestraints.txt
  } else {
    lappend extrab $ss
  }
	
  #make cispeptide
  if { $cis == ""} {
	  mol delete all
	  mol new [file join $workdir ${MOL}.psf]
	  mol addfile [file join $workdir ${MOL}-psfout.pdb]
	  cispeptide restrain -o [file join $workdir ${MOL}-cispeptide.txt]
    lappend extrab ${MOL}-cispeptide.txt
  } else {
    lappend extrab $cis
  }

	#chirality
  if { $chir == ""} {
	  chirality restrain -o [file join $workdir ${MOL}-chirality.txt]
    lappend extrab ${MOL}-chirality.txt
  } else {
    lappend extrab $chir
  }

	#fix pdb
	set all [atomselect top all]
	$all set occupancy 0.0
	set restraint [atomselect top "$fixed"]
	$restraint set occupancy 1.0
	$all writepdb [file join $workdir fixed.pdb]
	
  
  		
  mdff setup -o $jobname -psf ${MOL}.psf -pdb ${MOL}-psfout.pdb -griddx $potentials -gridpdb $gridpdbs -extrab $extrab -gscale $gscale -minsteps $minsteps -numsteps $numsteps -fixpdb fixed.pdb -dir $workdir -parfiles $parfiles
			
  
  set frpdb [open [file join $workdir mdff_template.namd] "r"]
  set spdb [read $frpdb]
  close $frpdb
  set fwpdb [open [file join $workdir mdff_template.namd] "w"]
  regsub -all "dcdfreq *\n" $spdb "dcdfreq ${dcdfreq}\n" spdb
  puts $fwpdb $spdb
  close $fwpdb
  
  puts "Starting NAMD with job $jobname"
  exec namd2 $namdargs [file join $workdir ${jobname}-step1.namd] > [file join $workdir mdff_${jobname}-step1.log]
  puts "NAMD finished"

	set inmol [mol new [file join $workdir ${MOL}-psfout.pdb]]
  set results [mol new [file join $workdir ${MOL}.psf]]
  mol addfile [file join $workdir ${jobname}-step1.dcd] waitfor all
  regen_chains $inmol $results
  set sel [atomselect $results all]
  $sel frame last
  format_pdb $sel [file join $workdir ${jobname}-step1-result.pdb]
}

proc ::MODELMAKER::regen_chains {INMOL OUTMOL} {
  set sel [atomselect $OUTMOL "all"]
  set segnamesout [$sel get segname]
  #get unique list
  foreach segname $segnamesout {dict set tmp2 $segname 1}
  set segnamesout [dict keys $tmp2]
  
  foreach segnameout $segnamesout {
    set outsel [atomselect $OUTMOL "segname $segnameout"]
    set insel [atomselect $INMOL "segname $segnameout"]
    $outsel set chain [lindex [$insel get chain] 0]
  }
}

proc ::MODELMAKER::regen_segnames {INMOL OUTMOL} {
  
  set sel [atomselect $INMOL "all"]
  set segnamesin [$sel get segname]
  #get unique list
  foreach segname $segnamesin {dict set tmp2 $segname 1}
  set segnamesin [dict keys $tmp2]

  
  set sel [atomselect $OUTMOL "all"]
  set segnamesout [$sel get segname]
  #get unique list
  foreach segname $segnamesout {dict set tmp2 $segname 1}
  set segnamesout [dict keys $tmp2]
  
  foreach segnamein $segnamesin segnameout $segnamesout {
    if {$segnamein != "" && $segnameout != ""} {
      set sel [atomselect $OUTMOL "segname $segnameout"]
      $sel set segname $segnamein
    }
  }
}

proc ::MODELMAKER::cccolor_usage { } {
  variable DefaultThreshold  
  puts "Usage: modelmaker cccolor -pdb <input pdb file> -density <input density file> \
    -res <resolution of density in Angstroms> ?options?"
  puts "Options:"
  puts "  -threshold    <ignore density below threshold value> (default: $DefaultThreshold)>"
  puts "  -spacing      <grid spacing of simulated map> (default: matches input map)>"
  puts "  -sel          <selection text to calculate correlations on> (default: protein and noh)>"
  puts "  -mode         <calculate correlations for secondary structure (ss) or individual residue (res)> (default: ss)>"
  puts "  -ressel       <selection text of resid range if performing per-resiude calculation (e.g., resid 1 to 20) > (default: all)>"
  puts "  -bbonly       <calculate per-residue correlation for backbone atoms only (options: on or off)> (default: off)>"

}

proc ::MODELMAKER::cccolor { args } {
  variable DefaultThreshold  
  
  global env

  set nargs [llength [lindex $args 0]]
  if {$nargs == 0} {
    cccolor_usage
    error ""
  }
  
  foreach {name val} $args {
    switch -- $name {
      -pdb { set arg(pdb) $val }
      -density { set arg(density) $val }
      -res { set arg(res) $val }
      -threshold { set arg(threshold) $val }
      -spacing { set arg(spacing) $val }
      -sel { set arg(sel) $val }
      -mode { set arg(mode) $val }
      -ressel { set arg(ressel) $val }
      -bbonly { set arg(bbonly) $val }
      default { puts "Unknown argument $name"; return  }
    }
  }
  
  if { [info exists arg(pdb)] } {
    set pdb $arg(pdb)
  } else {
    error "input pdb file required!"
  }
  
  if { [info exists arg(density)] } {
    set density $arg(density)
  } else {
    error "density map required!"
  }
  
  if { [info exists arg(res)] } {
    set res $arg(res)
  } else {
    error "density map resolution required!"
  }
  
  if { [info exists arg(threshold)] } {
    set threshold $arg(threshold)
  } else {
    set threshold $DefaultThreshold
  }
  
  if { [info exists arg(spacing)] } {
    set spacing $arg(spacing)
  } else {
    set spacing -1
  }
 
  if { [info exists arg(sel)] } {
    set sel $arg(sel)
  } else {
    set sel "protein and noh"
  }
  
  if { [info exists arg(mode)] } {
    set mode $arg(mode)
  } else {
    set mode "ss"
  }
  
  if { [info exists arg(ressel)] } {
    set ressel $arg(mode)
  } else {
    set ressel "all"
  }
  
  if { [info exists arg(bbonly)] } {
    set bbonly $arg(bbonly)
  } else {
    set bbonly 0
  }
  
  set ss_on 0
  set res_on 0
  switch $mode {
    "ss" { set ss_on 1  }
    "res" { set res_on 1 }
    default { puts "Unknown correlation mode $mode"; return }
  }
  ::CCColor::cccolor $pdb $density $res $threshold $spacing $ss_on $res_on $sel $ressel "/usr/tmp" $bbonly
}

proc ::MODELMAKER::get_empty_density_usage { } {
  variable DefaultMaskCutoff  
  puts "Usage: modelmaker get_empty_density -pdb <input pdb file> -density <input density file> ?options?"
  puts "Options:"
  puts "  -cutoff    <Distance cutoff in Angstroms> (default: $DefaultMaskCutoff)>"

}

proc ::MODELMAKER::get_empty_density { args } {
  variable DefaultMaskCutoff  
  
  set nargs [llength [lindex $args 0]]
  if {$nargs == 0} {
    get_empty_density_usage
    error ""
  }
  
  foreach {name val} $args {
    switch -- $name {
      -pdb { set arg(pdb) $val }
      -density { set arg(density) $val }
      -cutoff { set arg(cutoff) $val }
      default { puts "Unknown argument $name"; return  }
    }
  }
  
  if { [info exists arg(pdb)] } {
    set pdb $arg(pdb)
  } else {
    error "input pdb file required!"
  }
  
  if { [info exists arg(density)] } {
    set density $arg(density)
  } else {
    error "density map required!"
  }
  
  if { [info exists arg(cutoff)] } {
    set cutoff $arg(cutoff)
  } else {
    set cutoff $DefaultMaskCutoff
  }
  
  set mol [mol new $pdb]
  set sel [atomselect $mol "all"]

  volmap mask $sel -o mask.dx -cutoff $cutoff
  mdff griddx -i mask.dx -o mask_invert.dx
  volutil -mult $density mask_invert.dx -o "[file rootname [file tail $density]]_empty.dx"
  file delete -force mask.dx mask_invert.dx
  mol delete $mol
}

## Calling to the SSAnalysis package to analyze
## secondary structures of a sequence of frames
## Return: the persistence of each secondary structure 
##      * per residue (in percentage) and the representative frame of 
##        the most prevalent secondary structure (average ss)
##  
##      * dcd (and pdb) file containing frames with the average ss
##
##      * dcd (and pdb) file containing frames with the ss pattern searched 
##        by the user
##

proc ::MODELMAKER::ssanalysis_usage {} {

  puts "ssanalysis: secondary structure analysis"
  puts "Usage: ssanalysis \[args...\]"
  puts "Args:"
  puts "  -workdir    (optional) <working/project directory for job> (default: current folder)>"
  puts "  -mol        (optional if -molfile provided) VMD molecule containing the frames to be analyzed"
  puts "  -molfile    (optional if -mol provided) file containing the frames to be analyzed"
  puts "  -strtcfile  (mandatory file if -molfile provided as dcd) structure file containing structure\
   information (e.g. *.psf or *.pdb)"
  puts "  -avg        flag to calculate avg secondary structure <0 or 1> (default 0)"
  puts "  -sel        (mandatory if -avg 1)atom selection text defining region to be analyzed (default:all )"
  puts "  -plot       plot tool to be used to display average ss <multiplot or gnuplot> (default multiplot and gnuplot on text mode)"
  puts "  -showplot   (only when plot = multiplot)show plot with average secondary structure analysis <0 or 1> (default 1)"  
  puts "  -seqfind    flag to search for secondary sequence <0 or 1>"
  puts "  -seq        (mandatory if -seqfind 1) nested list of searching ss element and atom selection pairs. Search for\
    consecutive ss elements in the atom selection. - e.g. {{\"H\" \"resid 1 to 10\",\"C\" \"resid 15 to 22\"}} "
  puts "  -output     prefix to be used in the naming of the output files (default: output)"
 
  return
}

proc ::MODELMAKER::ssanalysis { args } {

  if {[expr fmod([llength $args],2)] > 0.0 || [llength $args] == 0} {
    ::MODELMAKER::ssanalysis_usage
    return
  }
  return [eval ss_analysis $args] 
}

proc ::MODELMAKER::model_usage { } {
  variable DefaultModelSel
  variable DefaultModelN
  puts "Usage: modelmaker model -template <input pdb file> -fasta <input fasta file> -n <number of models to generate> ?options?"
  puts "Options:"
  puts "  -sel      <selection text for modeling> (default: $DefaultModelSel)>"
  puts "  -helix    <list of start and end resids for helix restraints> (default: none)>"
  puts "  -n        <number of models to generate> (default: $DefaultModelN)>"
  puts "  -o        <name prefix of generated models> (default: Taken from template)>"
  puts "  -resstart <starting residue number to use for output> (default: Taken from template)> "
  puts "  -chain    <chain id of generated models> (default: Taken from template)>"
  puts "  -seg      <segname of generated models> (default: Taken from template)>"
  
}

proc ::MODELMAKER::model { args } {
  variable DefaultModelSel
  variable DefaultModelN
  
  set nargs [llength [lindex $args 0]]
  if {$nargs == 0} {
    model_usage
    error ""
  }
  
  foreach {name val} $args {
    switch -- $name {
      -template { set arg(template) $val }
      -fasta { set arg(fasta) $val }
      -sel { set arg(sel) $val }
      -helix { set arg(helix) $val }
      -n { set arg(n) $val }
      -o { set arg(o) $val }
      -resstart { set arg(resstart) $val }
      -chain { set arg(chain) $val }
      -seg { set arg(seg) $val }
      default { puts "Unknown argument $name"; return  }
    }
  }
  
  if { [info exists arg(template)] } {
    set template $arg(template)
  } else {
    error "input template pdb file required!"
  }
  
  if { [info exists arg(fasta)] } {
    set fasta $arg(fasta)
  } else {
    error "input fasta file required!"
  }
  
  if { [info exists arg(sel)] } {
    set sel $arg(sel)
  } else {
    set sel $DefaultModelSel
  }
  
  if { [info exists arg(n)] } {
    set n $arg(n)
  } else {
    set n $DefaultModelN
  }
  
  if { [info exists arg(o)] } {
    set outputname $arg(o)
  } else {
    set outputname [file rootname [file tail $arg(template)]]
  }
  
  if { [info exists arg(chain)] } {
    set outchain $arg(chain)
  } else {
    set outchain "none"
  }
  
  if { [info exists arg(seg)] } {
    set outseg $arg(seg)
  } else {
    set outseg "none"
  }
  
  if { [info exists arg(resstart)] } {
    set resstart $arg(resstart)
  } else {
    set resstart "none"
  }
  
  if { [info exists arg(helix)] } {
    set helix "-helix $arg(helix)"
  } else {
    set helix ""
  }

  
  set tempmol [mol new $template]
  #renaming info
  set temp_sel [atomselect $tempmol all]
  if {$resstart == "none"} {
    set resstart [lindex [lsort -integer [$temp_sel get resid]] 0]
  }
  if {$outchain == "none"} {
    set outchain [lindex [$temp_sel get chain] 0]
  }
  if {$outseg == "none"} {
    set outseg [lindex [$temp_sel get segname] 0]
  }

  #modeling info
  set worksel [atomselect $tempmol $sel]
  
  set resids [$worksel get resid]
  set startres [lindex $resids 0]
  set endres [lindex $resids end]
  $worksel writepdb "tmppdb.pdb"
  
  if {$sel != "all"} {
    modelmaker seqsub -i $fasta -o "tmp.fasta" -start $startres -end $endres
    set infile [open "tmp.fasta" r]
    set file_data [read $infile]
    close $infile
  } else {
    set infile [open $fasta]
    set file_data [read $infile]
    close $infile 
  }
 
  set outfile [open tmpseq.seq w]
  #add seq header
  puts $outfile ">P1;tmpseq" 
  puts $outfile "sequence:::::::::"
  #remove fasta header 
  foreach line $file_data  {
    if { [lindex $line 0] != ">" } {
      puts $outfile $line
    }
  }
  puts $outfile "*"
  close $outfile

  set chain [lindex [$worksel get chain] 0]

  exec python $::env(RosettaVMDDIR)/align2d.py -template "tmppdb.pdb" -sequence "tmpseq.seq" -chain $chain
  exec python $::env(RosettaVMDDIR)/model-single.py -template "tmppdb.pdb" -sequence "tmpseq.seq" -alignment "alignment.aln" \
     -n $n {*}$helix > modeller.log 
  
  set fp [open "modeller.log" r]
  while {[gets $fp line] >= 0} {
    if {[regexp "tmpseq.*.pdb" [lindex $line 0]]} {
     lappend scores $line 
    }
  }
  close $fp
  set ranked [lsort -index 1 $scores]
  for {set i 0} {$i < [llength $ranked]} {incr i} {
    set struct [lindex $ranked $i]
    set tmpmol [mol new [lindex $struct 0]]
    set tmpsel [atomselect $tmpmol all]
    if {$outchain != ""} {$tmpsel set chain $outchain} 
    if {$outseg != ""} {$tmpsel set segname $outseg}
    renumber $tmpsel $resstart
    $tmpsel writepdb "${outputname}-best-[expr $i+1].pdb"
    #file rename -force [lindex $struct 0] "${outputname}-best-[expr $i+1].pdb"
  }
  
}


proc ::MODELMAKER::fragments_usage { } {
  puts "Usage: modelmaker fragments -fasta <input fasta file> ?options?"
  puts "Options:"
  puts "  -o          <Prefix for fragment file output. Default taken from fasta.> "
  puts "  -np         <Number of processors to use. MPI version only> "
  
}

proc ::MODELMAKER::fragments { args } {
  variable rosettaEXE
  variable rosettaPath
  variable rosettadbpath
  
  set nargs [llength [lindex $args 0]]
  if {$nargs == 0} {
    fragments_usage
    error ""
  }
  
  foreach {name val} $args {
    switch -- $name {
      -fasta { set arg(fasta) $val }
      -np { set arg(np) $val }
      -o { set arg(o) $val }
      default { puts "Unknown argument $name"; return  }
    }
  }
  
  if { [info exists arg(fasta)] } {
    set fasta $arg(fasta)
  } else {
    error "input fasta file required!"
  }
  
  if { [info exists arg(np)] } {
    set mpinp $arg(np)
  } elseif { [string match "*mpi.*" $rosettaEXE] }  {
    error "number of processors (-np) must be specified when using MPI version of Rosetta!"
  }
  
  if { [info exists arg(o)] } {
    set output_prefix $arg(o)
  } else {
    set output_prefix [file rootname $fasta]
  }
  
  if { [string match "*mpi*" $rosettaEXE] } {
    set mpi_args "mpiexec -np $mpinp"      
  } else {
    set mpi_args ""
  }
  
  exec {*}$mpi_args [glob $rosettaPath/fragment_picker.$rosettaEXE] -in:file:fasta $fasta \
    -in:path:database $rosettadbpath \
    -in:file:vall ${rosettadbpath}/../../tools/fragment_tools/vall.apr24.2008.extended.gz \
    -out:file:frag_prefix $output_prefix \
    -overwrite >> fragments.log
  

}
