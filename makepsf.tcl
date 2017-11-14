# PSF GENERATION SCRIPT
# Maximilian Scheurer, April 2016
# bla
#########
# Input #
#########
namespace eval ::MakePsf {
# Set up Variable
  namespace export auto_makepsf
  set version 0.1
  set packageDescription "MakePsf Plugin for RosettaVMD protocol"
 
    # Variable for the path of the script
    variable home [file join [pwd] [file dirname [info script]]]
}
package provide MakePsf $MakePsf::version

proc auto_makepsf {args} \
{
  return [eval ::MakePsf::auto_makepsf $args]
}

proc ::MakePsf::auto_makepsf {MOL topnames ch_seg mutations} \
{
  ##########
  # !!!!!! #
  ##########
    puts "running auto_makepsf"
    mol new $MOL
    set molname  [file rootname [file tail $MOL]]
#    puts "$ch_seg"

    package require psfgen
    resetpsf

    foreach top $topnames {
      topology $top
    }


    #########################
    # Define common aliases #
    #########################
    pdbalias residue G GUA
    pdbalias residue C CYT
    pdbalias residue A ADE
    pdbalias residue T THY
    pdbalias residue U URA

    foreach bp { GUA CYT ADE THY URA } {
     pdbalias atom $bp "O5\*" O5'
     pdbalias atom $bp "C5\*" C5'
     pdbalias atom $bp "O4\*" O4'
     pdbalias atom $bp "C4\*" C4'
     pdbalias atom $bp "C3\*" C3'
     pdbalias atom $bp "O3\*" O3'
     pdbalias atom $bp "C2\*" C2'
     pdbalias atom $bp "O2\*" O2'
     pdbalias atom $bp "C1\*" C1'
   }

   pdbalias atom ILE CD1 CD
   pdbalias atom SER HG HG1
   pdbalias residue MSE MET
   pdbalias residue HIS HSD 

   #  pdbalias residue HIS HSD for delta protonated
   #  pdbalias residue HIS HSE for epsilon protonated
   #  pdbalias residue HIS HSP for dpuble protonated


   # Heme aliases
   pdbalias residue HEM HEME
   pdbalias atom HEME "N A" NA
   pdbalias atom HEME "N B" NB
   pdbalias atom HEME "N C" NC
   pdbalias atom HEME "N D" ND

   # Water aliases
   pdbalias residue HOH TIP3
   pdbalias atom TIP3 O OH2

   # Ion aliases
   pdbalias residue K POT
   pdbalias atom K K POT
   pdbalias residue ICL CLA
   pdbalias atom ICL CL CLA
   pdbalias residue INA SOD
   pdbalias atom INA NA SOD
   pdbalias residue CA CAL
   pdbalias atom CA CA CAL
   pdbalias residue ZN ZN2


   # Other aliases
   pdbalias atom LYS 1HZ HZ1
   pdbalias atom LYS 2HZ HZ2
   pdbalias atom LYS 3HZ HZ3

   pdbalias atom ARG 1HH1 HH11
   pdbalias atom ARG 2HH1 HH12
   pdbalias atom ARG 1HH2 HH21
   pdbalias atom ARG 2HH2 HH22

   pdbalias atom ASN 1HD2 HD21
   pdbalias atom ASN 2HD2 HD22

   #pdbalias residue PDBRESNAME TOPRESNAME changes to Resname from PDB into the one of TOPFILE
   #pdbalias atom RESNAME PDBATOMNAME TOPATOMNAME changes to atom name from PDB into the one of TOPFILE

   ####################################
   # Split input files for each chain #
   ####################################

   ### input file structure: name chain segname 
#   set inputfile [open "${MOL}_chains.txt" r]
#   set data [read $inputfile]
#   close $inputfile
#   set ch_seg [split $data "\n"]
   
   foreach var $ch_seg {
    puts "name:[lindex $var 0] chain:[lindex $var 1] segname:[lindex $var 2]" 
    set chainID [lindex $var 1]
    set segName [lindex $var 2]
    set sel1 [atomselect top "segname $segName"]
    # set sel1 [atomselect top "chain $chainID"]
    foreach selatom $sel1 {
      $selatom set segid $segName 
    }
    $sel1 writepdb $segName.pdb

    set segMutations [lsearch -index 0 -all -exact -inline $mutations $segName]

    ##################################
    # Add Hydrogens by Topology file #
    ##################################

    segment $segName {
      # first NTER
      # last CTER
      # auto angles
      # auto dihedrals
      pdb $segName.pdb
      foreach mutation $segMutations {
        mutate [lindex $mutation 1] [lindex $mutation 2]
      }
      # mutate 303 HSE
    }
    ##################################
    # Overwrite labels with new ones #
    ##################################
    coordpdb $segName.pdb $segName
    file delete $segName.pdb
  }

  ###############################################
  # Complete all missing atoms by Topology file #
  ###############################################
  guesscoord
  ####################################
  # Output modified pdb and psf file #
  ####################################
  #writepdb $MOL-psftemp.pdb

  writepdb $molname-psfout.pdb
  writepsf $molname.psf
  # Rafael first psf file before guesscoord???
  mol new $molname.psf
  mol addfile $molname-psfout.pdb

  foreach var $ch_seg {
    puts "name:[lindex $var 0] chain:[lindex $var 1] segname:[lindex $var 2]" 
    set chain [lindex $var 1]
    set sel [atomselect [molinfo top] "segname [lindex $var 2]"]
    $sel set chain $chain
  }
  set tot [atomselect [molinfo top] "all"]
  $tot writepdb $molname-psfout.pdb
}
