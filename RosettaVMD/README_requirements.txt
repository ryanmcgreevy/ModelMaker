Author: Maximilian Scheurer
Supervisor: Till Rudack

Software requirements:
VMD installation
NAMD installation
Situs 2.8 (http://situs.biomachina.org)
Resmap (for creation of resolution information maps, http://resmap.sourceforge.net)
gnuplot
Python 2.7 (for helix registration, currently not working)

Config file:
The config file is a tcl file, steering VMD in order to manage
Rosetta, NAMD for MDFF and VMD for system preparation and analysis.

A1) Linking the tcl package to the file is done by:
set packagePath "/path/to/package/index"
lappend ::auto_path $packagePath
package require RosettaVMD

A2) For highresolution MDFF (adapt gridpdb to local resolution), specify
package require DensitySelector

To use Rosetta, NAMD, VMD and all analysis tools correctly, the following variables ought to be set:

topdir: directory with topology files
topfiles: list of topology files in the topdir

parfiles: list of charmm parameter files in the topdir

path: path to the folder containing the namd2 binary
namdArgs: additional namd2 arguments, such as the number of cores "+p4" or "+setcpuaffinity" for the CUDA version
dcdfreq: dcd stride for NAMD

vmdexe: executable path of the given vmd installation
gnuplotexe: executable path of gnuplot

rosettapath: path to the Rosetta binary folder
platform: platform the Rosetta binaries have been compiled on: either "macosclangrelease" for Mac OS X or "linuxgccrelease" for Linux
rosettaDBpath: path to the Rosetta database

ch_seg: chain and segnames for psf generation (list, can be read from a file)
mutations: mutations for specified residues in the PDB (list, can be read from a file)
	if no mutations are to be included, specify mutations as an empty list as: "set mutations {}"




