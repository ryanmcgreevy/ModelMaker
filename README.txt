ModelMaker 0.1, a VMD plugin for creating complete models with Rosetta

ModelMaker is an all-in-one modeling solution combining the strength of established 
but often complicated to use modeling suits like Rosetta and Modeller with 
MDFF in an easy-to-use manner in VMD. Integrative modeling approaches usually 
aim to automate the process of structure analysis to avoid human bias, yet the 
experience of structural biologists may be a desirable factor in structure refinement. 
Therefore, ModelMaker uniquely allows for incorporation of user expertise by taking advantage of 
interactive MDFF, allowing for manual manipulation of the structure.

To use ModelMaker, add the following to your ~/.vmdrc:

set env(ROSETTADB) "/your/path/to/rosetta/main/database"
set env(ROSETTAPATH) "/your/path/to/rosetta/main/source/bin"
set auto_path [linsert $auto_path 0 /your/path/to/modelmaker]
package require modelmaker

Make sure to replace the paths with the correct ones specific to the location of
your Rosetta and ModelMaker folders.

In VMD's TkConsole, type 'modelmaker' to see a list of commands. You can also
type 'modelmaker command' where 'command' is one of the available ModelMaker commands
to see more information about it. An example of how to use ModelMaker to fold
protein termini from the tutorial http://www.ks.uiuc.edu/Training/Tutorials/science/rosetta-mdff/rosetta-mdff-tutorial-html/node4.html is below:


modelmaker abinitio -model rpn11_yeast_23-306_complete.pdb -anchor 1 -fragfiles {{rpn11_yeast_23-306_frag9 rpn11_yeast_23-306_frag3}} -sel "resid 196 to 284" -nstruct 1

modelmaker analyze -model rpn11_yeast_23-306_complete.pdb -nstruct 1 -align_template "resid 100 to 195" -comps {{ss 196 284 "A"}}
