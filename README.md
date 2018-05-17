# ModelMaker
Version 0.1, a VMD plugin for creating complete models with _Rosetta_

ModelMaker is an all-in-one modeling solution combining the strength of established
but often complicated to use modeling suits like _Rosetta_ and _Modeller_ with
MDFF in an easy-to-use manner in VMD. Integrative modeling approaches usually
aim to automate the process of structure analysis to avoid human bias, yet the
experience of structural biologists may be a desirable factor in structure refinement.
Therefore, ModelMaker uniquely allows for incorporation of user expertise by taking advantage of
interactive MDFF, allowing for manual manipulation of the structure.

To use ModelMaker, add the following to your `~/.vmdrc`:
```
set env(ROSETTADB) "/your/path/to/rosetta/main/database"
set env(ROSETTAPATH) "/your/path/to/rosetta/main/source/bin"
set auto_path [linsert $auto_path 0 /your/path/to/modelmaker]
package require modelmaker
```

Make sure to replace the paths with the correct ones specific to the location of
your Rosetta and ModelMaker folders.

You can also set the ``env(ROSETTAEXE)`` environment variable to use a specific
set of Rosetta binaries. This environment variable should be set to the binary prefix
you wish to use, but will be set by default depending on your platform. (e.g., on OSX,
this environment variable will usually be 
``static.macosclangrelease``)


If you want to use the parallel MPI version of Rosetta, you must compile the MPI binaries with
the following command:

```
./scons.py bin mode=release extras=mpi
```
For more information please see the [Rosetta Documentation](https://www.rosettacommons.org/docs/latest/rosetta_basics/MPI).

Then you will need to set your ``env(ROSETTAEXE)`` environment variable to the prefix
of the MPI binaries (usually mpi.\* where \* is the rest of the specific platform,
like ``macosclangrelease`` or ``linuxgccrelease``) and set the ``-np`` option in the ModelMaker commands
to specify the number of processors to use.

In VMD's TkConsole, type `modelmaker` to see a list of commands. You can also
type `modelmaker command` where 'command' is one of the available ModelMaker commands
to see more information about it. 


An example of how to obtain full length models with ModelMaker from the [tutorial](http://www.ks.uiuc.edu/Training/Tutorials/science/rosetta-mdff/rosetta-mdff-tutorial-html/node4.html) is below:

```
modelmaker full_length_model -template rpn11_yeast_4ocm.pdb -fragfiles {rpn11_yeast_23-306_frag9 rpn11_yeast_23-306_frag3} -fasta rpn11_yeast_23-306.fasta -resstart 23

```

An example of how to use ModelMaker to fold
protein termini from the [tutorial](http://www.ks.uiuc.edu/Training/Tutorials/science/rosetta-mdff/rosetta-mdff-tutorial-html/node4.html) is below:

```
modelmaker abinitio -model rpn11_yeast_4ocm.pdb_full_length.pdb -anchor "resid 23" -fragfiles {{rpn11_yeast_23-306_frag9 rpn11_yeast_23-306_frag3}} -sel "resid 213 to 306" -nstruct 5

modelmaker analyze -model rpn11_yeast_4ocm.pdb_full_length.pdb -template rpn11_yeast_4ocm.pdb_full_length.pdb -nstruct 5 -align_template "resid 23 to 212" -comps {{ss 210 306 "A"}}

```
An example of how to use ModelMaker to model amino acid insertions from the [tutorial](http://www.ks.uiuc.edu/Training/Tutorials/science/rosetta-mdff/rosetta-mdff-tutorial-html/node5.html) is below:

```
modelmaker insertion -model rpn11_yeast_4ocm.pdb_full_length.pdb -fragfiles {rpn11_yeast_23-306_frag9 rpn11_yeast_23-306_frag3} -sel "resid 160 to 179" -nstruct 5 -fasta rpn11_yeast_23-306.fasta -jobname rpn11_insertion

modelmaker analyze -model rpn11_yeast_4ocm.pdb_full_length.pdb -template rpn11_yeast_4ocm.pdb_full_length.pdb -nstruct 5 -align_template "resid 23 to 159 or resid 180 to 216" -comps  {{ss 160 179 "A"} {cluster 160 179 "A" 2}} -jobname rpn11_insertion -insertion yes
```

An example of how to use ModelMaker to refine a model against a mid-resolution density from the [tutorial](http://www.ks.uiuc.edu/Training/Tutorials/science/rosetta-mdff/rosetta-mdff-tutorial-html/node6.html) is below:

```
modelmaker refine -model rpn11_yeast.pdb -anchor "resid 23"  -sel {"resid 212 to 228" "resid 296 to 306"} -density rpn11_model_5_2594_density.mrc -res 7.7 -nstruct 1 -jobname rpn11_refine -score -0.3
```
