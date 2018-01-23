from modeller import *
from modeller.automodel import * # Load the automodel class
import argparse
import os.path

parser = argparse.ArgumentParser(description=("This script generates a complete homology model " +
"using MODELLER."))

parser.add_argument("-template", action="store", dest="TEMPLATE", required=True, 
help='Template structure (.pdb) - required')
parser.add_argument("-sequence", action="store", dest="SEQ", required=True, 
help='Sequence file (.seq) - required')
parser.add_argument("-alignment", action="store", dest="ALIGN", required=True, 
help='Alignment file (.aln) - required')
args = parser.parse_args()

template = os.path.splitext(args.TEMPLATE)[0]
seq = os.path.splitext(args.SEQ)[0]
align = args.ALIGN

log.verbose()
env = environ()

class MyModel(automodel):
    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms
#       Add some restraints from a file:
#       rsr.append(file='my_rsrs1.rsr')

#       Residues x through y should be an alpha helix:
        rsr.add(secondary_structure.alpha(self.residue_range('8:', '17:')))
        rsr.add(secondary_structure.alpha(self.residue_range('21:', '45:')))



a = MyModel(env, alnfile=align, # alignment filename
              knowns=template, # codes of the templates
              sequence=seq, # code of the target
              assess_methods=(	assess.DOPE,
                            	assess.DOPEHR,
                            	#soap_protein_od.Scorer(),
                            	assess.GA341))

a.starting_model = 1 	# index of the first model
a.ending_model = 10		# index of the last model
a.make()				# do homology modeling

#class MyModel(automodel):
#    def special_restraints(self, aln):
#        rsr = self.restraints
#        at = self.atoms
##       Add some restraints from a file:
##       rsr.append(file='my_rsrs1.rsr')
#
##       Residues x through y should be an alpha helix:
#        rsr.add(secondary_structure.alpha(self.residue_range('8:', '17:')))
#        rsr.add(secondary_structure.alpha(self.residue_range('21:', '45:')))



#a = MyModel(env, alnfile='%s_yeast_mouse.aln' % subunit, # alignment filename
#              knowns='%s_yeast'% (subunit), # codes of the templates
#              sequence='%s_mouse'% (subunit), # code of the target
#              assess_methods=(	assess.DOPE,
#                            	assess.DOPEHR,
#                            	#soap_protein_od.Scorer(),
#                            	assess.GA341))

#a.starting_model = 1 	# index of the first model
#a.ending_model = 10		# index of the last model
#a.make()				# do homology modeling
