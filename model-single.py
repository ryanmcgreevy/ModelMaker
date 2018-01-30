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
parser.add_argument("-helix", action="store", dest="HELIX", nargs='+', required=False, default='', 
help='List of residue ranges for alpha helix restraints. - optional')
parser.add_argument("-n", action="store", dest="N", type=int, required=False, default=1, 
help='Number of models to generate. Default: 1')
args = parser.parse_args()

template = os.path.splitext(args.TEMPLATE)[0]
seq = os.path.splitext(args.SEQ)[0]
align = args.ALIGN
helix = args.HELIX
n = args.N
log.verbose()
env = environ()

class MyModel(automodel):
    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms
#       Add some restraints from a file:
#       rsr.append(file='my_rsrs1.rsr')

#       Residues x through y should be an alpha helix:
        if helix:
          i = 0
          while i < len(helix):
            rsr.add(secondary_structure.alpha(self.residue_range('%s:' % helix[i], '%s:' % helix[i+1])))
            i += 2

a = MyModel(env, alnfile=align, # alignment filename
              knowns=template, # codes of the templates
              sequence=seq, # code of the target
              assess_methods=(	assess.DOPE,
                            	assess.DOPEHR,
                            	#soap_protein_od.Scorer(),
                            	assess.GA341))
a.starting_model = 1 	# index of the first model
a.ending_model = n		# index of the last model
a.make()				# do homology modeling
