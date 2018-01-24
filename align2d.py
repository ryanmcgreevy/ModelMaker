from modeller import *
#import sys
import os.path
import argparse

parser = argparse.ArgumentParser(description=("This script performs a 2d alignment based on " +
"a given structure and sequence using MODELLER."))

parser.add_argument("-template", action="store", dest="TEMPLATE", required=True, 
help='Template structure (.pdb) - required')
parser.add_argument("-sequence", action="store", dest="SEQ", required=True, 
help='Sequence file (.seq) - required')
parser.add_argument("-chain", action="store", dest="CHAIN", required=True, 
help='Chain ID - required')
args = parser.parse_args()

inmodel = args.TEMPLATE
chain = args.CHAIN
sequence = args.SEQ


model_name = os.path.splitext(inmodel)[0]
sequence_name = os.path.splitext(sequence)[0]

env = environ()
aln = alignment(env)
mdl = model(env, file=inmodel, model_segment=('FIRST:%s' % chain,'LAST:%s' % chain))
aln.append_model(mdl, align_codes=model_name, atom_files=inmodel)
aln.append(file=sequence, align_codes=sequence_name)
aln.align2d()
aln.write(file='alignment.aln', alignment_format='PIR')
aln.write(file='alignment.pap', alignment_format='PAP')
