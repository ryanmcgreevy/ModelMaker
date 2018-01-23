from modeller import *
import sys
import os.path

inmodel = sys.argv[1]
chain = sys.argv[2]
sequence = sys.argv[3]


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
