#Python script for structural clustering of MD trajectory.
#Written by Keith Cassidy <ckcassidy1@gmail.com>
#On TCBG machines: source /home/ccassid2/.bashrc

import numpy
from numpy import vstack,array
import csv
import os
import glob
import time
import argparse
from scipy.spatial.distance import squareform

parser = argparse.ArgumentParser(description=("This script uses K-medoids to cluster structures " +
"from an MD trajectory based on pairwise RMSD. The user must provide a .psf file (-psf) and "     +
"a .dcd trajectory file (-dcd). Optionally the user may specify an atom selection (-sel). In "    +
"addition, the user may provide a precomputed, square distance matrix (-dist) or write the "      +
"distance matrix to a text file (-o).  For K-medoids, the user may specify the min (-kmin) and "  + 
"max (-kmax) number of clusters to consider for the silhouette calculation as well as the number "+
"of passes for each k-medoids iteration (-p). If kmin = kmax, the silhouette score will not be "  +     
"computed. This script returns text files containing (1) the silhouette score, if computed, for " +
"each cluster (clu_silhouette.txt), (2) the structures corresponding to each of the K medoids "   +
"(clu_medoids.txt), (3) the cluster label for all structures (clu_labels.txt), and (4) the "      +
"structures contained in each of the k individual clusters (clu_{K}_frames.txt). Note, this "     +
"script overwrites output from previous clustering runs and requires the Numpy, Scipy, "          +
"MDAnalysis, Pycluster, and scikit-learn packages."))
parser.add_argument("-psf", action="store", dest="PSF", required=True, 
help='Structure file (.psf) - required')
parser.add_argument("-dcd", action="store", dest="DCD", required=True, 
help='Trajectory file (.dcd) - required')
parser.add_argument("-sel", nargs="?", default="name CA", type=str, dest="SEL", 
help='Atom selection - default: "name CA"')
parser.add_argument("-dist", nargs="?", default=None, type=str, dest="input", 
help='Provide square distance matrix - default: None')
parser.add_argument("-kmin", nargs="?", default=2, type=int,dest="kmin", 
help='Min number of clusters for silhouette calc. - default: 2')
parser.add_argument("-kmax", nargs="?", default=10, type=int,dest="kmax", 
help='Max number of clusters for silhouette calc. - default: 10')
parser.add_argument("-p", nargs="?", default=2000, type=int, dest="passes", 
help='Number of iterations per k-medoids run - default: 2000')
parser.add_argument("-o", action="store_true", dest="write_dist", required=False, 
help='Write distance matrix to text file? - optional')
args = parser.parse_args()

import matplotlib.pyplot as plt
from MDAnalysis.analysis.align import *
import MDAnalysis.lib.qcprot as qcp
from MDAnalysis import *
from Pycluster import kmedoids
from sklearn.metrics import silhouette_samples, silhouette_score

PSF = args.PSF
DCD = args.DCD
SEL = args.SEL
input_dist = args.input
clus_min = int(args.kmin)
clus_max = int(args.kmax)
passes = int(args.passes)
write_dist = args.write_dist

print input_dist

#Delete cluster assignment files from previous clusterings.
for filePath in glob.glob("clu_*.txt"):
	if os.path.isfile(filePath):
		os.remove(filePath)

print "\n>>> Python script for RMSD-based structural clustering of MD trajectory."
print "\n>>> Direct questions to Keith Cassidy: <ckcassidy1@gmail.com>"
print "\n>>> Please cite:"
print "\tStone JE, Perilla JR, Cassidy CK, Schulten K,"
print "\t\"GPU-accelerated molecular dynamics clustering analysis with openACC\","
print "\tIn Parallel Programming with openACC, Elsevier, 2016."

#Record runtime
start_time = time.time()

#Get trajectory information.
print "\n>>> Loading structure and trajectory files..."
traj = Universe(PSF,DCD)
traj2 = Universe(PSF,DCD)

atoms = traj.select_atoms("%s" % SEL)
N = len(atoms)
print "\tNo. atoms: %s" % N

frames = traj.trajectory
nframes = len(frames)
print "\tNo. frames: %s" % nframes

pairs = nframes*(nframes-1)/2
data = numpy.zeros(pairs)

if input_dist == None:
#Compute condensed distance matrix.
	print "\n>>> Computing pairwise RMSD distance matrix..."
	count = 0
	count2 = 0
	for js in traj2.trajectory[0:nframes-1:1]: 
		ref = traj2
		select = '%s' % SEL
		selections = {'reference': select, 'target': select}

		ref_atoms = ref.select_atoms(selections['reference'])
		traj_atoms = traj.select_atoms(selections['target'])
		natoms = len(traj_atoms)

		ref_com = ref_atoms.center_of_mass()
		ref_positions = ref_atoms.positions - ref_com
		traj_positions = traj_atoms.positions.copy()

		for ts in traj.trajectory[int(count+1):int(nframes):1]:
		    x_com = traj_atoms.center_of_mass()
		    traj_positions[:] = traj_atoms.positions - x_com
		    R = np.zeros((9,), dtype=np.float64)
		    a = ref_positions.T.astype('float64')
		    b = traj_positions.T.astype('float64')
		    data[count2] = qcp.CalcRMSDRotationalMatrix(a, b, natoms, R, None)
		    count2 += 1
		count += 1
	print "\tDone."

	#Create redundant distance matrix for K-Medoids function
	data = squareform(data)

else:
	print "\n>>> Loading provided distance matrix..."
	data = 	np.loadtxt('%s' % input_dist)

if write_dist:
	print "\n>>> Writing distance_matrix.txt..."
	np.savetxt('distance_matrix.txt', data,fmt='%10.3f')

if clus_min != clus_max:
	#Compute silhouette score for k=kmin -> k=kmax:
	print "\n>>> Determining optimal number of clusters using silhouette score..."
	sil = numpy.zeros(clus_max-clus_min+1)
	for i in range(clus_min,clus_max+1):
		idx,error,nfound = kmedoids(data,nclusters=i,npass=passes)
		sil[i-clus_min] = silhouette_score(data,idx, metric='precomputed',
		sample_size=None,random_state=None)

	clus_num = numpy.argmax(sil)+clus_min
	print "\tOptimal number of clusters: %s" % clus_num
else:
	sil = 	sil = numpy.zeros(clus_max-clus_min+1)	
	clus_num = clus_min
	print "\n>>> Requested %s clusters" % clus_num

#K-Medoids clustering
idx,error,nfound = kmedoids(data,nclusters=clus_num,npass=passes)
medoids = numpy.unique(idx)
for i in range(0,clus_num):
	idx[idx==medoids[i]] = i

print "\n>>> K-medoids clustering results:"
for i in range(0,clus_num):
	print "\tcluster_%s: %s pts (%1.2f percent), medoid: %s" %(i,len(data[idx==i,0]),
	100*float(len(data[idx==i,0]))/len(idx),medoids[i])

#Write output files.
print "\n>>> Writing output files..."
#print "\tsilhouette.txt\n\tclu_labels.txt\n\tclu_medoids.txt\n\tclu_{k}_frames.txt"

#Write silhouette scores for tested clusters.
if len(sil) > 1:
	with open("clu_silhouette.txt", 'w+') as f:
		for j in range(clus_min,len(sil)+clus_min):
			writer = csv.writer(f,delimiter=" ")	
			writer.writerow([j,sil[j-clus_min]])

#Write files containing individual cluster assignments.
num = numpy.arange(0,len(idx),1)
txt = (numpy.vstack((idx,num))).T
for j in range(0,clus_num):
	with open("clu_%s_frames.txt" % j, 'w+') as f:
		writer = csv.writer(f)
		for k in range(0, len(idx)):
			if int(txt[k][0]) == j:
				writer.writerow([txt[k][1]])

#Write all cluster assignments to single file.
numpy.savetxt('clu_labels.txt',txt,fmt='%i')

#Write medoid frames
num = numpy.arange(0,clus_num,1)
med = (numpy.vstack((num,medoids))).T
numpy.savetxt('clu_medoids.txt',med,fmt='%i')

print "\n>>> Finished calculation. Runtime: %1.2f seconds.\n" % (time.time() - start_time)
