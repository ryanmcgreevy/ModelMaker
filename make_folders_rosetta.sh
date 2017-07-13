#!/bin/bash

MOL="rpt5_human_413-439"

# external inputs needed:
# full length sequence fasta
# fragments files
# incomplete pdb file

# Step 1 make rosetta input file if not a file from pdb database
# make_rosetta_input_pdb.tcl

# Step 2 get full model for each of the chains
# cd full_length_model
# use write_chain.tcl to get the relevant chains of PDB databas structure
# define MOL and MOL_template
# ./run_full_length_model.sh for each chain alone if you have more then one chain!
# cp together all chains in one pdb
##################
# !!!!!!!!!!!!!! #
##################
# Attention: fragment file must have exactly the same amino acids as pdb file
# Attention all Heteroatoms have to be deleted Rosetta only works for standard amino acids
# Attention: if the input file is from an MD simulation run make_rosetta_input_pdb.tcl

# Step 3
# define in fold_separate_chains.xml
# - the number of fixed AA <Index name="A_fix" resnums="1-319" />
# - the chains <Chain name="ChainA" chains="A" />
# - the fixed parts <Or name="full_fix" selectors="A_fix,ChainB" />
# - the fragment files for chains in which something is going to be fold <Fragments large_frags="aat000_09_05.200_v1_3" small_frags="aat000_03_05.200_v1_3" selector="ChainA" />
# the template $MOL <RigidChunkCM name="chunk" region_selector="full_fix" template="AB_filled.pdb" selector="full_fix"
#                    apply_to_template="centroid" />v

#Step 4
# cd rosetta_output
#./test_run_fold.sh <JOBNAME> <input pdb>

#Step 5
#./run_cluster_rosetta_fold_termini.sh <JOBNAME> <input pdb>

# Step 6
#aligne
#cd ../full_length_model/
#create a template.pdb
#asigne the rigid region in the template and in the modelled structure that should be aligned
# Attention: Numbering could be different
# Attention: both selection must have exactly the same atoms
# Attention: be awere of atoms with occupancy different to 1
#Dont forget to exclude all other missing pieces from rosetta outputs for modeling
#cd ../scripts/
#vmdis align_rosetta.tcl

#Step 7
#Rosetta Scoring
# Take best 50% of structures
# set MOL name and Num Structures
# vmdis call_rosetta_scoring.tcl

#Step 8
#Cluster
# make the cluster input dcd and remove the translations: +2AA on both sides around the modeled part
#vmdis make_cluster_input.tcl
# identify optimal number of clusters:
# cluster_number.sh <Start Residid> <End Residid> <Max Structure Number>
# get the number of clusters from cluster_<Max Structure Number>_out1.log
#get clusters:
#./clustering.sh <Start Residid> <End Residid> <Max Structure Number> <Number of Clusters>
# write cluster to pdb and representative pdb
# edit mol and max_structure in read_cluster_file.tcl
# vmdis read_cluster_file.tcl

#Step9
#Secondary Structure Analysis

# mkdir scripts
mkdir full_length_model
mkdir mdff_step1_abinitio
mkdir refinment
# mkdir rosetta_input
# mkdir analysis
# mkdir rosetta_output
# mkdir rosetta_output/OUTPUT_FILES
# mkdir rosetta_output/pdb_out
# mkdir rosetta_output/sc_out

###################
# rosetta scripts #
###################

#prepare full length model
# cp /home/trudack/scripts/rosetta/make_rosetta_input_pdb.tcl ./scripts/
cp /home/trudack/scripts/rosetta/make_rosetta_input_pdb.tcl ./full_length_model/
# cp /home/trudack/scripts/rosetta/run_full_length_model.sh  ./scripts/
cp /home/trudack/scripts/rosetta/run_full_length_model.sh ./full_length_model/
# cp /home/trudack/scripts/set_chainID.tcl ./scripts/
cp /home/trudack/scripts/set_chainID.tcl ./full_length_model/
# cp /home/trudack/scripts/rosetta/write_chain.tcl  ./scripts/
cp /home/trudack/scripts/rosetta/write_chain.tcl ./full_length_model/
cp /home/trudack/scripts/rosetta/auto_mdff/* ./mdff_step1_abinitio/


#get structures
cp /home/trudack/scripts/rosetta/rosetta_vmd_cluster_aa/*  ./
cp /home/trudack/scripts/rosetta/rosetta_vmd_cluster/*  ./refinment/
# cp /home/trudack/scripts/rosetta/test_run_fold_termini.sh ./rosetta_output/


#analysis
# cp /home/trudack/scripts/rosetta/call_rosetta_scoring.tcl ./scripts/
# cp /home/trudack/scripts/rosetta/call_rosetta_scoring.tcl ./analysis/
# cp /home/trudack/scripts/rosetta/rosetta_scoring.tcl ./scripts/
# cp /home/trudack/scripts/rosetta/rosetta_scoring.tcl ./analysis/
# cp /home/trudack/scripts/rosetta/align_rosetta.tcl ./scripts/
# cp /home/trudack/scripts/rosetta/cluster_fold.tcl  ./scripts/
# cp /home/trudack/scripts/rosetta/make_cluster_input.tcl ./scripts/
# cp /home/trudack/scripts/rosetta/make_cluster_input.tcl ./analysis/
# cp /home/trudack/scripts/rosetta/read_cluster_file.tcl ./scripts/
# cp /home/trudack/scripts/rosetta/read_cluster_file.tcl ./analysis/
# cp /home/trudack/scripts/rosetta/cluster_number.sh ./scripts/
# cp /home/trudack/scripts/rosetta/cluster_number.sh ./analysis/
# cp /home/trudack/scripts/rosetta/clustering.sh ./scripts/
# cp /home/trudack/scripts/rosetta/clustering.sh ./analysis/
# cp /home/trudack/scripts/rosetta/ss_analysis.tcl ./scripts/
# cp /home/trudack/scripts/rosetta/ss_analysis.tcl ./analysis/

#######################
# rosetta input files #
#######################

#fragment files and fasta
cp /Scr/scr-test-trudack/modelling_inputs/robetta_fragments/${MOL}/${MOL}* ./
cp /Scr/scr-test-trudack/modelling_inputs/fasta/${MOL}.fasta ./
# cp /home/trudack/scripts/rosetta/fold_termini.xml ./rosetta_input/
