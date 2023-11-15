#!/usr/sfw/bin/python

#Author: Justyna Dopierala
#Date: Nov 2023

#################################################################################
### import packages ###
#################################################################################

import re
import os
import sys
import numpy as np
import pandas as pd

#################################################################################
### extract project path from command line arguments and check if exists ###
#################################################################################

dir_path = sys.argv[1]

# check if directory path provided as argument exists

isExist = os.path.exists(dir_path)
if not isExist:
  raise Exception("The path provided as command line argument does not exist")

#################################################################################
### define paths ###
#################################################################################

fa_file_path = dir_path + "data/fa_input/"
output_path = dir_path + "analysis/3_input_library_inspect/"

inspect_fa_out_file = open(output_path + "1_input_fa_INSPECT_output.txt", "w")

#################################################################################
### define re ###
#################################################################################

tab = re.compile("\t")
entry = re.compile(">")
NegCtrl = re.compile("NegCtrl")
NegCtrlIn = re.compile("NegCtrlintron")
PosCtrl = re.compile("PosCtrl")

#################################################################################
### lists, counters and hash tables ###
#################################################################################

gene_hash ={}
unique_gene_list = []
sequence_list = []
duplicates = []
wrong_seq_len_list = []
wrong_n_of_guides_gene_list = []

NegCtrl_counter = 0
NegCtrlIn_counter = 0
PosCtrl_counter = 0

#################################################################################
### open and inspect fa file ###
#################################################################################

fa_file = open(fa_file_path + "library.fa")           # open library file
check_points = [1000, 5000, 25000, 50000]

line_n = 0
for line in fa_file:
    line_n = line_n +1
    line = line.strip()
    
    if line_n in check_points:
      print("... Parsed line " + str(line_n))
    
    # find id line
    id_line  = entry.search(line)
    if id_line:
      line_split = line.split('|')                    # split line
      guide_id  = line_split[1]                       # extract guide id
      gene_name = line_split[2]                       # extract gene name
      
      pos_ctrl = PosCtrl.search(gene_name)
      if pos_ctrl:                                    # find and count positive control guides
        PosCtrl_counter = PosCtrl_counter + 1
      elif gene_name == "NegCtrl":                    # find and count negative contrtol guides
        NegCtrl_counter = NegCtrl_counter +1
      elif gene_name == "NegCtrlIntron":              # find and count negative control intron
        NegCtrlIn_counter = NegCtrlIn_counter + 1
      else:                                           # if not id line (i.e. sequence line)
        if gene_name not in gene_hash:                # put gene name in a hash table (key), if not there aready
          gene_hash[gene_name] = []
          gene_hash[gene_name].append(guide_id)       # append guide to a list of guides in a hash (as value)
          unique_gene_list.append(gene_name)          # append gene name to a list of unique genes
        else:
          gene_hash[gene_name].append(guide_id)
    else:                                             # if not sequence id line
      seq_length = len(line)
      if seq_length != 20:                            # check sequence length 
        wrong_seq_len_list.append(guide_id)
      if line not in sequence_list:                   # add sequnce to a list, if not there already
        sequence_list.append(line)
      else:                                           # find duplicate sequences
        duplicates.append(line)

fa_file.close()

#################################################################################
# check n of guides per gene
#################################################################################

for i in gene_hash:
  guide_list = gene_hash[i]
  n_of_guides = len(guide_list)
  if n_of_guides != 4:
    wrong_n_of_guides_gene_list.append(i)

# make a string of genes with wrong number of guides

genes_wrong_guide_n_string = ""

w_gene_count = 0
for gene in wrong_n_of_guides_gene_list:
  w_gene_count = w_gene_count + 1
  if w_gene_count ==1:
    genes_wrong_guide_n_string = genes_wrong_guide_n_string + gene
  else:
    genes_wrong_guide_n_string = genes_wrong_guide_n_string + ", " + gene

#################################################################################
# print output results to screen
#################################################################################

print("n of sequences in the library = " + str(len(sequence_list)))
print("n of duplicate sequences = " + str(len(duplicates)))
print("n of sequences with length other than 20 = " + str(len(wrong_seq_len_list)))

print("n of genes in library = " + str(len(unique_gene_list)))
print("NegCtrl_counter = " + str(NegCtrl_counter))
print("NegCtrlIn_counter = " + str(NegCtrlIn_counter))      
print("PosCtrl_counter = " + str(PosCtrl_counter))
print("genes where guide n other than 4: " + str(genes_wrong_guide_n_string))

#################################################################################
# write out results to file 
#################################################################################

inspect_fa_out_file.write("n of sequences in the library = " + str(len(sequence_list)) + "\n")
inspect_fa_out_file.write("n of duplicate sequences = " + str(len(duplicates)) + "\n")
inspect_fa_out_file.write("n of sequences with length other than 20 = " + str(len(wrong_seq_len_list)) + "\n")
inspect_fa_out_file.write("\n")
inspect_fa_out_file.write("N of genes in library = " + str(len(unique_gene_list)) + "\n")
inspect_fa_out_file.write("N of NegCtrl guides = " + str(NegCtrl_counter) + "\n")
inspect_fa_out_file.write("N of NegCtrlIn guides = " + str(NegCtrlIn_counter) + "\n")
inspect_fa_out_file.write("N of PosCtrl guides = " + str(PosCtrl_counter) + "\n")
inspect_fa_out_file.write("genes where guide n other than 4: " + genes_wrong_guide_n_string)

inspect_fa_out_file.close()

#################################################################################
### PART 2: load TCGA matrix and select genes ###
#################################################################################

full_TCGA_matrix_path = dir_path +"analysis/1_TCGA_matrix_full/TCGA_matrix_full.tsv"
full_TCGA_arr = open(full_TCGA_matrix_path)
selected_genes_TCGA_matrix_file = open(dir_path + "analysis/2_TCGA_matrix_selected/selected_genes_TCGA_matrix.tsv", "w")

genes_not_in_TCGA_matrix_list = []

line_counter = 0
gene_counter = 0
gene_mapped_counter = 0
gene_mapped_list = []

for line in full_TCGA_arr:
  line = line.strip()
  line_counter = line_counter + 1
  if line_counter == 1:
    selected_genes_TCGA_matrix_file.write(line + "\n")          # write header to output file
  else:
    gene_counter = gene_counter + 1
    line_split = tab.split(line)
    annot_part = line_split[0]
    annot_split = annot_part.split(" ")
    gene_name = annot_split[1].strip()
    if gene_name in unique_gene_list:                           # if gene in library.fa file write line to output file
      gene_mapped_counter = gene_mapped_counter + 1
      gene_mapped_list.append(gene_name)
      selected_genes_TCGA_matrix_file.write(line + "\n")
    else:
      genes_not_in_TCGA_matrix_list.append(gene_name)
  
selected_genes_TCGA_matrix_file.close()
full_TCGA_arr.close()

# get the list of genes in the library that don't have a match in TCGA matrix

NOT_mapped_genes = []

for i in unique_gene_list:
  if i not in gene_mapped_list:
    NOT_mapped_genes.append(i)

print(len(NOT_mapped_genes))
print(NOT_mapped_genes[1:10])


