import copy
import glob
import itertools
import os
# from config import PARIWISE_QA_SCRIPT,Q_SCORE,TM_SCORE_PATH
import sys

import numpy as np
import time
import config as config_path
import eva_utils
from eva_utils import specific_filename_reader
# from pdb_cleaning import pdb_filtering

is_monomer_scoring_done = False
is_dimer_scoring_done = False
PARIWISE_QA_SCRIPT = config_path.config.PARIWISE_QA_SCRIPT
TM_SCORE_PATH = config_path.config.TM_SCORE_PATH
Q_SCORE = config_path.config.Q_SCORE
DOCK_Q_PATH = config_path.config.DOCK_Q_PATH
MM_ALIGN = config_path.config.MM_ALIGN_PATH
GLINTER_DIR = config_path.config.GLINTER_DIR
total_time = time.perf_counter()
pdb_profile_dict = {}
eva_utils.check_path_exists(PARIWISE_QA_SCRIPT, TM_SCORE_PATH, Q_SCORE, DOCK_Q_PATH, MM_ALIGN, GLINTER_DIR)
import eva_utils as eva_util
import pdb_cleaning as pdb_c


#
# monomer_sequences_dir = "/home/bdmlab/MM_Eva/casp14/H1036/H1036.fasta"
# # # input_dir = "/home/bdmlab/hetero_test/lite_test/concatenated_pdb/"
# input_dir = "/home/bdmlab/MM_Eva/casp14/test_h1036_stuct/"
# stoichiometry = "A3B3C3"
# CPU_COUNT=10
# # # output_dir = "/home/bdmlab/hetero_test/lite_test/out/"
# output_dir = "/home/bdmlab/MM_Eva/casp14/test_h1036/"
# predicted_structures_AF2 = "/home/bdmlab/H1036/H1036_af2/"
#python multi_eva_modified.py ../data/fasta_casp14/casp_capri_fasta/T1032.fasta /home/rsr3gt/programs/Multi_Eva/Multimet_evatest_samples/predictions/T1032_lite/ A2 /home/rsr3gt/programs/Multi_Eva/data/pdbs_casp_alphafold/T1032/ /home/rsr3gt/programs/Multi_Eva/output/Qs_T1032/


#
# monomer_sequences_dir = "/home/bdmlab/MM_Eva/casp14/H1036/H1036.fasta"
# input_dir ="/home/bdmlab/MM_Eva/casp14/H1036/"
# stoichiometry = "A3B3C3"
# predicted_structures = "/home/bdmlab/MM_Eva/casp14/H1036/H1036_af2"
# output_dir = "/home/bdmlab/H1036/output_new/"
# CPU_COUNT=10
# predicted_structures_AF2 = "/home/bdmlab/H1036/H1036_af2/"
#
monomer_sequences_dir = sys.argv[1]
input_dir = sys.argv[2]
stoichiometry = sys.argv[3]
predicted_structures_AF2 = sys.argv[4]
CPU_COUNT = sys.argv[5].strip()
output_dir = sys.argv[6]

if os.path.isfile(monomer_sequences_dir):
    print(str(monomer_sequences_dir) + " Found")
if os.path.exists(input_dir):
    print(str(input_dir) + " Found")
print("stoichiometry " + str(stoichiometry))
if os.path.exists(predicted_structures_AF2):
    print(str(predicted_structures_AF2) + " Found")


#
#
# monomer_sequences_dir = "/home/bdmlab/hetero_test/H1045_new/h1045_dummt.fasta"
#
# # input_dir = "/home/bdmlab/hetero_test/lite_test/concatenated_pdb/"
# input_dir = "/home/bdmlab/hetero_test/H1045_new/concatenated_pdb/"
# stoichiometry = "A1B1"
# # output_dir = "/home/bdmlab/hetero_test/lite_test/out/"
# output_dir = "/home/bdmlab/hetero_test/H1045_new/output/"



TARGET_NAME = os.path.basename(monomer_sequences_dir).replace(".fasta", "")
fasta_dir = eva_utils.dir_maker(output_dir + "fasta/")
score_dir = eva_utils.dir_maker(output_dir + "score/")
monomer_score_dir = eva_utils.dir_maker(score_dir + "monomer/")
predicted_monomer_dir_AF2 = eva_utils.dir_maker(output_dir + "predicted_monomer/")
pred_structures = eva_util.specific_filename_reader(predicted_structures_AF2, "")
predicted_cmap_dir = eva_util.dir_maker(output_dir + "predicted_cmaps/")
fasta_stoic_dict = eva_util.multi_fasta_reader(_seq_file=monomer_sequences_dir)
eva_util.save_multi_fasta(fasta_dir, fasta_stoic_dict)
filtered_input_pdb = eva_util.dir_maker(output_dir + "filter_input_dir/")

pdb_c.pdb_filtering(_input_dir=input_dir,_output_dir=filtered_input_pdb)

warning_file = output_dir+"warning.log"
for values in pred_structures:
    # get_fasta_from_pdb_array
    temp_fasta_values = eva_utils.get_fasta_from_pdb_array(
        eva_utils.contents_to_info(eva_utils.read_pdb(predicted_structures_AF2 + values + ".pdb")))
    temp_seq = eva_util.sequence_finder(_seq_fasta_dict=fasta_stoic_dict, _fasta_string=temp_fasta_values)
    cmd = "cp " + predicted_structures_AF2 + str(values) + ".pdb " + str(predicted_monomer_dir_AF2) + "/sequence_" + str(
        temp_seq) + "_A.pdb"
    os.system(cmd)
    print(cmd)
    cmd = "cp " + predicted_structures_AF2 + str(values) + ".pdb " + str(predicted_monomer_dir_AF2) + "/sequence_" + str(
        temp_seq) + "_B.pdb"
    os.system(cmd)
    print(cmd)

stoi_details = eva_util.get_stoichiometry_details(stoichiometry)
#########
predicted_pdb_files = eva_util.specific_filename_reader(filtered_input_pdb, "")
TOTAL_SUBMISSION = len(predicted_pdb_files)
predicted_monomer_dir = eva_utils.dir_maker(output_dir + "monomer/")
all_chains_discovered = []

print("number of pdb is " + str(len(predicted_pdb_files)))
for monomer in predicted_pdb_files:
    temp_monomer_name = filtered_input_pdb + monomer
    for values in eva_util.monomer_pdb_filtering(_pdb=temp_monomer_name, _dir=predicted_monomer_dir):
        all_chains_discovered.append(values)

all_chains_discovered = list(dict.fromkeys(all_chains_discovered))

# print(all_chains_discovered)
all_monomer_files = []
predicted_monomer_chains_dir = eva_utils.dir_maker(output_dir + "monomer_chains/")
print(" Started seperating of Chains")
for pdb in predicted_pdb_files:
    temp_predicted_pdb_profile = eva_util.predicted_pdb_profile()
    temp_predicted_pdb_profile.monomers_chains = []
    all_pdb_skeleton = {}
    all_pdb_chain = {}
    # temp_cluster = {}
    temp_predicted_pdb_profile.name = pdb
    for monomer_chain in all_chains_discovered:
        a_monomer_chain_skeleton = predicted_monomer_dir + str(pdb) + "/" + str(
            pdb) + "_chain_" + monomer_chain + ".pdb"
        a_monomer_chain_fasta = predicted_monomer_dir + str(pdb) + "/" + str(pdb) + "_chain_" + monomer_chain + ".fasta"
        if os.path.exists(a_monomer_chain_skeleton):
            temp_predicted_pdb_profile.monomers_chains.append(copy.deepcopy(monomer_chain))
            monomer_chain_skeleton = eva_util.read_skeleton(a_monomer_chain_skeleton)
            all_pdb_skeleton[monomer_chain] = monomer_chain_skeleton
            all_pdb_chain[monomer_chain] = eva_util.read_fasta(a_monomer_chain_fasta)

    temp_predicted_pdb_profile.chain_skeleton_CA = all_pdb_skeleton
    temp_predicted_pdb_profile.chain_fasta = all_pdb_chain
    pdb_profile_dict[pdb] = temp_predicted_pdb_profile
print(" Ended seperating of Chains")
print("Multimer scoring started")
multimer_score_file = score_dir+"multimer_score.txt"
if not os.path.exists(multimer_score_file):
    for pdb_1 in predicted_pdb_files:
        # print(pdb_1)
        temp_MM_score_command = []
        mm_valie = 0
        for pdb_2 in predicted_pdb_files:
            if pdb_1 != pdb_2:
                # mm_valie = eva_util.get_MM_score(input_dir + "/" + pdb_1, input_dir + "/" + pdb_2, MM_ALIGN)
                temp_MM_score_command.append([filtered_input_pdb + "/" + pdb_1, filtered_input_pdb + "/" + pdb_2, MM_ALIGN])
        mm_valie =  eva_utils.get_MM_score_parallel_submit(temp_MM_score_command,CPU_COUNT)
        print(mm_valie)
        # print(str(np.average(temp_MM_score)))
        pdb_profile_dict.get(pdb_1).multimer_scoring =mm_valie

    eva_utils.save_mm_score(pdb_profile_dict,multimer_score_file)

print("Mapping of chains to clusters has been done")
#######NOW PUT ALL OF THEM IN ONE DIR OF MONOMER

for pdb in predicted_pdb_files:
    all_pdb_skeleton = {}
    for monomer_chain in pdb_profile_dict.get(pdb).monomers_chains:
        a_monomer_chain_skeleton = predicted_monomer_dir + str(pdb) + "/" + str(
            pdb) + "_chain_" + monomer_chain + ".pdb"
        if os.path.exists(a_monomer_chain_skeleton):
            monomer_chain_skeleton = eva_util.read_skeleton(a_monomer_chain_skeleton)
            all_pdb_skeleton[monomer_chain] = monomer_chain_skeleton
    pdb_profile_dict[pdb].chain_skeleton_CA = all_pdb_skeleton

all_true_sequence = fasta_stoic_dict.keys()

mm_score_dict = {}
if os.path.exists(multimer_score_file):
    mm_score_dict = eva_util.read_mm_score(_path=multimer_score_file)

for values in mm_score_dict:
    temp_pdb_profile = pdb_profile_dict.get(values)
    if mm_score_dict.get(values) != None:
        temp_pdb_profile.multimer_scoring = mm_score_dict.get(values)

    else:
        temp_pdb_profile.multimer_scoring = 0


print("Generating Final scoring part")
eva_util.print_final_data_mmalign(_file_name=output_dir + "/" + str(os.path.basename(monomer_sequences_dir).split(".")[0]) + ".csv",
                              _file_data=pdb_profile_dict,
                              _chain_data=fasta_stoic_dict)
#
# print("GRAND TOTAL TIME "+str(time.perf_counter()-total_time)+"\n")
