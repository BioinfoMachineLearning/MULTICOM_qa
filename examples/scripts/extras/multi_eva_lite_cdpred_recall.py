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
CDPRED_DIR = config_path.config.CDPRED_PATH
total_time = time.perf_counter()
pdb_profile_dict = {}
eva_utils.check_path_exists(PARIWISE_QA_SCRIPT, TM_SCORE_PATH, Q_SCORE, DOCK_Q_PATH, MM_ALIGN, GLINTER_DIR)
import eva_utils as eva_util
import pdb_cleaning as pdb_c


#
# monomer_sequences_dir = "/media/bdmlab/Data/DATA_RAJ/QA_PAPER/CASP15/CASP_fasta/H1129.fasta"
# input_dir ="/media/bdmlab/Data/DATA_RAJ/QA_PAPER/CASP15/input/H1129//"
# stoichiometry = "A1B1"
# predicted_structures_AF2 = "/media/bdmlab/Data/DATA_RAJ/QA_PAPER/CASP15/CASP_AF2_ts/H1129/"
# CPU_COUNT=10
# output_dir = "/media/bdmlab/Data/DATA_RAJ/QA_PAPER/CASP15/recall_output/H1129/"

# predicted_structures_AF2 = "/home/bdmlab/MM_Eva/casp14/H1036/H1036_af2/"
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

print("Multimer scoring Done")
print("Mapping chains to clusters")
#######chain cluster mapper
for pdb in pdb_profile_dict:
    # print(pdb)
    temp_predicted_pdb_profile = copy.deepcopy(pdb_profile_dict.get(pdb))
    temp_chain_cluster = {}
    temp_cluster_chain = {}
    for values in temp_predicted_pdb_profile.chain_fasta:
        temp = eva_util.sequence_finder(fasta_stoic_dict, temp_predicted_pdb_profile.chain_fasta.get(values))
        temp_chain_cluster[values] = str(temp)
        if temp_cluster_chain.get(str(temp)) != None:
            temp_cluster_chain[str(temp)].append(values)
            # update
        else:
            temp_cluster_chain[str(temp)] = [values]

    pdb_profile_dict.get(pdb).chain_cluster = copy.deepcopy(temp_chain_cluster)
    pdb_profile_dict.get(pdb).cluster_chain = copy.deepcopy(temp_cluster_chain)
print("Mapping of chains to clusters has been done")

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
print("Monomer scoring Ended")
#
#
# #### 1 #### FIND COMBINATION of all PAIRS
# we will get all heterodimer connections here
all_dimer_combination = list(itertools.combinations(all_true_sequence, 2))
# all_dimer_combination = filter(lambda x: (x[0] == x[1]), all_dimer_combination)
##############ADJUSTING FOR HOMODIMERS
for dimers in fasta_stoic_dict:
    all_dimer_combination.append([dimers, dimers])

### check for homomdimer
### check for heterodimer
all_pdb_dimers_contact = []
# #### 2 #### FIND If CONTACT
for pdb in pdb_profile_dict:
    temp = pdb_profile_dict.get(pdb)
    temp_dimer = []
    temp_all_dimer_combination = copy.deepcopy(all_dimer_combination)
    for dimers in all_dimer_combination:
        if dimers[0] == dimers[1]:
            if temp.cluster_chain.get(str(dimers[0])) != None:
                len_a_dimers = len(temp.cluster_chain.get(str(dimers[0])))
            
                if len_a_dimers == 1:
                    temp_all_dimer_combination.remove(dimers)

    for dimers in temp_all_dimer_combination:
        first_chain_list = temp.cluster_chain.get(str(dimers[0]))
        second_chain_list = temp.cluster_chain.get(str(dimers[1]))

        if first_chain_list != None and second_chain_list !=None:
            for first_chain in first_chain_list:
                for second_chain in second_chain_list:
                    if first_chain != second_chain:
                        first_chain_pdb = temp.chain_skeleton_CA[first_chain]
                        second_chain_pdb = temp.chain_skeleton_CA[second_chain]
                        if eva_util.if_contact(first_chain_pdb, second_chain_pdb):
                            temp_dimer.append(str(first_chain) + str(second_chain))
                            all_pdb_dimers_contact.append(str(pdb) + "_" + str(dimers[0]) + str(dimers[1]))
    rr_removed_temp_dimer = []
    for values in temp_dimer:
        temp_values = []
        for char in values:
            temp_values.append(char)

        temp_values.sort()
        # print(str(temp_values[0])+str(temp_values[1]))
        rr_removed_temp_dimer.append(str(temp_values[0]) + str(temp_values[1]))
    rr_removed_temp_dimer = list(dict.fromkeys(rr_removed_temp_dimer))
    temp.dimers = copy.deepcopy(rr_removed_temp_dimer)

#################chain_hit.count('-')
all_pdb_dimers_contact = list(dict.fromkeys(all_pdb_dimers_contact))
valid_dimer_in_contact = []
print("Finding the valide dimers")
for values in all_pdb_dimers_contact:
    temp_arr = values.split("_")
    valid_dimer_in_contact.append(temp_arr[len(temp_arr) - 1])

valid_dimer_combos = []
for _pdbs in all_dimer_combination:
    str_pdb = str(_pdbs[0]) + str(_pdbs[1])
    number = valid_dimer_in_contact.count(str_pdb)
    dimer_treshold_consideration = int(len(pdb_profile_dict) * 0.2)
    if number >= dimer_treshold_consideration:
        valid_dimer_combos.append(str_pdb)
print("Completed looking for the valide dimer")


predicted_structures_dimer_cmap_dir = eva_util.dir_maker(output_dir + "struct_dimer_cmaps/")
dimer_strcutures_dir = eva_util.dir_maker(output_dir + "dimer_structures_pdb/")

for values in valid_dimer_combos:
    eva_util.dir_maker(dimer_strcutures_dir + str("sequence_") + values + "/")
    eva_util.dir_maker(predicted_structures_dimer_cmap_dir + str("sequence_") + values + "/")
print("Generating the cmpas of the predicted structures")
for pdb in pdb_profile_dict:
    temp_pdb_profile = pdb_profile_dict.get(pdb)
    dimer_for_cmaps = eva_util.dimer_for_cmaps(valid_dimer_combos, temp_pdb_profile)
    for dimer in dimer_for_cmaps:
        cluster_1 = temp_pdb_profile.chain_cluster.get(dimer[0])
        cluster_2 = temp_pdb_profile.chain_cluster.get(dimer[1])
        chain_concat = cluster_1 + cluster_2
        if int(cluster_1) > int(cluster_2):
            chain_concat = cluster_2 + cluster_1

        monomer_a = temp_pdb_profile.chain_skeleton_CA.get(dimer[0])
        monomer_b = temp_pdb_profile.chain_skeleton_CA.get(dimer[1])
        dimer_cmap_file_name = predicted_structures_dimer_cmap_dir + "sequence_" + str(chain_concat) + "/" + str(
            pdb) + "_chain_" + dimer + ".cmap"

        dimer_pdb_file_name = dimer_strcutures_dir + "sequence_" + str(chain_concat) + "/" + str(
            pdb) + "_chain_" + dimer + ".pdb"
        if not os.path.exists(dimer_cmap_file_name):
            temp_cmaps = eva_util.get_CA_cmaps(_first_chain=monomer_a, _second_chain=monomer_b)
            np.savetxt(dimer_cmap_file_name, temp_cmaps)
        if not os.path.exists(dimer_pdb_file_name):
            temp_dimer_array = eva_util.chain_replacer(copy.deepcopy(monomer_a), "A") + eva_util.chain_replacer(
                copy.deepcopy(monomer_b), "B")
            eva_util.pdb_from_array(_pdb=temp_dimer_array, _filename=dimer_pdb_file_name)

print("Concluded generating the cmpas of the predicted structures")


print("generating the cmpas using the predictors")


extra_cmaps = eva_util.dir_maker(predicted_cmap_dir + "/extras/")
for _pred_cmap_candidate in valid_dimer_combos:
    _is_homodimer = False
    first_chain = ""
    second_chain = ""
    expected_cmaps_name = predicted_cmap_dir
    if _pred_cmap_candidate[0] == _pred_cmap_candidate[1]:
        _is_homodimer = True
        first_chain = predicted_monomer_dir_AF2 + "/sequence_" + _pred_cmap_candidate[0] + "_A.pdb"
        second_chain = predicted_monomer_dir_AF2 + "/sequence_" + _pred_cmap_candidate[0] + "_B.pdb"
        expected_cmaps_name = expected_cmaps_name + "/sequence_" + str(_pred_cmap_candidate[0]) + "_A:sequence_" + str(
            _pred_cmap_candidate[1]) + "_B.cmap"
    else:
        first_chain = predicted_monomer_dir_AF2 + "/sequence_" + _pred_cmap_candidate[0] + "_A.pdb"
        second_chain = predicted_monomer_dir_AF2 + "/sequence_" + _pred_cmap_candidate[1] + "_A.pdb"
        expected_cmaps_name = expected_cmaps_name + "/sequence_" + str(_pred_cmap_candidate[0]) + "_A:sequence_" + str(
            _pred_cmap_candidate[1]) + "_A.cmap"

    if not os.path.exists(expected_cmaps_name):
        print("expected_cmaps " + str(expected_cmaps_name) + " not found \n")
        print("predicting the interactions \n")
        eva_utils.cdpred_runner(first_chain, second_chain, extra_cmaps, _is_homodimer, expected_cmaps_name,CDPRED_DIR)
    else:
        print("expected_cmaps " + str(expected_cmaps_name) + " found \n")


print("Concluded generating the cmpas using the predictors")


print("Started icps and recall scoring part")
for pdb in pdb_profile_dict:
    print(pdb)
    temp_pdb_profile = pdb_profile_dict.get(pdb)
    # temp_icps_target = {}
    temp_recall_target = {}
    for dimer in valid_dimer_combos:
        ## get all the dimer cmaps
        ##for each calculate check if empty and also take the maximum
        dimer_cmap_dir = predicted_structures_dimer_cmap_dir + "sequence_" + str(dimer) + "/"
        print(dimer_cmap_dir)
        if dimer[0] != dimer[1]:
            predicted_cmap = predicted_cmap_dir + "sequence_" + str(dimer[0]) + "_A:" + "sequence_" + str(
                dimer[1]) + "_A.cmap"
        else:
            predicted_cmap = predicted_cmap_dir + "sequence_" + str(dimer[0]) + "_A:" + "sequence_" + str(
                dimer[1]) + "_B.cmap"
        ## get all the cmaps
        all_cmap_same_type = eva_util.specific_filename_reader(dimer_cmap_dir, pdb)
        # temp_list_icps = [0]
        temp_list_recall = [0]
        for _cmaps in all_cmap_same_type:
            dimer_cmap_file_name = dimer_cmap_dir + str(_cmaps) + ".cmap"
            if os.path.exists(predicted_cmap) and os.path.exists(dimer_cmap_file_name):
                transpose = False
                chains = _cmaps.split("_")[-1]
                temp_transpose_check = str(temp_pdb_profile.chain_cluster.get(chains[0])) + str(
                    temp_pdb_profile.chain_cluster.get(chains[1]))
                if temp_transpose_check != dimer:
                    transpose = True
                possible = True
                # try:
                a_recall_score = copy.deepcopy( eva_util.get_recall(_struct_cmap=dimer_cmap_file_name, _pred_cmap=predicted_cmap,
                                            _transpose=transpose))
                # except:
                # a_recall_score = 0

                temp_list_recall.append(a_recall_score)

        if len(temp_recall_target) > 0:
            temp_recall_target.update({dimer: np.max(temp_list_recall)})
        else:
            temp_recall_target[dimer] = np.max(temp_list_recall)

    pdb_profile_dict.get(pdb).recall = temp_recall_target
print("Completed icps and recall scoring part")
monomer_score_dict = {}

mm_score_dict = {}
if os.path.exists(multimer_score_file):
    mm_score_dict = eva_util.read_mm_score(_path=multimer_score_file)

for values in mm_score_dict:
    temp_pdb_profile = pdb_profile_dict.get(values)
    if mm_score_dict.get(values) != None:
        temp_pdb_profile.multimer_scoring = mm_score_dict.get(values)

    else:
        temp_pdb_profile.multimer_scoring = 0


for pdb_values in pdb_profile_dict:
    temp_pdb_profile = pdb_profile_dict.get(pdb_values)
    # prev_ms_score = copy.deepcopy(pdb_profile_dict.get(pdb_values).ms_scores)
    temp_mono_score_dict = {}
    for values in fasta_stoic_dict:
        temp_ms_chainwise = []
        # a_monomer_score = monomer_score_dict.get(str(values)).get(pdb_values)
        monomer_scores_targetwise = monomer_score_dict.get(str(values))
        if monomer_scores_targetwise != None:
            for ms in monomer_scores_targetwise:
                if pdb_values in ms:
                    temp_ms_chainwise.append(monomer_scores_targetwise.get(ms))
        print(temp_ms_chainwise)
        print(pdb_values)
        if len(temp_ms_chainwise)>0:
            temp_mono_score_dict[values] = np.max(temp_ms_chainwise)
        else:
             temp_mono_score_dict[values] = 0

    pdb_profile_dict.get(pdb_values).ms_scores = temp_mono_score_dict

print("Generating Final scoring part")
eva_util.print_final_data_new_lite_recall(_file_name=output_dir + "/" + str(os.path.basename(monomer_sequences_dir).split(".")[0]) + "_recall.csv",
                              _file_data=pdb_profile_dict,
                              _chain_data=fasta_stoic_dict, _dimer_data=valid_dimer_combos)
#
# print("GRAND TOTAL TIME "+str(time.perf_counter()-total_time)+"\n")
