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

is_monomer_scoring_done = True
PARIWISE_QA_SCRIPT = config_path.config.PARIWISE_QA_SCRIPT
TM_SCORE_PATH = config_path.config.TM_SCORE_PATH
Q_SCORE = config_path.config.Q_SCORE
total_time = time.perf_counter()
pdb_profile_dict = {}
import eva_utils as eva_util

monomer_sequences_dir = "/home/bdmlab/multimer_test/input/H1036.fasta"
# input_dir = "/home/bdmlab/hetero_test/lite_test/concatenated_pdb/"
input_dir = "/home/bdmlab/multimer_test/input/H1036/"
stoichiometry = "A3B3C3"
# output_dir = "/home/bdmlab/hetero_test/lite_test/out/"
output_dir = "/home/bdmlab/multimer_test/output/"

# monomer_sequences_dir = sys.argv[1]
# # input_dir = "/home/bdmlab/hetero_test/lite_test/concatenated_pdb/"
# input_dir =sys.argv[2]
# stoichiometry = sys.argv[3]
# # output_dir = "/home/bdmlab/hetero_test/lite_test/out/"
# output_dir = sys.argv[4]


TARGET_NAME = os.path.basename(monomer_sequences_dir).replace(".fasta", "")
fasta_dir = eva_utils.dir_maker(output_dir + "fasta/")
score_dir = eva_utils.dir_maker(output_dir + "score/")
monomer_score_dir = eva_utils.dir_maker(score_dir + "monomer/")

fasta_stoic_dict = eva_util.multi_fasta_reader(_seq_file=monomer_sequences_dir)
eva_util.save_multi_fasta(fasta_dir, fasta_stoic_dict)

stoi_details = eva_util.get_stoichiometry_details(stoichiometry)
predicted_pdb_files = eva_util.specific_filename_reader(input_dir, "")
TOTAL_SUBMISSION = len(predicted_pdb_files)
predicted_monomer_dir = eva_utils.dir_maker(output_dir + "monomer/")
all_chains_discovered = []
start = time.perf_counter()
for monomer in predicted_pdb_files:
    temp_monomer_name = input_dir + monomer
    for values in eva_util.monomer_pdb_filtering(_pdb=temp_monomer_name, _dir=predicted_monomer_dir):
        all_chains_discovered.append(values)

all_chains_discovered = list(dict.fromkeys(all_chains_discovered))

# print(all_chains_discovered)
all_monomer_files = []
predicted_monomer_chains_dir = eva_utils.dir_maker(output_dir + "monomer_chains/")
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

#######chain cluster mapper
for pdb in pdb_profile_dict:
    temp_predicted_pdb_profile = copy.deepcopy(pdb_profile_dict.get(pdb))
    temp_chain_cluster = {}
    temp_cluster_chain = {}
    for values in temp_predicted_pdb_profile.chain_fasta:
        temp = eva_util.sequence_finder(fasta_stoic_dict, temp_predicted_pdb_profile.chain_fasta.get(values))
        temp_chain_cluster[values] = temp
        if temp_cluster_chain.get(temp) != None:
            temp_cluster_chain[temp].append(values)
            # update
        else:
            temp_cluster_chain[temp] = [values]

    pdb_profile_dict.get(pdb).chain_cluster = copy.deepcopy(temp_chain_cluster)
    pdb_profile_dict.get(pdb).cluster_chain = copy.deepcopy(temp_cluster_chain)
#######NOW PUT ALL OF THEM IN ONE DIR OF MONOMER
for true_squence in fasta_stoic_dict:
    current_dir_name = predicted_monomer_chains_dir + "/sequence_" + str(true_squence) + "/"
    eva_util.dir_maker(current_dir_name)
    for _pdb in pdb_profile_dict:
        temp_pdb = pdb_profile_dict.get(_pdb).cluster_chain.get(true_squence)
        for chains in temp_pdb:
            monomer_pdb_name = predicted_monomer_dir + str(_pdb) + "/" + str(_pdb) + "_chain_" + str(chains) + ".pdb"
            if os.path.exists(monomer_pdb_name):
                os.system("cp " + monomer_pdb_name + " " + current_dir_name)

if not is_monomer_scoring_done:
    ##################### MONOMER SCORING PART #################################
    for chain_value in fasta_stoic_dict:
        temp_chain_dir = predicted_monomer_chains_dir + "sequence_" + str(chain_value) + "/"
        # print(predicted_monomer_dir + "/**/*"+"_chain_"+str(chain_value))
        all_monomer_chained_files = glob.glob(predicted_monomer_dir + "/**/*" + "_chain_" + str(chain_value) + ".pdb",
                                              recursive=True)
        cmd = "perl " + PARIWISE_QA_SCRIPT + " " + temp_chain_dir + " " + fasta_dir + "sequence_" + str(
            chain_value) + "_A.fasta" + " " + Q_SCORE + " " + TM_SCORE_PATH + " " + chain_value + " " + monomer_score_dir
        print(cmd)
        os.system(cmd)
    #################### MONOMER SCORING PART #################################
    end_time_start = time.perf_counter()
    print("qA score time " + str(end_time_start - start) + "\n")

# ### WARNING
# # removing unmapped files
# modified_all_chains_discovered= copy.deepcopy(all_chains_discovered)
# for values in modified_all_chains_discovered :
#     if  a_multimer.chain_fasta.get(values) == None:
#         all_chains_discovered.remove(values)
#
#
# start = time.perf_counter()
#
#
# # ALL CHAINS CA
# start = time.perf_counter()
for pdb in predicted_pdb_files:
    # temp_predicted_pdb_profile = eva_util.predicted_pdb_profile()
    # temp_predicted_pdb_profile.monomers_chains = []
    all_pdb_skeleton = {}
    # temp_predicted_pdb_profile.name = pdb
    for monomer_chain in pdb_profile_dict.get(pdb).monomers_chains:
        a_monomer_chain_skeleton = predicted_monomer_dir + str(pdb) + "/" + str(
            pdb) + "_chain_" + monomer_chain + ".pdb"
        if os.path.exists(a_monomer_chain_skeleton):
            # temp_predicted_pdb_profile.monomers_chains.append(copy.deepcopy(monomer_chain))
            monomer_chain_skeleton = eva_util.read_skeleton(a_monomer_chain_skeleton)
            all_pdb_skeleton[monomer_chain] = monomer_chain_skeleton
    # temp_predicted_pdb_profile.chain_skeleton_CA = all_pdb_skeleton
    pdb_profile_dict[pdb].chain_skeleton_CA = all_pdb_skeleton

all_true_sequence = fasta_stoic_dict.keys()

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
            len_a_dimers = len(temp.cluster_chain.get(dimers[0]))
            if len_a_dimers == 1:
                temp_all_dimer_combination.remove(dimers)

    for dimers in temp_all_dimer_combination:
        first_chain_list = temp.cluster_chain.get(dimers[0])
        second_chain_list = temp.cluster_chain.get(dimers[1])
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

print("here")
predicted_structures_dimer_cmap_dir = eva_util.dir_maker(output_dir + "struct_dimer_cmaps/")
dimer_strcutures_dir = eva_util.dir_maker(output_dir + "dimer_structures_pdb/")
# eva_util.dir_maker(dimer_strcutures_dir + str("all") + "/")
# eva_util.dir_maker(predicted_structures_dimer_cmap_dir + str("all") + "/")
for values in valid_dimer_combos:
    eva_util.dir_maker(dimer_strcutures_dir + str("sequence_") + values + "/")
    eva_util.dir_maker(predicted_structures_dimer_cmap_dir + str("sequence_") + values + "/")

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
            temp_dimer_array = eva_util.chain_replacer(copy.deepcopy(monomer_a), "A") + eva_util.chain_replacer(copy.deepcopy(monomer_b), "B")
            eva_util.pdb_from_array(_pdb=temp_dimer_array, _filename=dimer_pdb_file_name)

#             ####THEN CONCAT DIME


# ##introduce that 20%
# #### 5 #### GENERATE CMAPS
#
#
# print(valid_dimer_combindations)
# ############################ TBD ##########################
# #### 4 #### GLINTER RUN
# #### FOR NOW PRECOMPUTED
# ############################ TBD ##########################
#
#
#
# ################QUESION ALPHAFOLD STRUCTURE ??
# #then run glinter
# ####################ALREADY DONE CURRENTLY USING PRE_COMPUTED
#####target  wise
for pdb in pdb_profile_dict:
    print(pdb)
    temp_pdb_profile = pdb_profile_dict.get(pdb)
    a_dimer_score_dict = {}
    #####possible dimer wise

    for values in valid_dimer_combos:
        print(values)

        # all_specific_dimer =glob.glob(dimer_strcutures_dir + "sequence_" + str(values) + "/*" + ".pdb",
        #           recursive=True)
        all_specific_dimer = eva_util.specific_filename_reader(  _input_dir=dimer_strcutures_dir + "sequence_" + str(values),            _extension="pdb")
        all_specific_target_dimer = eva_util.specific_filename_reader(
            _input_dir=dimer_strcutures_dir + "sequence_" + str(values), _extension=pdb)
        temp_same_dimer_wise = []
        ##### same type of dimer wise
        for first_chain in all_specific_target_dimer:

            temp_dimer_chain_wise = []
            chain_first = dimer_strcutures_dir + "sequence_" + str(values) + "/" + str(first_chain) + ".pdb"
            print(chain_first)
            #####specific dimer wise
            for predicted_dimer in all_specific_dimer:
                chain_second = dimer_strcutures_dir + "sequence_" + str(values) + "/" + predicted_dimer + ".pdb"
                print(chain_second)
                if chain_first != chain_second:
                    temp_dimer_chain_wise.append(eva_util.get_dock_q_score(_true=chain_first, _current=chain_second))
            temp_same_dimer_wise.append(np.average(temp_dimer_chain_wise))
        if len(temp_same_dimer_wise) != 0:
            a_dimer_score_dict[values] = np.max(temp_same_dimer_wise)
        else:
            a_dimer_score_dict[values] = 0
    # else:
    #     a_dimer_score_dict[values] = 0.0
    pdb_profile_dict.get(pdb).ds_scores = a_dimer_score_dict
print("here")
# #load glinter cmap
# end_time_start = time.perf_counter()
# print("DS SCORING PART "+str(end_time_start-start)+"\n")
# print("RECALL CALULCATOR")
# start = time.perf_counter()
# # print("DS SCORING PART "+str(end_time_start-start)+"\n")
# predicted_cmap_dir= output_dir +"predicted_cmaps/"
# #
# # #get icps score
# for pdb in pdb_profile_dict:
#     print(pdb)
#     temp_pdb_profile = pdb_profile_dict.get(pdb)
#
#     for dimer in valid_dimer_combindations:
#
#         if dimer in temp_pdb_profile.dimers:
#             dimer_cmap_file_name = predicted_structures_dimer_cmap_dir+str(pdb)+"_chain_"+dimer+".cmap"
#             #GLINTER FORMAT
#             pred_file_name = predicted_cmap_dir+TARGET_NAME+"_"+dimer[0]+":"+TARGET_NAME+"_"+dimer[1]+".cmap"
#             prev_icps =copy.deepcopy( pdb_profile_dict.get(pdb).icps_scores)
#             prev_recall =copy.deepcopy( pdb_profile_dict.get(pdb).recall)
#             if os.path.exists(pred_file_name) and os.path.exists(dimer_cmap_file_name):
#                 a_icps_score=copy.deepcopy(eva_util.get_icps_score(_struct_cmap=dimer_cmap_file_name,_pred_cmap=pred_file_name))
#                 if len(prev_icps)>0:
#                     prev_icps.update({dimer,a_icps_score})
#                 else:
#                     prev_icps[dimer] = a_icps_score
#                 pdb_profile_dict.get(pdb).icps_scores = prev_icps
#
#
#                 a_recall_score = copy.deepcopy( eva_util.get_recall(_struct_cmap=dimer_cmap_file_name,_pred_cmap=pred_file_name))
#                 if len(prev_recall)> 0:
#                     prev_recall.update({dimer,a_recall_score})
#                 else:
#                     prev_recall[dimer] = a_recall_score
#
#                 pdb_profile_dict.get(pdb).recall = prev_recall
#         else:
#             prev_icps = copy.deepcopy(pdb_profile_dict.get(pdb).icps_scores)
#             prev_recall = copy.deepcopy(pdb_profile_dict.get(pdb).recall)
#             if len(prev_icps) > 0:
#                 prev_icps.update({dimer, 0})
#             else:
#                 prev_icps[dimer] = 0
#             pdb_profile_dict.get(pdb).icps_scores = prev_icps
#             if len(prev_recall) > 0:
#                 prev_recall.update({dimer, 0})
#             else:
#                 prev_recall[dimer] = 0
#
#             pdb_profile_dict.get(pdb).icps_scores = prev_icps
#             pdb_profile_dict.get(pdb).recall = prev_recall
#
# monomer_score_dict ={}
# end_time_start = time.perf_counter()
# print("recall and icps SCORING PART "+str(end_time_start-start)+"\n")
# print("MONOMER SCORE  CALULCATOR")
# #
#
# for monomers in all_chains_discovered:
#     monomer_score_file =  monomer_score_dir+"/"+monomers+str(".tm")
#     monomer_score_dict[monomers] = eva_util.read_monomer_score(_path = monomer_score_file)
# print("FINAL SCORE  CALULCATOR")
# for pdb_values in pdb_profile_dict:
#     temp_pdb_profile = pdb_profile_dict.get(pdb_values)
#     prev_ms_score = copy.deepcopy(pdb_profile_dict.get(pdb_values).ms_scores)
#     mono_score_dict = {}
#     for values in all_chains_discovered:
#         a_monomer_score  =  monomer_score_dict.get(values).get(pdb_values)
#         mono_score_dict[values]=a_monomer_score
#
#     pdb_profile_dict.get(pdb_values).ms_scores = mono_score_dict
#
# print("here")
#
# eva_util.print_final_data(_file_name="/home/bdmlab/result_time_test.csv",_file_data=pdb_profile_dict,_chain_data =all_chains_discovered)
#
# print("GRAND TOTAL TIME "+str(time.perf_counter()-total_time)+"\n")
