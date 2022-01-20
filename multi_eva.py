import copy
import glob
import itertools
import os
# from config import PARIWISE_QA_SCRIPT,Q_SCORE,TM_SCORE_PATH
import numpy as np

import  config as config_path
PARIWISE_QA_SCRIPT = config_path.config.PARIWISE_QA_SCRIPT
TM_SCORE_PATH=config_path.config.TM_SCORE_PATH
Q_SCORE = config_path.config.Q_SCORE
# print(config_path.config.PARIWISE_QA_SCRIPT)

pdb_profile_dict  = {}
######REINDEXING SHOULD BE DONE


import eva_utils as eva_util
# monomer_sequences_dir = "/home/bdmlab/T1032o.fasta"
# input_dir = "/home/bdmlab/Multimet_evatest_samples/predictions/T1032o_lite/"
# stoichiometry = "A2"
# output_dir = "/home/bdmlab/multi_eva_test/"

monomer_sequences_dir = "/home/bdmlab/Multimet_evatest_samples/casp_fasta/T1038.fasta"
input_dir = "/home/bdmlab/Multimet_evatest_samples/predictions/T1038o_lite/"
stoichiometry = "A2"
output_dir = "/home/bdmlab/multi_eva_test/T1038_LITE/"


score_dir = output_dir+"score/"
os.system("mkdir -p "+score_dir)
monomer_score_dir = score_dir+"monomer"
os.system("mkdir -p "+monomer_score_dir)

stoi_details = eva_util.get_stoichiometry_details(stoichiometry)
predicted_pdb_files = eva_util.specific_filename_reader(input_dir,"")
# make the monomer dir
predicted_monomer_dir = output_dir+"monomer/"
os.system("mkdir -p "+predicted_monomer_dir)
all_chains_discovered= []

#get all the chains
for monomer in predicted_pdb_files :
    temp_monomer_name = input_dir+monomer
    for values in eva_util.monomer_pdb_filtering(_pdb= temp_monomer_name,_dir = predicted_monomer_dir):
        all_chains_discovered.append(values)

all_chains_discovered = list(dict.fromkeys(all_chains_discovered))
# print(all_chains_discovered)
all_monomer_files = []
predicted_monomer_chains_dir = output_dir+"monomer_chains/"
os.system("mkdir -p " + predicted_monomer_chains_dir)
##Monomer Finding
#copying all the files of the same type on one folder
# a_multimer = eva_util.multimer
a_multimer = eva_util.fasta_to_chain_mapper(_fasta_file=monomer_sequences_dir, _stoi=stoichiometry, _chains=copy.deepcopy(all_chains_discovered))
fasta_dir = output_dir+"fasta/"
os.system("mkdir -p "+fasta_dir)
##lets make fasta dir and put fasta files there
for chain_name  in a_multimer.chains:
    print(chain_name)
    print(a_multimer.chain_fasta[chain_name])
    fasta_content= ">"+str(chain_name)+"\n"+a_multimer.chain_fasta[chain_name]
    eva_util.write2File(_filename=fasta_dir+str(chain_name)+".fasta",_cont=fasta_content)

##fasta_to_chain_mapper(_fasta_file="/home/bdmlab/T1032o.fasta", _stoi="A2", _chains=["A", "B"]):

#	die "need six parameters: input model dir, fasta sequence file, pairwise_QA path, tm_score path, target name (output name), output dir.\n";
# for chain_values in all_chains_discovered
##################### MONOMER SCORING PART #################################
for chain_value in all_chains_discovered:
    temp_chain_dir =  predicted_monomer_chains_dir+str(chain_value)
    os.system("mkdir -p " +temp_chain_dir)
    print(predicted_monomer_dir + "/**/*"+"_chain_"+str(chain_value))
    all_monomer_chained_files = glob.glob(predicted_monomer_dir + "/**/*"+"_chain_"+str(chain_value)+".pdb", recursive = True)
    for values in all_monomer_chained_files:
        os.system("cp "+values+" "+temp_chain_dir)

    # command for pairwise_qa
    # perl pairwise_model_eva.pl /home/rajroy/multi_eva_test/monomer_chains/A/ /home/rajroy/multi_eva_test/Multimet_evatest_samples/casp_fasta/H1036A.fasta /home/rajroy/pairwiseQA/q_score /home/rajroy/Downloads/tools/TMscore A /home/rajroy/q_A/

    # if chain_value !="A":
    # cmd = PARIWISE_QA_SCRIPT +" "+predicted_monomer_chains_dir+str(chain_value)+" "+monomer_sequences_dir+ " "+Q_SCORE+ " "+TM_SCORE_PATH+" "+chain_value + " "+output_dir
    cmd = PARIWISE_QA_SCRIPT +" "+predicted_monomer_chains_dir+str(chain_value)+" "+fasta_dir+str(chain_value)+".fasta"+ " "+Q_SCORE+ " "+TM_SCORE_PATH+" "+chain_value + " "+monomer_score_dir
    print(cmd)
    os.system(cmd)
#################### MONOMER SCORING PART #################################

# ALL CHAINS CA
for pdb in predicted_pdb_files:
    temp_predicted_pdb_profile = eva_util.predicted_pdb_profile()
    temp_predicted_pdb_profile.monomers_chains =[]
    all_pdb_skeleton = {}
    temp_predicted_pdb_profile.name = pdb
    for monomer_chain in all_chains_discovered:
        a_monomer_chain_skeleton = predicted_monomer_dir + str(pdb) + "/" + str(pdb) + "_chain_" + monomer_chain + ".pdb"
        if os.path.exists(a_monomer_chain_skeleton):
            temp_predicted_pdb_profile.monomers_chains.append(copy.deepcopy(monomer_chain))
            monomer_chain_skeleton = eva_util.read_skeleton(a_monomer_chain_skeleton)
            all_pdb_skeleton[monomer_chain]=monomer_chain_skeleton

    temp_predicted_pdb_profile.chain_skeleton_CA = all_pdb_skeleton
    pdb_profile_dict[pdb] = temp_predicted_pdb_profile


#### 1 #### FIND COMBINATION of all PAIRS
all_dimer_combination = list (itertools.combinations(all_chains_discovered, 2))
#remove duplicate
# all_dimer_combination = filter(lambda x: (x[0] == x[1]), all_dimer_combination)

all_pdb_dimers_contact = []
#### 2 #### FIND If CONTACT
for pdb in pdb_profile_dict:
    temp =     pdb_profile_dict.get(pdb)
    temp_dimer = []
    for dimers in all_dimer_combination:
        first_chain = dimers[0]
        first_chain_pdb =temp.chain_skeleton_CA[first_chain]
        second_chain = dimers[1]
        second_chain_pdb = temp.chain_skeleton_CA[second_chain]
        if eva_util.if_contact(first_chain_pdb,second_chain_pdb) ==True:
            temp_dimer.append(str(first_chain)+str(second_chain))
            all_pdb_dimers_contact.append(str(first_chain)+str(second_chain))

    temp.dimers=copy.deepcopy(temp_dimer)
#find that 20% threshold for dimers in contact
#### 3 #### FIND 20% overlapping pairs
valid_dimer_combindations = []
dimer_treshold_consideration = int(len(pdb_profile_dict) *0.2)
for values in all_dimer_combination:
    dimer_str  = str(values[0])+str(values[1])
    if (all_pdb_dimers_contact.count(dimer_str)) >= dimer_treshold_consideration:
        valid_dimer_combindations.append(str(values[0])+str(values[1]))


#### 5 #### GENERATE CMAPS


print(valid_dimer_combindations)
############################ TBD ##########################
#### 4 #### GLINTER RUN
#### FOR NOW PRECOMPUTED
############################ TBD ##########################



################QUESION ALPHAFOLD STRUCTURE ??
#then run glinter
####################ALREADY DONE CURRENTLY USING PRE_COMPUTED
#MAKE CMAPS
predicted_structures_dimer_cmap_dir = output_dir+"struct_dimer_cmaps/"
dimer_strcutures_dir  = output_dir+"dimer_structures_pdb/"
if not os.path.exists(predicted_structures_dimer_cmap_dir):
    os.system("mkdir -p "+predicted_structures_dimer_cmap_dir)
if not os.path.exists(dimer_strcutures_dir):
    os.system("mkdir -p "+dimer_strcutures_dir)
for pdb in pdb_profile_dict:
    temp_pdb_profile  = pdb_profile_dict.get(pdb)
    for dimer in temp_pdb_profile.dimers:
        if dimer in valid_dimer_combindations:
            monomer_a = temp_pdb_profile.chain_skeleton_CA.get(dimer[0])
            monomer_b = temp_pdb_profile.chain_skeleton_CA.get(dimer[1])
            dimer_cmap_file_name = predicted_structures_dimer_cmap_dir+str(pdb)+"_chain_"+dimer+".cmap"
            dimer_pdb_file_name = dimer_strcutures_dir + str(pdb) + "_chain_" + dimer + ".pdb"
            if not os.path.exists(dimer_cmap_file_name):
                temp_cmaps = eva_util.get_CA_cmaps(_first_chain=monomer_a,_second_chain=monomer_b)

                np.savetxt(dimer_cmap_file_name,temp_cmaps)
                temp_dimer_array = copy.deepcopy(monomer_a)+copy.deepcopy(monomer_b)
                eva_util.pdb_from_array(_pdb=temp_dimer_array,_filename=dimer_pdb_file_name)

            ####THEN CONCAT DIMER PDB AND SAVE
#then calculate the other score

for pdb in pdb_profile_dict:
    temp_pdb_profile  = pdb_profile_dict.get(pdb)
    a_dimer_score_dict = {}

    for dimer in  valid_dimer_combindations:
        dimer_chain_scores = [ ]
        if dimer in temp_pdb_profile.dimers:
            chain_first =dimer_strcutures_dir+pdb+"_chain_"+str(dimer)+".pdb"
            all_specific_dimer_chain = eva_util.specific_filename_reader(_input_dir=dimer_strcutures_dir,
                                                                          _extension="_chain_"+str(dimer))
            for predicted_dimer in all_specific_dimer_chain:
                chain_second = dimer_strcutures_dir + predicted_dimer + ".pdb"
                if chain_second != chain_first:
                    dimer_chain_scores.append(eva_util.get_dock_q_score(_true=chain_first,_current=chain_second))
            a_dimer_score_dict[dimer] = np.average(dimer_chain_scores)
        else:
            a_dimer_score_dict[dimer] = 0.0
    pdb_profile_dict.get(pdb).ds_scores = a_dimer_score_dict
print("HERE")
            # eva_util.





#### 6 #### COMPARE WITH PREDICTION THEN JOB DONE
#### 7 ####TRUE EVALUATIONS
