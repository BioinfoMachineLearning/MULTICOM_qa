import copy
import glob
import os
# from config import PARIWISE_QA_SCRIPT,Q_SCORE,TM_SCORE_PATH
import  config as config_path
PARIWISE_QA_SCRIPT = config_path.config.PARIWISE_QA_SCRIPT
TM_SCORE_PATH=config_path.config.TM_SCORE_PATH
Q_SCORE = config_path.config.Q_SCORE
# print(config_path.config.PARIWISE_QA_SCRIPT)

import eva_utils as eva_util
monomer_sequences_dir = "/home/bdmlab/T1032o.fasta"
input_dir = "/home/bdmlab/Multimet_evatest_samples/predictions/T1032o_lite/"
stoichiometry = "A2"
output_dir = "/home/bdmlab/multi_eva_test/"
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

# for chain_values in all_chains_discovered
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
##fasta_to_chain_mapper(_fasta_file="/home/bdmlab/T1032o.fasta", _stoi="A2", _chains=["A", "B"]):

#	die "need six parameters: input model dir, fasta sequence file, pairwise_QA path, tm_score path, target name (output name), output dir.\n";
# print(all_monomer_files)

#diff chains name needed
#break the multimers into monomer and get a dir for them

#lets read all of them

