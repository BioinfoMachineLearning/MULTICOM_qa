import glob
import os
import config as c_path

import eva_utils as eva_util
monomer_sequences_dir = ""
input_dir = "/home/rajroy/Multimet_evatest_samples/predictions/T1032o/"
stoichiometry = "A2"
output_dir = "/home/rajroy/multi_eva_test/"



stoi_details = eva_util.get_stoichiometry_details(stoichiometry)
predicted_pdb_files = eva_util.specific_filename_reader(input_dir,"")

predicted_monomer_dir = output_dir+"monomer/"
os.system("mkdir -p "+predicted_monomer_dir)
all_chains_discovered= []
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
for chain_values in all_chains_discovered:
    temp_chain_dir =  predicted_monomer_chains_dir+str(chain_values)
    os.system("mkdir -p " +temp_chain_dir)
    print(predicted_monomer_dir + "/**/*"+"_chain_"+str(chain_values))
    all_monomer_chained_files = glob.glob(predicted_monomer_dir + "/**/*"+"_chain_"+str(chain_values)+".pdb", recursive = True)
    for values in all_monomer_chained_files:
        os.system("cp "+values+" "+temp_chain_dir)

    if chain_values !="A":
        cmd = c_path.PARIWISE_QA_SCRIPT +" "+predicted_monomer_chains_dir+str(chain_values)+" "

#command for pairwise_qa

# for chain_values in all_chains_discovered

# perl pairwise_model_eva.pl /home/rajroy/multi_eva_test/monomer_chains/A/ /home/rajroy/multi_eva_test/Multimet_evatest_samples/casp_fasta/H1036A.fasta /home/rajroy/pairwiseQA/q_score /home/rajroy/Downloads/tools/TMscore A /home/rajroy/q_A/



#	die "need six parameters: input model dir, fasta sequence file, pairwise_QA path, tm_score path, target name (output name), output dir.\n";
# print(all_monomer_files)

#diff chains name needed
#break the multimers into monomer and get a dir for them

#lets read all of them

