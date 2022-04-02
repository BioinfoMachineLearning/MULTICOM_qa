import copy
import os
import eva_utils

def specific_filename_reader(_input_dir, _extension):
    file_names = []
    for root, directories, files in os.walk(_input_dir):
        for file in files:
            if _extension in file:

                if not file in file_names:
                    file_names.append(file)
    return file_names

def pdb_filtering(_input_dir,_output_dir):

    input_dir = _input_dir

    output_dir = _output_dir
    pdb_list = specific_filename_reader(_input_dir=input_dir, _extension="")
    for values in pdb_list:
        temp_pdb = input_dir + values

        temp_fasta_values = eva_utils.contents_to_info(eva_utils.read_pdb(temp_pdb))
        # print(len(temp_fasta_values))
        CA_first_chain = list(filter(lambda x: (x.atom_name == "CA"), copy.deepcopy(temp_fasta_values)))
        # print(len(CA_first_chain))

        str_pdb = eva_utils.pdb_from_array(_filename=output_dir + values.split(".")[0], _pdb=CA_first_chain)
# input_dir = "/home/bdmlab/H1036/H1036_pred/"
# output_dir = "/home/bdmlab/H1036/filter_pred/"
# pdb_filtering(input_dir,output_dir)
