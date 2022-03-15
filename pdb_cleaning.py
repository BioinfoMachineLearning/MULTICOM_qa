import copy
import sys

import eva_utils


def pdb_filtering(_input_dir,_output_dir):

    input_dir = _input_dir

    output_dir = _output_dir
    pdb_list = eva_utils.specific_filename_reader(_input_dir=input_dir, _extension="")
    for values in pdb_list:
        temp_pdb = input_dir + values

        temp_fasta_values = eva_utils.contents_to_info(eva_utils.read_pdb(temp_pdb))
        print(len(temp_fasta_values))
        CA_first_chain = list(filter(lambda x: (x.atom_name == "CA"), copy.deepcopy(temp_fasta_values)))
        print(len(CA_first_chain))

        str_pdb = eva_utils.pdb_from_array(_filename=output_dir + values, _pdb=CA_first_chain)
# input_dir = "/home/bdmlab/H1036/H1036_pred/"
# output_dir = "/home/bdmlab/H1036/filter_pred/"
# pdb_filtering(input_dir,output_dir)
