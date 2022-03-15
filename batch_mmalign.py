import csv
import subprocess

import numpy as np

import eva_utils


def file_reader(input_dir):
    contents = ""
    f = open(input_dir, "r")
    if f.mode == 'r':
        contents = f.read().splitlines()
        f.close()
    return contents

def get_MM_score(_true, _current, _MM_ALIGN):
    MM_ALIGN_PATH = _MM_ALIGN
    contents = subprocess.check_output([MM_ALIGN_PATH, _true, _current])
    tm_list = []

    for item in contents.decode("utf-8").split("\n"):

        if "TM-score=" in item:
            tm_list.append(float(item.strip().split(" ")[1].strip()))

    return np.min(tm_list)



output_dir = "/home/bdmlab/mm_align_output/"
MM_ALIGN_PATH = "/home/bdmlab/tools/MMalign"
list_of_file ="/home/bdmlab/predictions_cleaned/casp_14_oligomers.txt"

list_file = file_reader(list_of_file)

for Title in list_file:
    print(Title)
    true_dir = "/home/bdmlab/pdbs_cleaned/"
    temp_input_dir = "/home/bdmlab/predictions_cleaned/"+str(Title)+"/"
    true_pdb = true_dir + str(Title) + "/" + str(Title) + ".pdb"
    temp_pdb_list = eva_utils.specific_filename_reader(_input_dir=temp_input_dir, _extension="")
    score_array = []
    for values in temp_pdb_list:
        temp_pdb = temp_input_dir + str(values)
        score_array.append([str(values), float(get_MM_score(true_pdb, temp_pdb, MM_ALIGN_PATH))])

    name_of_output_file = output_dir+str(Title)+".csv"
    _header_row=['Target','MMalign_Score']
    with open(name_of_output_file, 'w') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(_header_row)
        for data in score_array:
            filewriter.writerow(data)


