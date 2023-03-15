import csv
import os
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
    rm_list = []
    for item in contents.decode("utf-8").split("\n"):

        if "TM-score=" in item:
            tm_list.append(float(item.strip().split(" ")[1].strip()))
        if "RMSD=" in item:
            rm_list.append(float(item.strip().split("RMSD=")[1].split(",")[0].strip()))
    print(os.path.basename(_current))
    print(tm_list)
    return np.min(tm_list),rm_list[0]

def specific_filename_reader(_input_dir, _extension):
    file_names = []
    for root, directories, files in os.walk(_input_dir):
        for file in files:
            if _extension in file:
                file = file.split(".")[0]
                if not file in file_names:
                    file_names.append(file)
    return file_names

#
output_dir = "/home/bdmlab/dokcG_new/DOKC_G_RMSD/"
MM_ALIGN_PATH = "/home/bdmlab/tools/MMalign"
# list_of_file ="/home/bdmlab/MM_Eva/casp14/predictions_cleaned/casp_14_oligomers.txt"
list_of_file ="/home/bdmlab/dock_G_valid_all.txt"

# output_dir = "/home/bdmlab/BMM55_new_true/"
# MM_ALIGN_PATH = "/home/bdmlab/tools/MMalign"
# list_of_file ="//home/bdmlab/DPROQ_NATIVE/bm_55_targets.txt"

list_file = file_reader(list_of_file)
# list_file =["H1111"]
for Title in list_file:
    print(Title)
    # true_dir = "//home/bdmlab/MM_Eva/casp14/pdbs_cleaned/"
    true_dir = "/home/bdmlab/MM_Eva/DG/DockG_native/"
    temp_input_dir = "/home/bdmlab/tidy_pdb/" + str(Title) + "/"
    # temp_input_dir = "/home/bdmlab/MM_Eva/DG/tidy_pdb_2/"
    true_pdb = true_dir + "/"  + str(Title) + ".pdb"

    #
    # true_dir = "/home/bdmlab/DPROQ_NATIVE/after_cleaning/"
    # temp_input_dir = "/home/bdmlab/Downloads/DproQ_benchmark/DProQ_benchmark/BM55-AF2/decoy/"+str(Title)+"/"
    # true_pdb = true_dir +"/" +str(Title)+".pdb"
    temp_pdb_list = specific_filename_reader(_input_dir=temp_input_dir, _extension="")
    score_array = []
    for values in temp_pdb_list:
        temp_pdb = temp_input_dir + str(values)+".pdb"
        value = get_MM_score(true_pdb, temp_pdb, MM_ALIGN_PATH)
        tm = value[0]
        rmsd =value[1]
        # score_array.append([str(values), get_MM_score(true_pdb, temp_pdb, MM_ALIGN_PATH)])
        score_array.append([str(values), tm,rmsd])

    name_of_output_file = output_dir+str(Title)+".csv"
    _header_row=['Target','MMalign_Score','RMSD']
    with open(name_of_output_file, 'w') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(_header_row)
        for data in score_array:
            filewriter.writerow(data)


