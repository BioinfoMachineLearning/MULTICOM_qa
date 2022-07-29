import csv
import os
import subprocess

import numpy as np

import eva_utils
DOCK_Q_PATH ="/home/bdmlab/tools/DockQ/DockQ.py"

def file_reader(input_dir):
    contents = ""
    f = open(input_dir, "r")
    if f.mode == 'r':
        contents = f.read().splitlines()
        f.close()
    return contents

class dockq_profile:
    # Monomer Score (MS)
    # Dimer Score (DS)
    # Interchain contact probability scores (ICPS):
    # Recall
    name = ""
    dockq = 0.0
    Fnat = 0.0
    Fnonnat = 0.0
    iRMS = 0.0
    LRMS=0.0

    pass


def get_dock_q_score(_true="/home/bdmlab/multi_eva_test/T1038/dimer_structures_pdb/T1038TS029_1o_chain_AB.pdb",
                     _current="/home/bdmlab/multi_eva_test/T1038/dimer_structures_pdb/T1038TS062_3o_chain_AB.pdb"):
    contents = subprocess.check_output([DOCK_Q_PATH, _true, _current])
    dock_q_score = 0
    temp = dockq_profile()
    temp.name = os.path.basename(_current)
    for item in contents.decode("utf-8").split("\n"):
        not_first_dock = True
        if "DockQ " in item:
            if len(item.strip().split(" ")) == 2:
                dock_q_score = item.strip().split(" ")[1].strip()
                # return float(dock_q_score)
                temp.dockq=float(dock_q_score)
        if "Fnat" in item and "(Fnat)" not in item:
            temp.Fnat = float(item.split(" ")[1])
        if "iRMS" in item and "(iRMS)" not in item:
            temp.iRMS = float(item.split(" ")[1])
        if "Fnonnat" in item :
            temp.Fnonnat = float(item.split(" ")[1])
        if "LRMS" in item:
            temp.LRMS = float(item.split(" ")[1])
    return [temp.Fnat,temp.Fnonnat,temp.iRMS,temp.LRMS,temp.dockq]

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
output_dir = "/home/bdmlab/HAF2_dockQ/"
MM_ALIGN_PATH = "/home/bdmlab/tools/MMalign"
# list_of_file ="/home/bdmlab/MM_Eva/casp14/predictions_cleaned/casp_14_oligomers.txt"
list_of_file ="/home/bdmlab/haf2.list"

# output_dir = "/home/bdmlab/BMM55_new_true/"
# MM_ALIGN_PATH = "/home/bdmlab/tools/MMalign"
# list_of_file ="//home/bdmlab/DPROQ_NATIVE/bm_55_targets.txt"

list_file = file_reader(list_of_file)

for Title in list_file:
    print(Title)
    # true_dir = "//home/bdmlab/MM_Eva/casp14/pdbs_cleaned/"
    true_dir = "/home/bdmlab/Downloads/DproQ_benchmark_uncompresed.tgz/DProQ_benchmark/HAF2/native/"
    temp_input_dir = "/home/bdmlab/Downloads/DproQ_benchmark_uncompresed.tgz/DProQ_benchmark/HAF2/decoy/" + str(Title) + "/pdb/"
    true_pdb = true_dir + "/"  + str(Title) + ".pdb"

    #
    # true_dir = "/home/bdmlab/DPROQ_NATIVE/after_cleaning/"
    # temp_input_dir = "/home/bdmlab/Downloads/DproQ_benchmark/DProQ_benchmark/BM55-AF2/decoy/"+str(Title)+"/"
    # true_pdb = true_dir +"/" +str(Title)+".pdb"
    temp_pdb_list = specific_filename_reader(_input_dir=temp_input_dir, _extension="pdb")
    score_array = []
    for values in temp_pdb_list:
        temp_pdb = temp_input_dir + str(values)+".pdb"
        value = get_dock_q_score(true_pdb, temp_pdb)
        Fnat=value[0]
        Fnonnat=value[1]
        iRMSm=value[2]
        LRMS=value[3]
        dockq=value[4]
        # score_array.append([str(values), get_MM_score(true_pdb, temp_pdb, MM_ALIGN_PATH)])
        score_array.append([str(values), Fnat,Fnonnat,iRMSm,LRMS,dockq])

    name_of_output_file = output_dir+str(Title)+".csv"
    _header_row=['Target','Fnat','Fnonnat','iRMS','LRMS','dockq']
    with open(name_of_output_file, 'w') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(_header_row)
        for data in score_array:
            filewriter.writerow(data)


