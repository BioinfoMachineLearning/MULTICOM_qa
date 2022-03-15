import os
import subprocess
import numpy as np

tmscore_path = "/home/bdmlab/Documents/tools/TMalign"
input_dir ="/home/bdmlab/concatenated_pdb/"
native_pdb = "/home/bdmlab/H1045.pdb"
output_dir ="/home/rajroy/tmscore/"


def specific_filename_reader(_input_dir, _extension):
    file_names = []
    for root, directories, files in os.walk(_input_dir):
        for file in files:
            if _extension in file:
                file = file.split(".")[0]
                if not file in file_names:
                    file_names.append(file.split(".")[0])
    print(file_names)
    return file_names


proteinNames = []  # initializing array
# reads all the folder in a directory
proteinNames =specific_filename_reader(input_dir,"")
value = []
DOCK_Q_PATH = "/home/bdmlab/tools/DockQ/DockQ.py"

def get_dock_q_score(_true="/home/bdmlab/multi_eva_test/T1038/dimer_structures_pdb/T1038TS029_1o_chain_AB.pdb",
                     _current="/home/bdmlab/multi_eva_test/T1038/dimer_structures_pdb/T1038TS062_3o_chain_AB.pdb"):
    contents = subprocess.check_output([DOCK_Q_PATH, _true, _current])
    dock_q_score = 0
    for item in contents.decode("utf-8").split("\n"):
        not_first_dock = True
        if "DockQ " in item:
            if len(item.strip().split(" ")) == 2:
                dock_q_score = item.strip().split(" ")[1].strip()
                return float(dock_q_score)

for proteins in proteinNames:


    pred = input_dir + proteins
    native_file =native_pdb
    if os.path.isfile(pred) and os.path.isfile(native_file):
        cmd = tmscore_path + " " + pred+ " "+native_file
        # print(cmd)
        # os.system(cmd + " &> " + output_dir+proteins+ ".txt" + "\n")
        # contents = subprocess.check_output([tmscore_path, pred, native_file])
        # tmscore = ""
        # for item in contents.decode("utf-8").split("\n"):
        #     if "TM-score=" in item:
        #         tmscore = item.strip().split("=")[1].split("(")[0]
        #         value.append(float(tmscore.strip()))
        #         print(tmscore)
        # print(proteins)
        if "XNCX" in proteins:
            print (proteins+str(",")+str(0))
        else:
            print (proteins+str(",")+str(get_dock_q_score(_true=native_file,_current=pred)))

# print(np.average(np.array(value)))
