import numpy as np
import copy
import csv
import math
import os
import re
import subprocess
import time
from Bio import pairwise2
# from PIL import Image as im
import concurrent.futures


class multimer:
    unique_monomers_chain = []
    chains = [],
    chain_fasta = {},
    stoi = [],
    # grp_id = ""
    scores = {}

    pass


fasta_3_to_1_code = {'ala': 'A', 'arg': 'R', 'asn': 'N', 'asp': 'D', 'asx': 'B', 'cys': 'C', 'glu': 'E', 'gln': 'Q',
                     'glx': 'Z', 'gly': 'G', 'his': 'H', 'ile': 'I', 'leu': 'L', 'lys': 'K', 'met': 'M', 'phe': 'F',
                     'pro': 'P', 'ser': 'S', 'thr': 'T', 'trp': 'W', 'tyr': 'Y', 'val': 'V'}


def write2File(_filename, _cont):
    with open(_filename, "w") as f:
        f.writelines(_cont)
        f.close()


def get_stoichiometry_details(_stoi):
    details = {}
    for i in range(int(len(_stoi) / 2)):
        value = i * 2
        details[str(_stoi[value])] = _stoi[value + 1]
    return details


def specific_filename_reader(_input_dir, _extension):
    file_names = []
    for root, directories, files in os.walk(_input_dir):
        for file in files:
            if _extension in file:
                file = file.split(".")[0]
                if not file in file_names:
                    file_names.append(file.split(".")[0])
    return file_names


def space_returner(_input):
    i = 0
    space = ""
    while i < _input:
        space = space + " "
        i = i + 1
    return space


def added_warning_logs(_file, _msg):
    # Open a file with access mode 'a'
    file_object = open(_file, 'a')
    # Append 'hello' at the end of file
    file_object.write(_msg)
    # Close the file
    file_object.close()
    return


def find_lowest_gap(_target, _hit):
    aln_val = pairwise2.align.globalms(_target, _hit, 5, -4, -1, -0.1)
    chain_target = list(aln_val[0][0])
    chain_hit = list(aln_val[0][1])
    # print(chain_target)
    # print(chain_hit)

    return chain_hit.count('-') / len(chain_hit)


def convert_to_pdb(_pdb, _name):
    content = ''
    for x in _pdb:
        content += correct_format(x) + '\n'
    f = open(_name, "w")
    f.write(content)
    f.close()
    return _pdb


def closest_key(_seq_fasta_dict, _fasta_string):
    val = []
    for key in _seq_fasta_dict:
        val.append(find_lowest_gap(_seq_fasta_dict.get(key), _fasta_string))
    seq = min(val)
    index_closest = val.index(seq)
    return index_closest


def sequence_finder(_seq_fasta_dict, _fasta_string):
    for key in _seq_fasta_dict:
        temp_fasta = _seq_fasta_dict.get(key)
        if temp_fasta == _fasta_string:
            return key
    seq_ = closest_key(_seq_fasta_dict, _fasta_string)
    # print(" closest_key ")
    return seq_


class predicted_pdb_profile:
    # Monomer Score (MS)
    # Dimer Score (DS)
    # Interchain contact probability scores (ICPS):
    # Recall
    name = ""
    dimers = []
    multimer_scoring = 0.0
    monomers_chains = []
    chain_skeleton_CA = []
    chain_fasta = []
    ds_scores = {}
    ms_scores = {}
    icps_scores = {}
    icps_rank = 0
    mm_align_rank = 0
    final_rank = 0
    recall = {}
    cluster_chain = {}
    chain_cluster = {}

    pass


def dimer_for_cmaps(_valid_dimer_combos, _temp_pdb_profile):
    list_dimer = []
    # _temp_list = [ ]
    # for values in _temp_pdb_profile.dimers:
    #     print(values)
    #     _temp_list.append(values)
    #
    # for values in _valid_dimer_combos:
    #     list_1 = _temp_pdb_profile.cluster_chain.get(values[0])
    #     list_2 = _temp_pdb_profile.cluster_chain.get(values[1])
    #     for l1_value in list_1:
    #         for l2_value in list_2:
    #             if list_1!=list_2:
    #                 list_dimer.append(str(l1_value)+str(l2_value))
    for values in _temp_pdb_profile.dimers:
        chain_1 = _temp_pdb_profile.chain_cluster.get(values[0])
        chain_2 = _temp_pdb_profile.chain_cluster.get(values[1])
        literal_chain_value = str(chain_1) + str(chain_2)
        if literal_chain_value in _valid_dimer_combos:
            list_dimer.append(values)
        else:
            literal_chain_value = str(chain_2) + str(chain_1)
            if literal_chain_value in _valid_dimer_combos:
                list_dimer.append(values)

    return list_dimer


def dir_maker(_dir_name):
    if not os.path.exists(_dir_name):
        os.system("mkdir -p " + _dir_name)
        return _dir_name
    else:
        print("Already exists ")
        return _dir_name


class pdb_lines:
    atom = ''
    serial = ''
    atom_name = ''
    alt_loc = ''
    res_name = ''
    chain = ''
    res_num = ''
    icode = ''
    x = ''
    y = ''
    z = ''
    occupancy = ''
    temp_fact = ''
    element = ''
    charge = ''

    pass


def split_line_to_tuple(line):
    a_pdb_line = pdb_lines()

    a_pdb_line.atom = line[0:6].strip()
    a_pdb_line.serial = line[6:12].strip()
    a_pdb_line.atom_name = line[12:16].strip()
    a_pdb_line.alt_loc = line[16].strip()
    a_pdb_line.res_name = line[17:20].strip()
    a_pdb_line.chain = line[20:22].strip()
    a_pdb_line.res_num = line[22:26].strip()
    a_pdb_line.icode = line[26:30].strip()
    a_pdb_line.x = line[30:38].strip()
    a_pdb_line.y = line[38:46].strip()
    a_pdb_line.z = line[46:54].strip()
    a_pdb_line.occupancy = line[54:60].strip()
    # a_pdb_line.temp_fact = line[60:76].strip()
    a_pdb_line.temp_fact = line[60:66].strip()
    a_pdb_line.element = line[76:78].strip()
    a_pdb_line.charge = line[78:80].strip()

    return a_pdb_line


def string_array_from_pdb_array(_pdb_row):
    _pdb_copy = copy.deepcopy(_pdb_row)
    # https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
    _pdb_copy.atom = _pdb_copy.atom  # 1-4
    _pdb_copy.serial = space_returner(4 - len(str(_pdb_copy.serial))) + str(_pdb_copy.serial)  # 7-11
    _pdb_copy.atom_name = _pdb_copy.atom_name + space_returner(3 - len(_pdb_copy.atom_name))  # 13-16
    _pdb_copy.alt_loc = space_returner(1 - len(_pdb_copy.alt_loc)) + _pdb_copy.alt_loc  # 17
    _pdb_copy.res_name = space_returner(3 - len(_pdb_copy.res_name)) + _pdb_copy.res_name  # 18-20
    _pdb_copy.chain = space_returner(1 - len(_pdb_copy.chain)) + _pdb_copy.chain  # 22
    _pdb_copy.res_num = space_returner(4 - len(_pdb_copy.res_num)) + _pdb_copy.res_num  # 23-26
    # _pdb_copy.icode = space_returner(2 - len(_pdb_copy.chain)) + _pdb_copy.icode  # 27
    _pdb_copy.icode = space_returner(1 - len(_pdb_copy.icode)) + _pdb_copy.icode  # 27
    _pdb_copy.x = space_returner(8 - len(_pdb_copy.x)) + _pdb_copy.x  # 31-38
    _pdb_copy.y = space_returner(8 - len(_pdb_copy.y)) + _pdb_copy.y  # 39-46
    _pdb_copy.z = space_returner(8 - len(_pdb_copy.z)) + _pdb_copy.z  # 47-54
    _pdb_copy.occupancy = space_returner(6 - len(_pdb_copy.occupancy)) + _pdb_copy.occupancy  # 55-60
    _pdb_copy.temp_fact = space_returner(6 - len(_pdb_copy.temp_fact)) + _pdb_copy.temp_fact  # 61-66
    _pdb_copy.element = space_returner(4 - len(_pdb_copy.element)) + _pdb_copy.element  # 73-76
    _pdb_copy.charge = space_returner(2 - len(_pdb_copy.charge)) + _pdb_copy.charge  # 77-78

    content = _pdb_copy.atom + space_returner(7 - len(_pdb_copy.serial)) + _pdb_copy.serial
    if len(_pdb_copy.atom_name) < 4:
        content = content + space_returner(2) + _pdb_copy.atom_name
    elif len(_pdb_copy.atom_name) == 4:
        content = content + " " + _pdb_copy.atom_name

    content = content + _pdb_copy.alt_loc + _pdb_copy.res_name + space_returner(
        1) + _pdb_copy.chain + _pdb_copy.res_num + _pdb_copy.icode + space_returner(
        3) + _pdb_copy.x + _pdb_copy.y + _pdb_copy.z + _pdb_copy.occupancy + _pdb_copy.temp_fact + space_returner(
        8) + _pdb_copy.element + _pdb_copy.charge
    return content

    # content = _pdb_copy.atom + space_returner(3) + _pdb_copy.serial + space_returner(
    #     2) + _pdb_copy.atom_name + _pdb_copy.alt_loc + _pdb_copy.res_name + space_returner(
    #     1) + _pdb_copy.chain + _pdb_copy.res_num + _pdb_copy.icode + space_returner(
    #     3) + _pdb_copy.x + _pdb_copy.y + _pdb_copy.z + _pdb_copy.occupancy + _pdb_copy.temp_fact + space_returner(
    #     6) + _pdb_copy.element
    # return content


def read_fasta(_fasta):
    file = open(_fasta, "r")
    output_array = []
    if file.mode == 'r':
        output_array = file.read().strip().splitlines()
        file.close()

    return output_array[1]


def read_pdb(pdb):
    contents = []
    with open(pdb, "r") as f:
        for line in f:
            if (line.startswith("ATOM")):
                # pass
                contents.append(line)
    return contents


def chain_replacer(_pdb_file, _new_chain_name):
    _temp = copy.deepcopy(_pdb_file)
    for value in _temp:
        value.chain = _new_chain_name
    return _temp


def pdb_from_array(_pdb, _filename):
    array = []
    content = ''
    number = 1
    for x in _pdb:
        x.serial = number
        val = string_array_from_pdb_array(x)
        array.append(val)
        content = content + val + '\n'
        number = number + 1
    f = open(_filename, "w")
    f.write(content + 'END')
    f.close()
    return array


def contents_to_info(contents):  # reads the ATOM line. Then splits the info into respective frames and returns the data
    split_contents = []
    for lines in contents:
        if lines.startswith("ATOM"):
            pdb_line = split_line_to_tuple(lines.strip())
            split_contents.append(pdb_line)
    return split_contents


def separate_by_chain(_pdb, _name):
    # print(_pdb)
    result = list(filter(lambda x: (x.chain == _name), _pdb))
    return result


def fix_serial(_array, _no=1):
    number = _no
    for x in _array:
        x.serial = number
        number = number + 1
    return _array


# MM_ALIGN_PATH = "/home/bdmlab/Documents/tools/MMalign"
# def get_MM_score(_true="/home/bdmlab/multi_eva_test/T1038/dimer_structures_pdb/T1038TS029_1o_chain_AB.pdb",
#                  _current="/home/bdmlab/multi_eva_test/T1038/dimer_structures_pdb/T1038TS062_3o_chain_AB.pdb",
#                  _MM_ALIGN="/home/bdmlab/Documents/tools/MMalign"):
#     MM_ALIGN_PATH = _MM_ALIGN
#     contents = subprocess.check_output([MM_ALIGN_PATH, _true, _current])
#     tm_list = []
#
#     for item in contents.decode("utf-8").split("\n"):
#
#         if "TM-score=" in item:
#             tm_list.append(float(item.strip().split(" ")[1].strip()))
#
#     return np.min(tm_list)


def get_MM_score(_arr):
    MM_ALIGN_PATH = _arr[2]
    # print(MM_ALIGN_PATH, _arr[0], _arr[1])
    contents = subprocess.check_output([MM_ALIGN_PATH, _arr[0], _arr[1]])
    try:
        tm_list = []

        for item in contents.decode("utf-8").split("\n"):

            if "TM-score=" in item:
                tm_list.append(float(item.strip().split(" ")[1].strip()))

        return np.min(tm_list)
    except:
        return 0.0


def get_MM_score_parallel_submit(_array, _CPU_COUNT):
    all_value = []
    worker = int(_CPU_COUNT)
    with concurrent.futures.ThreadPoolExecutor(max_workers=worker) as executor:
        result_futures = list(map(lambda x: executor.submit(get_MM_score, x), _array))
        for future in concurrent.futures.as_completed(result_futures):
            try:
                # print('resutl is', future.result())
                all_value.append(future.result())
                # print(type(future.result()))
            except Exception as e:
                print('e is', e, type(e))
                all_value.append(0)

    return np.average(all_value)


def correct_format(_pdb_row):
    _pdb_copy = copy.deepcopy(_pdb_row)
    # https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
    _pdb_copy.atom = _pdb_copy.atom  # 1-4
    _pdb_copy.serial = space_returner(5 - len(str(_pdb_copy.serial))) + str(_pdb_copy.serial)  # 7-11
    _pdb_copy.atom_name = _pdb_copy.atom_name + space_returner(3 - len(_pdb_copy.atom_name))  # 13-16
    _pdb_copy.alt_loc = space_returner(1 - len(_pdb_copy.alt_loc)) + _pdb_copy.alt_loc  # 17
    _pdb_copy.res_name = space_returner(3 - len(_pdb_copy.res_name)) + _pdb_copy.res_name  # 18-20
    _pdb_copy.chain = space_returner(1 - len(_pdb_copy.chain)) + _pdb_copy.chain  # 22
    _pdb_copy.res_num = space_returner(4 - len(_pdb_copy.res_num)) + _pdb_copy.res_num  # 23-26
    # _pdb_copy.icode = space_returner(2 - len(_pdb_copy.chain)) + _pdb_copy.icode  # 27
    _pdb_copy.icode = space_returner(1 - len(_pdb_copy.icode)) + _pdb_copy.icode  # 27
    _pdb_copy.x = space_returner(8 - len(_pdb_copy.x)) + _pdb_copy.x  # 31-38
    _pdb_copy.y = space_returner(8 - len(_pdb_copy.y)) + _pdb_copy.y  # 39-46
    _pdb_copy.z = space_returner(8 - len(_pdb_copy.z)) + _pdb_copy.z  # 47-54
    _pdb_copy.occupancy = space_returner(6 - len(_pdb_copy.occupancy)) + _pdb_copy.occupancy  # 55-60
    _pdb_copy.temp_fact = space_returner(6 - len(_pdb_copy.temp_fact)) + _pdb_copy.temp_fact  # 61-66
    _pdb_copy.element = space_returner(4 - len(_pdb_copy.element)) + _pdb_copy.element  # 73-76
    _pdb_copy.charge = space_returner(2 - len(_pdb_copy.charge)) + _pdb_copy.charge  # 77-78
    content = _pdb_copy.atom + space_returner(2) + _pdb_copy.serial

    if len(_pdb_copy.atom_name) < 4:
        content = content + space_returner(2) + _pdb_copy.atom_name
    elif len(_pdb_copy.atom_name) == 4:
        content = content + " " + _pdb_copy.atom_name

    content = content + _pdb_copy.alt_loc + _pdb_copy.res_name + space_returner(
        1) + _pdb_copy.chain + _pdb_copy.res_num + _pdb_copy.icode + space_returner(
        3) + _pdb_copy.x + _pdb_copy.y + _pdb_copy.z + _pdb_copy.occupancy + _pdb_copy.temp_fact + space_returner(
        8) + _pdb_copy.element + _pdb_copy.charge

    return content
def get_unique_chains(_inp_details):
    chain_array = []
    for val in _inp_details:
        chain_array.append(val.chain)
    return list(dict.fromkeys(chain_array))


def save_multi_fasta(_dir, _map):
    for key in _map:
        fasta_value = _map.get(key)
        # print(key, fasta_value)
        fasta_value = ">sequence_" + key + "\n" + fasta_value
        write2File(_dir + "sequence_" + key + "_A.fasta", fasta_value)
        write2File(_dir + "sequence_" + key + "_B.fasta", fasta_value)


def multi_fasta_reader(_seq_file):
    file = open(_seq_file, "r")
    stoi_fasta_dict = {}
    output_array = []
    if file.mode == 'r':
        output_array = file.read().strip().splitlines()
        file.close()
    counter = 0
    temp_name = ""
    for values in output_array:
        # print(str(counter) + " " + values)
        if counter % 2 == 1:
            stoi_fasta_dict[temp_name] = values
        else:
            temp_name = str(int(counter / 2))
        counter = counter + 1

    return stoi_fasta_dict


def get_fasta_from_pdb_array(_pdb):
    pdb_a = _pdb
    index_tracker_a = []
    for val in pdb_a:
        index_tracker_a.append(str(val.res_num) + "_" + str(val.res_name))
    index_tracker_a = list(dict.fromkeys(index_tracker_a))
    # print(index_tracker_a)
    fasta_string = ''
    for values in index_tracker_a:
        three_code = values.split('_')[1].lower()
        fasta_string = fasta_string + str(fasta_3_to_1_code.get(three_code))
    # print(len(fasta_string))
    return fasta_string


def monomer_pdb_filtering(_pdb, _dir):
    tar_name = os.path.basename(_pdb)
    tar_dir = _dir + "/" + tar_name
    os.system("mkdir -p " + tar_dir)
    full_pdb = contents_to_info(read_pdb(_pdb))
    chain_finder = get_unique_chains(full_pdb)
    # print(chain_finder)
    for chain in chain_finder:
        temp_monomer_pdb = separate_by_chain(full_pdb, chain)
        tar_monomer_file = tar_dir + "/" + tar_name + "_chain_" + str(chain) + ".pdb"
        monomer_string = pdb_from_array(_pdb=temp_monomer_pdb, _filename=tar_monomer_file)

        fasta_name = get_fasta_from_pdb_array(temp_monomer_pdb)
        fasta_value = ">sequence_" + chain + "\n" + fasta_name
        fasta_file_name = tar_dir + "/" + tar_name + "_chain_" + str(chain) + ".fasta"
        write2File(_filename=fasta_file_name, _cont=fasta_value)

    return chain_finder


def chain_pdb_combination_generator(_stoi, _chains, _fasta_stoic_dict):
    stoi_repeat_dict = {}
    stoi_chain_map_dict = {}
    counter = 0
    temp_value = ""
    string1 = _stoi
    pattern = r'[0-9]'
    only_chain_string = re.sub(pattern, '', string1)
    a_multimer = multimer()
    a_multimer.chains = copy.deepcopy(_chains)
    a_multimer.stoi = _stoi
    temp_chain_dict = {}
    string2 = _stoi
    pattern = r'[A-Za-z]'
    only_subunits = re.sub(pattern, ',', string2)
    total_subunits = 0
    counter = 0
    for values in only_subunits.split(","):
        if values.strip() != "":
            stoi_repeat_dict[only_chain_string[counter]] = int(values)
            total_subunits = total_subunits + int(values)
            counter = counter + 1
    unique_counter = 0
    while len(_chains) > 0:
        for val in stoi_repeat_dict:
            # print(val)
            if not val in stoi_chain_map_dict:
                # stoi_chain_map_dict[val] = [_fasta_stoic_dict[str(unique_counter)],chains[0]]
                stoi_chain_map_dict[val] = [_fasta_stoic_dict[str(unique_counter)], _chains[0]]
                stoi_repeat_dict[val] = int(stoi_repeat_dict[val]) - 1
                temp_chain = str(_chains[0])
                temp_fasta = str(_fasta_stoic_dict[str(unique_counter)])
                temp_chain_dict[temp_chain] = temp_fasta
                a_multimer.unique_monomers_chain.append(_chains[0])
                del _chains[0]
                _fasta_stoic_dict.pop(str(unique_counter))
                unique_counter = unique_counter + 1
            else:
                if stoi_repeat_dict[val] != 0:
                    stoi_chain_map_dict[val].append(_chains[0])
                    prev_fasta = temp_chain_dict[val]
                    temp_chain_dict[_chains[0]] = prev_fasta
                    # print(_chains[0])
                    stoi_repeat_dict[val] = int(stoi_repeat_dict[val]) - 1
                    del _chains[0]
            ##workaround for shortcut
            if len(temp_chain_dict) == total_subunits:
                _chains = []
    a_multimer.chain_fasta = temp_chain_dict
    return a_multimer


# def fasta_to_chain_mapper(_fasta_file="/home/bdmlab/H1060.fasta",
#                           chains=["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q",
#                                   "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "a"], _stoi="A6B3C12D6"):
def fasta_to_chain_mapper(_fasta_file="/home/bdmlab/T1032o.fasta", _stoi="A2", _chains=["A", "B"]):
    fasta_stoic_dict = multi_fasta_reader(_seq_file=_fasta_file)
    chain_mapper = chain_pdb_combination_generator(_stoi, _chains, fasta_stoic_dict)

    return chain_mapper


def read_skeleton(_pdb_path="/home/bdmlab/Multimet_evatest_samples/true_monomer/H1036/H1036_A1.pdb"):
    full_pdb = contents_to_info(read_pdb(_pdb_path))
    # return list(filter(lambda x: (x.atom_name == "CA"), full_pdb))
    return full_pdb


def distance(coord1, coord2):
    x1 = float(coord1["x"])
    y1 = float(coord1["y"])
    z1 = float(coord1["z"])
    x2 = float(coord2["x"])
    y2 = float(coord2["y"])
    z2 = float(coord2["z"])
    d = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)

    return d


def if_contact(_first_chain_path, _second_chain_path):
    first_chain_CA = list(filter(lambda x: (x.atom_name == "CA"), _first_chain_path))
    second_chain_CA = list(filter(lambda x: (x.atom_name == "CA"), _second_chain_path))
    for a_cord in first_chain_CA:
        for b_cord in second_chain_CA:
            # can make it a little more optimized
            if np.sqrt((float(a_cord.x) - float(b_cord.x)) ** 2 + (float(a_cord.y) - float(b_cord.y)) ** 2 + (
                    float(a_cord.z) - float(b_cord.z)) ** 2) <= 8:
                return True
    return False


def get_CA_cmaps(_first_chain, _second_chain):
    CA_first_chain = list(filter(lambda x: (x.atom_name == "CA"), copy.deepcopy(_first_chain)))
    CA_second_chain = list(filter(lambda x: (x.atom_name == "CA"), copy.deepcopy(copy.deepcopy(_second_chain))))
    chain_len_a = len(CA_first_chain)
    chain_len_b = len(CA_second_chain)
    cmap_array = np.zeros((chain_len_a, chain_len_b))
    for a_cord in range(chain_len_a):
        for b_cord in range(chain_len_b):
            # can make it a little more optimized
            dist = np.sqrt((float(CA_first_chain[a_cord].x) - float(CA_second_chain[b_cord].x)) ** 2 + (
                    float(CA_first_chain[a_cord].y) - float(CA_second_chain[b_cord].y)) ** 2 + (
                                   float(CA_first_chain[a_cord].z) - float(CA_second_chain[b_cord].z)) ** 2)
            if dist <= 8:
                cmap_array[a_cord][b_cord] = 1

    return cmap_array


# DOCK_Q_PATH = "/home/bdmlab/Documents/DockQ/DockQ.py"


# def get_dock_q_score(_true="/home/bdmlab/multi_eva_test/T1038/dimer_structures_pdb/T1038TS029_1o_chain_AB.pdb",
#                      _current="/home/bdmlab/multi_eva_test/T1038/dimer_structures_pdb/T1038TS062_3o_chain_AB.pdb",
#                      _DOCK_Q_PATH="/home/bdmlab/Documents/DockQ/DockQ.py"):
#     DOCK_Q_PATH = _DOCK_Q_PATH
#     contents = subprocess.check_output([DOCK_Q_PATH, _true, _current])
#     dock_q_score = 0
#     try:
#         for item in contents.decode("utf-8").split("\n"):
#             not_first_dock = True
#             if "DockQ " in item:
#                 if len(item.strip().split(" ")) == 2:
#                     dock_q_score = item.strip().split(" ")[1].strip()
#                     return float(dock_q_score)
#     except:
#         return 0.0


def get_dock_q_score(_inp):
    DOCK_Q_PATH = _inp[2]
    dock_q_score = 0
    try:
        contents = subprocess.check_output([DOCK_Q_PATH, _inp[0], _inp[1]])
        for item in contents.decode("utf-8").split("\n"):
            not_first_dock = True
            if "DockQ " in item:
                if len(item.strip().split(" ")) == 2:
                    dock_q_score = item.strip().split(" ")[1].strip()
                    return float(dock_q_score)
    except:
        try:
            contents = subprocess.check_output([DOCK_Q_PATH, _inp[1], _inp[0]])
            for item in contents.decode("utf-8").split("\n"):
                not_first_dock = True
                if "DockQ " in item:
                    if len(item.strip().split(" ")) == 2:
                        dock_q_score = item.strip().split(" ")[1].strip()
                        return float(dock_q_score)
        except:
            return 0


def get_dock_q_score_parallel_submit(_array, _CPU_COUNT):
    all_value = []
    worker = int(_CPU_COUNT)
    with concurrent.futures.ThreadPoolExecutor(max_workers=worker) as executor:
        result_futures = list(map(lambda x: executor.submit(get_dock_q_score, x), _array))
        for future in concurrent.futures.as_completed(result_futures):
            try:
                # print('resutl is', future.result())
                all_value.append(future.result())
                print(future.result())
                # print(type(future.result()))
            except Exception as e:
                print('e is', e, type(e))
                all_value.append(0)

    return np.average(all_value)


def get_icps_score(_struct_cmap, _pred_cmap, _transpose):
    first_cmap_copy = np.loadtxt(_struct_cmap)
    second_cmap_copy = np.loadtxt(_pred_cmap)
    # second_cmap_copy = np.load(_pred_cmap, allow_pickle=True)
    # if _transpose:
    #     first_cmap_copy = np.transpose(first_cmap_copy)

    len_a, len_b = second_cmap_copy.shape
    icps_list = []
    con_number = np.count_nonzero(first_cmap_copy)
    for i in range(con_number):
        (x, y) = np.unravel_index(np.argmax(first_cmap_copy, axis=None), first_cmap_copy.shape)
        first_cmap_copy[x][y] = 0
        if x < len_a and y < len_b:
            icps_list.append(second_cmap_copy[x][y])

    return np.average(icps_list)


# def show_cmap_image(data, _name):
#     data = im.fromarray(data)
#
#     data = data.convert("L")
#
#     data.save(_name)
#     return


def get_recall(_struct_cmap, _pred_cmap, _transpose):
    struct_cmap = np.loadtxt(_struct_cmap)
    pred_cmap = np.loadtxt(_pred_cmap)
    # pred_cmap = np.load(_pred_cmap, allow_pickle=True)
    len_a, len_b = pred_cmap.shape
    if _transpose:
        struct_cmap = np.transpose(struct_cmap)
    con_number = int(min(len_b, len_a) / 5)
    true_positive = 0
    s_len_a, s_len_b = struct_cmap.shape
    i = 0
    while i < con_number:
        (x, y) = np.unravel_index(np.argmax(pred_cmap, axis=None), pred_cmap.shape)
        pred_cmap[x][y] = 0
        if x < s_len_a and y < s_len_b:
            i = i + 1
            if int(struct_cmap[x][y]) == 1:
                true_positive = true_positive + 1
                # print(true_positive)
    return true_positive / con_number


def read_monomer_score(_path="/home/bdmlab/multi_eva_test/T1038_LITE/score/monomer/A.tm"):
    chain_name = os.path.basename(_path).replace(".tm", "")
    out_dict = {}
    file = open(_path, "r")
    output_array = []
    if file.mode == 'r':
        output_array = file.read().strip().splitlines()
        file.close()
    # for values in output_array[4:len(output_array) - 1]:
    for values in output_array:
        name = values.split(" ")[0].replace("_chain_" + str(chain_name) + ".pdb", "")
        out_dict[name] = float(values.split(" ")[1].strip())
    return out_dict


def read_mm_score(_path="/home/bdmlab/multi_eva_test/T1038_LITE/score/monomer/A.tm"):
    out_dict = {}
    file = open(_path, "r")
    output_array = []
    if file.mode == 'r':
        output_array = file.read().strip().splitlines()
        file.close()
    for values in output_array:
        name = values.split(" ")[0].strip()
        out_dict[name] = float(values.split(" ")[1].strip())
    return out_dict


def report_individual_target(_header_row, _data_array, _file_name):
    data_array = copy.deepcopy(_data_array)

    name_of_output_file = _file_name

    with open(name_of_output_file, 'w') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        # filewriter.writerow(    ['Name', 'TOP_S_L10', 'TOP_S_L5'])
        filewriter.writerow(_header_row)
        for data in data_array:
            filewriter.writerow(data)
    # print(output_dir + name_of_output_file)
    # print(name_of_output_file)


def get_preci_val(_x):
    return "{:.5f}".format(_x)


def replace_nan(_number):
    if math.isnan(_number):
        return 0
    else:
        return _number


def print_final_data(_file_name, _file_data, _chain_data):
    _data_array = []
    all_chains_discovered = copy.deepcopy(_chain_data)
    file_data = copy.deepcopy(_file_data)
    row_string = ""
    data_row = []
    for values in file_data:
        temp = file_data.get(values)

        temp_ms_score = []
        temp_icps = []
        temp_recall = []
        temp_ds_score = []
        for monomers in all_chains_discovered:
            temp_ms_score.append(replace_nan(float(temp.ms_scores.get(monomers))))
        ms = replace_nan(np.average(temp_ms_score))

        for dimers in temp.dimers:
            temp_ds_score.append(replace_nan(float(temp.ds_scores.get(dimers))))
        ds = replace_nan(np.average(temp_ds_score))

        for dimers in temp.dimers:
            temp_icps.append(replace_nan(float(temp.icps_scores.get(dimers))))
        is_c = replace_nan(np.average(temp_icps))

        for dimers in temp.dimers:
            temp_recall.append(replace_nan(float(temp.recall.get(dimers))))

        rec = replace_nan(np.average(temp_recall))

        final_score = (ms + ds + is_c + rec) / 4
        row_string = row_string + str(values) + "," + str(ms) + "," + str(ds) + "," + str(is_c) + "," + str(
            rec) + "," + str(final_score) + "\n"
        data_row.append([values, get_preci_val(ms), get_preci_val(ds), get_preci_val(is_c), get_preci_val(rec),
                         get_preci_val(final_score)])
    head_row = ['Name', 'Monomer_score', 'Dimer_score', 'ICP_score', 'recall_score', 'final_score']
    report_individual_target(_header_row=head_row, _file_name=_file_name, _data_array=data_row)


def get_header_string_lite(_stoic, _dimer):
    head_string = ["Name"]

    for _dimers in _dimer:
        head_string.append("ICPS_" + str(_dimers))

    head_string.append("average_ICPS")

    head_string.append("average_MMS")
    head_string.append("Final_score")
    head_string.append("ICPS_Rank")
    head_string.append("MM_Rank")
    head_string.append("Final_Rank")
    return head_string


def get_header_string(_stoic, _dimer):
    head_string = ["Name"]
    for monomer in _stoic:
        head_string.append("MS_" + str(monomer))

    for _dimers in _dimer:
        head_string.append("DS_" + str(_dimers))

    for _dimers in _dimer:
        head_string.append("ICPS_" + str(_dimers))

    for _dimers in _dimer:
        head_string.append("R_" + str(_dimers))
    head_string.append("average_MS")
    head_string.append("average_DS")
    head_string.append("average_ICPS")
    head_string.append("average_R")
    head_string.append("average_MMS")
    head_string.append("Final_score")
    return head_string


def get_casp_file(_data, _casp_file):
    with open(os.path.join(_casp_file.replace(".csv", "_casp_format.txt")), 'w') as f:
        f.writelines('PFRMAT QA' + '\n')
        f.writelines(f'TARGET {os.path.basename(_casp_file).split(".")[0]}' + '\n')
        f.writelines('AUTHOR MULTICOM_qa' + '\n')
        f.writelines('METHOD Multimer_eva' + '\n')
        f.writelines('MODEL 1'+ '\n')
        f.writelines('QMODE 1'+ '\n')
        for record in _data:
            # f.writelines(str(record[0]) + " " + str(record[-4])[0:7] + " " + '{:.6f}'.format(record[-6]) + "\n")
            f.writelines(str(record[0]) + " " + '{:.4f}'.format(record[-4])  + " " + '{:.6f}'.format(record[-6]) + "\n")
        f.writelines('END')

def print_final_data_mmalign(_file_name, _file_data, _chain_data):
    _data_array = []

    file_data = copy.deepcopy(_file_data)
    data_row = []
    for values in file_data:
        temp = file_data.get(values)
        data_row.append([values,temp.multimer_scoring] )
    head_row = ["Name","MMalign score"]
    # head_row = ['Name', 'Monomer_score', 'Dimer_score', 'ICP_score', 'recall_score', 'final_score']
    report_individual_target(_header_row=head_row, _file_name=_file_name, _data_array=data_row)


def print_final_data_icps(_file_name, _file_data, _chain_data, _dimer_data):
    _data_array = []
    dimer_interaction_discover = copy.deepcopy(_dimer_data)
    file_data = copy.deepcopy(_file_data)
    data_row = []

    # icps_scores = {}
    # icps_rank = {}
    # mm_align_rank = {}
    # final_rank = {}
    # recall = {}
    for values in file_data:
        temp = file_data.get(values)
        temp_icps = []
        total_values = []

        for dimers in dimer_interaction_discover:
            if temp.icps_scores.get(dimers) != None:
                temp_icps.append(replace_nan(temp.icps_scores.get(dimers)))
                total_values.append(replace_nan(temp.icps_scores.get(dimers)))
            else:
                total_values.append(0)
                temp_icps.append(0)

        total_values.append(np.average(temp_icps))
        total_values.append(temp.multimer_scoring)
        final_score =  np.average(temp_icps)
        total_values.append(final_score)
        data_row.append([values] + total_values)
    # icps -- >sorted(data_row, key=lambda x: x[-3], reverse=True)
    icps_counter = 1
    temp_copy = copy.deepcopy(data_row)
    for ic_values in sorted(temp_copy, key=lambda x: x[-3], reverse=True):
        file_data.get(ic_values[0]).icps_rank = icps_counter
        icps_counter = icps_counter + 1
    temp_copy = copy.deepcopy(data_row)
    mm_counter = 1
    for mm_values in sorted(temp_copy, key=lambda x: x[-2], reverse=True):
        file_data.get(mm_values[0]).mm_align_rank = mm_counter
        mm_counter = mm_counter + 1
    temp_copy = copy.deepcopy(data_row)
    for _values in temp_copy:
        temp = file_data.get(_values[0])
        temp.final_rank = float(0.6 * int(temp.mm_align_rank) + int(temp.icps_rank) * 0.4)
    rank_list = []
    for final_values in data_row:
        temp = file_data.get(final_values[0])
        final_values.append(temp.icps_rank)
        final_values.append(temp.mm_align_rank)
        final_values.append(temp.final_rank)

        # mmalign -->sorted(data_row, key=lambda x: x[-2], reverse=True)
    # data_row = sorted(data_row, key=lambda x: x[-3], reverse=True)
    # true_top_1_name = sorted(file_data.items(), key=lambda x: x[1].icps_scores, reverse=True)[0]
    # true_top_1_score = true_top_1_name[1].icps_scores

    head_row = get_header_string_lite(_chain_data, _dimer_data)
    # head_row = ['Name', 'Monomer_score', 'Dimer_score', 'ICP_score', 'recall_score', 'final_score']
    report_individual_target(_header_row=head_row, _file_name=_file_name, _data_array=data_row)
    top_5_arr = sorted(data_row, key=lambda x: x[-1])
    get_casp_file(_data=top_5_arr, _casp_file=_file_name)


def print_final_data_new_lite(_file_name, _file_data, _chain_data, _dimer_data):
    _data_array = []
    dimer_interaction_discover = copy.deepcopy(_dimer_data)
    file_data = copy.deepcopy(_file_data)
    data_row = []

    # icps_scores = {}
    # icps_rank = {}
    # mm_align_rank = {}
    # final_rank = {}
    # recall = {}
    for values in file_data:
        temp = file_data.get(values)
        temp_icps = []
        total_values = []

        for dimers in dimer_interaction_discover:
            if temp.icps_scores.get(dimers) != None:
                temp_icps.append(replace_nan(temp.icps_scores.get(dimers)))
                total_values.append(replace_nan(temp.icps_scores.get(dimers)))
            else:
                total_values.append(0)
                temp_icps.append(0)

        total_values.append(np.average(temp_icps))
        total_values.append(temp.multimer_scoring)
        final_score = 0.4 * np.average(temp_icps) + 0.6 * temp.multimer_scoring
        total_values.append(final_score)
        data_row.append([values] + total_values)
    # icps -- >sorted(data_row, key=lambda x: x[-3], reverse=True)
    icps_counter = 1
    temp_copy = copy.deepcopy(data_row)
    for ic_values in sorted(temp_copy, key=lambda x: x[-3], reverse=True):
        file_data.get(ic_values[0]).icps_rank = icps_counter
        icps_counter = icps_counter + 1
    temp_copy = copy.deepcopy(data_row)
    mm_counter = 1
    for mm_values in sorted(temp_copy, key=lambda x: x[-2], reverse=True):
        file_data.get(mm_values[0]).mm_align_rank = mm_counter
        mm_counter = mm_counter + 1
    temp_copy = copy.deepcopy(data_row)
    for _values in temp_copy:
        temp = file_data.get(_values[0])
        temp.final_rank = float(0.6 * int(temp.mm_align_rank) + int(temp.icps_rank) * 0.4)
    rank_list = []
    for final_values in data_row:
        temp = file_data.get(final_values[0])
        final_values.append(temp.icps_rank)
        final_values.append(temp.mm_align_rank)
        final_values.append(temp.final_rank)

        # mmalign -->sorted(data_row, key=lambda x: x[-2], reverse=True)
    # data_row = sorted(data_row, key=lambda x: x[-3], reverse=True)
    # true_top_1_name = sorted(file_data.items(), key=lambda x: x[1].icps_scores, reverse=True)[0]
    # true_top_1_score = true_top_1_name[1].icps_scores

    head_row = get_header_string_lite(_chain_data, _dimer_data)
    # head_row = ['Name', 'Monomer_score', 'Dimer_score', 'ICP_score', 'recall_score', 'final_score']
    report_individual_target(_header_row=head_row, _file_name=_file_name, _data_array=data_row)
    top_5_arr = sorted(data_row, key=lambda x: x[-1])
    get_casp_file(_data=top_5_arr, _casp_file=_file_name)


def print_final_data_new(_file_name, _file_data, _chain_data, _dimer_data):
    _data_array = []
    all_chains_discovered = copy.deepcopy(_chain_data)
    dimer_interaction_discover = copy.deepcopy(_dimer_data)
    file_data = copy.deepcopy(_file_data)
    row_string = ""
    data_row = []
    for values in file_data:
        temp = file_data.get(values)
        temp_ms_score = []
        temp_icps = []
        temp_recall = []
        temp_ds_score = []
        total_values = []
        # total_values.append(values)
        for monomers in all_chains_discovered:
            if temp.ms_scores.get(monomers) != None:
                temp_ms_score.append(replace_nan(float(temp.ms_scores.get(monomers))))
                total_values.append(replace_nan(float(temp.ms_scores.get(monomers))))
            else:
                temp_ms_score.append(0)
                total_values.append(0)

        # ms = replace_nan(temp_ms_score)

        for dimers in dimer_interaction_discover:
            if temp.ds_scores.get(dimers) != None:
                temp_ds_score.append(replace_nan(temp.ds_scores.get(dimers)))
                total_values.append(replace_nan(float(temp.ds_scores.get(dimers))))
            else:
                total_values.append(0)
                temp_ds_score.append(0)

        # ds = replace_nan(temp_ds_score)

        for dimers in dimer_interaction_discover:
            if temp.icps_scores.get(dimers) != None:
                temp_icps.append(replace_nan(temp.icps_scores.get(dimers)))
                total_values.append(replace_nan(temp.icps_scores.get(dimers)))
            else:
                total_values.append(0)
                temp_icps.append(0)

        for dimers in dimer_interaction_discover:
            if temp.recall.get(dimers) != None:
                temp_recall.append(replace_nan(temp.recall.get(dimers)))
                total_values.append(replace_nan(temp.recall.get(dimers)))
            else:
                total_values.append(0)
                temp_recall.append(0)
            # temp_recall.append(temp.recall.get(dimers))
        total_values.append(np.average(temp_ms_score))
        total_values.append(np.average(temp_ds_score))
        total_values.append(np.average(temp_icps))
        total_values.append(np.average(temp_recall))
        total_values.append(temp.multimer_scoring)
        final_score = np.average(
            [np.average(temp_ms_score), np.average(temp_ds_score), np.average(temp_icps), np.average(temp_recall),
             temp.multimer_scoring])
        total_values.append(final_score)

        data_row.append([values] + total_values)

    head_row = get_header_string(_chain_data, _dimer_data)
    # head_row = ['Name', 'Monomer_score', 'Dimer_score', 'ICP_score', 'recall_score', 'final_score']
    report_individual_target(_header_row=head_row, _file_name=_file_name, _data_array=data_row)


# print(get_recall(_struct_cmap="/home/bdmlab/hetero_test/multi/struct_dimer_cmaps/H1045TS285_3_chain_AB.cmap", _pred_cmap="/home/bdmlab/test/.cmap"))
# print(get_recall(_struct_cmap="/home/bdmlab/hetero_test/multi/struct_dimer_cmaps/H1045TS285_3_chain_AB.cmap", _pred_cmap="/home/bdmlab/true/.cmap"))
# GLINTER_DIR ="/home/rsr3gt/anaconda3/envs/multi_eva/"
def glinter_runner(_first_pdb, _second_pdb, _out_dir, _is_homodimer, expected_cmaps_name, _glinter):
    GLINTER_DIR = _glinter
    envs = GLINTER_DIR + "scripts/set_env.sh"
    print(envs)
    os.system("export MKL_SERVICE_FORCE_INTEL=1")
    os.system("source " + str(envs))
    #    os.system("cd "+GLINTER_DIR)
    #    os.system("export MKL_SERVICE_FORCE_INTEL=1")
    name_1_list = os.path.basename(_first_pdb).split(".")[0]
    name_2_list = os.path.basename(_second_pdb).split(".")[0]
    #    os.system("cd " + GLINTER_DIR)
    if _is_homodimer == True:
        cmd = GLINTER_DIR + "/scripts/build_homo.sh " + str(_first_pdb) + " " + str(_second_pdb) + " " + str(
            _out_dir) + " " + str(name_2_list)
        print(os.system(cmd))
    else:
        cmd = GLINTER_DIR + "/scripts/build_hetero.sh " + str(_first_pdb.replace("//", "/")) + " " + str(
            _second_pdb.replace("//", "/")) + " " + str(_out_dir.replace("//", "/"))
        print(os.system(cmd))
    name = str(name_1_list) + ":" + str(name_2_list)
    cmap_file = _out_dir + name + "/score_mat.pkl"

    #    content =  np.load(cmap_file, allow_pickle=True)
    print(cmd)
    print(cmap_file)
    if os.path.exists(cmap_file):
        content = np.load(cmap_file, allow_pickle=True)
        _out_dir = _out_dir.replace("//", "/")
        dest_file = _out_dir.replace("extras/", "") + name + ".cmap"
        print(dest_file)
        np.savetxt(expected_cmaps_name, content)
        # cmd=  "cp "+cmap_file+" "+dest_file
    #        os.system(cmd)

    return


def cdpred_runner(_first_pdb, _second_pdb, _out_dir, _is_homodimer, expected_cmaps_name, _cdpred_dir):
    CDPRED_DIR = _cdpred_dir
    name_1_list = os.path.basename(_first_pdb).split(".")[0]
    name_2_list = os.path.basename(_second_pdb).split(".")[0]
    _name_full = str(name_1_list) + "_" + str(name_2_list)
    #    os.system("cd " + GLINTER_DIR)
    if _is_homodimer == True:
        # python lib/Model_predict.py -n T1084A_T1084B -p ./example/T1084A_T1084B.pdb -a ./example/T1084A_T1084B.a3m -m homodimer -o ./output/T1084A_T1084B/
        cmd = "sh " + CDPRED_DIR + " -n " + _name_full + " -p " + str(_first_pdb) + " -m homodimer -o " + str(
            _out_dir) + "/" + str(_name_full)
        print(os.system(cmd))
    else:
        # python lib/Model_predict.py -n H1017A_H1017B -p ./example/H1017A.pdb ./example/H1017B.pdb -a ./example/H1017A_H1017B.a3m -m heterodimer -o ./output/H1017A_H1017B/
        cmd = "sh " + CDPRED_DIR + " -n " + _name_full + " -p '" + str(_first_pdb) + " " + str(
            _second_pdb) + "' -m heterodimer -o " + str(_out_dir) + "/" + str(_name_full)
        print(os.system(cmd))

    cmap_file = _out_dir + _name_full + "/predmap/" + _name_full + ".htxt"

    cmd = cmd.replace("//", "/")
    cmap_file = cmap_file.replace("//", "/")
    print(cmd)
    print(cmap_file)
    if os.path.exists(cmap_file):
        content = np.loadtxt(cmap_file)
        _out_dir = _out_dir.replace("//", "/")
        dest_file = _out_dir.replace("extras/", "") + str(name_1_list) + ":" + str(name_2_list) + ".cmap"
        print(dest_file)
        np.savetxt(expected_cmaps_name, content)

    return


def check_single_exists(_file):
    if os.path.exists(_file):
        return True
    else:
        print("This file is missing" + str(_file))
        return False


def check_path_exists(_PARIWISE_QA_SCRIPT, _TM_SCORE_PATH, _Q_SCORE, _DOCK_Q_PATH, _MM_ALIGN, _GLINTER_DIR):
    print(_PARIWISE_QA_SCRIPT)
    print(_TM_SCORE_PATH)
    print(_Q_SCORE)
    print(_DOCK_Q_PATH)
    print(_MM_ALIGN)
    print(_GLINTER_DIR)
    if check_single_exists(_PARIWISE_QA_SCRIPT) and check_single_exists(_TM_SCORE_PATH) and check_single_exists(
            _Q_SCORE) and check_single_exists(_DOCK_Q_PATH) and check_single_exists(_MM_ALIGN) and check_single_exists(
        _GLINTER_DIR):
        return True
    else:
        return False


def save_mm_score(_pdb_dict, _file_name):
    mm_string = ""
    for values in copy.deepcopy(_pdb_dict):
        temp_string = str(values) + " " + str(_pdb_dict.get(values).multimer_scoring) + "\n"
        mm_string = mm_string + temp_string

    write2File(_file_name, mm_string)


def get_TM_score(_array):
    TM_Score_PATH = _array[2]
    _true = _array[0]
    _pred = _array[1]
    # contents = subprocess.check_output([MM_ALIGN_PATH, _true, _current])
    try:
        contents = subprocess.check_output([TM_Score_PATH, _pred, _true])
        for item in contents.decode("utf-8").split("\n"):
            if "TM-score    =" in item:
                tm_score_value = item.replace("TM-score    =", "").strip().split(" ")[0]
                return float(tm_score_value)
    except:
        return 0


def get_tm_score_parallel_submit(_array, _CPU_COUNT):
    all_value = []
    worker = int(_CPU_COUNT)
    with concurrent.futures.ThreadPoolExecutor(max_workers=worker) as executor:
        result_futures = list(map(lambda x: executor.submit(get_TM_score, x), _array))
        for future in concurrent.futures.as_completed(result_futures):
            try:
                all_value.append(future.result())
            except Exception as e:
                print('e is', e, type(e))
                all_value.append(0)

    return np.average(all_value)
