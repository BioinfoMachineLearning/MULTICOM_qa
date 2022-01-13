import copy
import os
import re


class multimer:
    unique_monomers_chain =[]
    chains = [],
    chain_fasta = {},
    stoi = [],
    # grp_id = ""
    scores = {}

    pass

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


def convert_to_pdb(_pdb, _name):
    content = ''
    for x in _pdb:
        content += correct_format(x) + '\n'
    f = open(_name, "w")
    f.write(content)
    f.close()
    return _pdb


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
    _pdb_copy.icode = space_returner(2 - len(_pdb_copy.chain)) + _pdb_copy.icode  # 27
    _pdb_copy.x = space_returner(8 - len(_pdb_copy.x)) + _pdb_copy.x  # 31-38
    _pdb_copy.y = space_returner(8 - len(_pdb_copy.y)) + _pdb_copy.y  # 39-46
    _pdb_copy.z = space_returner(8 - len(_pdb_copy.z)) + _pdb_copy.z  # 47-54
    _pdb_copy.occupancy = space_returner(6 - len(_pdb_copy.occupancy)) + _pdb_copy.occupancy  # 55-60
    _pdb_copy.temp_fact = space_returner(6 - len(_pdb_copy.temp_fact)) + _pdb_copy.temp_fact  # 61-66
    _pdb_copy.element = space_returner(4 - len(_pdb_copy.element)) + _pdb_copy.element  # 73-76
    _pdb_copy.charge = space_returner(2 - len(_pdb_copy.charge)) + _pdb_copy.charge  # 77-78
    content = _pdb_copy.atom + space_returner(3) + _pdb_copy.serial + space_returner(
        2) + _pdb_copy.atom_name + _pdb_copy.alt_loc + _pdb_copy.res_name + space_returner(
        1) + _pdb_copy.chain + _pdb_copy.res_num + _pdb_copy.icode + space_returner(
        3) + _pdb_copy.x + _pdb_copy.y + _pdb_copy.z + _pdb_copy.occupancy + _pdb_copy.temp_fact + space_returner(
        6) + _pdb_copy.element
    return content


def read_pdb(pdb):
    contents = []
    with open(pdb, "r") as f:
        for line in f:
            if (line.startswith("ATOM")):
                # pass
                contents.append(line)
    return contents


def pdb_from_array(_pdb, _filename):
    array = []
    content = ''
    for x in _pdb:
        val = string_array_from_pdb_array(x)
        array.append(val)
        content = content + val + '\n'
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
    _pdb_copy.icode = space_returner(2 - len(_pdb_copy.chain)) + _pdb_copy.icode  # 27
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


def pdb_from_array(_pdb, _filename):
    array = []
    content = ''
    for x in _pdb:
        val = string_array_from_pdb_array(x)
        array.append(val)
        content = content + val + '\n'
    f = open(_filename, "w")
    f.write(content + 'END')
    f.close()
    return array


def get_unique_chains(_inp_details):
    chain_array = []
    for val in _inp_details:
        chain_array.append(val.chain)
    return list(dict.fromkeys(chain_array))


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
        print(str(counter) + " " + values)
        if counter % 2 == 1:
            stoi_fasta_dict[temp_name] = values
        else:
            temp_name = str(int(counter / 2))
        counter = counter + 1

    return stoi_fasta_dict


def monomer_pdb_filtering(_pdb, _dir):
    tar_name = os.path.basename(_pdb)
    tar_dir = _dir + "/" + tar_name
    os.system("mkdir -p " + tar_dir)
    full_pdb = contents_to_info(read_pdb(_pdb))
    chain_finder = get_unique_chains(full_pdb)
    # print(chain_finder)
    for chain in chain_finder:
        temp_monomer_pdb = []
        temp_monomer_pdb = separate_by_chain(full_pdb, chain)

        tar_monomer_file = tar_dir + "/" + tar_name + "_chain_" + str(chain) + ".pdb"
        monomer_string = pdb_from_array(_pdb=temp_monomer_pdb, _filename=tar_monomer_file)
        # convert_to_pdb(_pdb=monomer_string,_name=tar_monomer_file)

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
                    print(_chains[0])
                    stoi_repeat_dict[val] = int(stoi_repeat_dict[val]) - 1
                    del _chains[0]
    a_multimer.chain_fasta = temp_chain_dict
    return a_multimer


# def fasta_to_chain_mapper(_fasta_file="/home/bdmlab/H1060.fasta",
#                           chains=["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q",
#                                   "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "a"], _stoi="A6B3C12D6"):
def fasta_to_chain_mapper(_fasta_file="/home/bdmlab/T1032o.fasta", _stoi="A2", _chains=["A", "B"]):
    fasta_stoic_dict = multi_fasta_reader(_seq_file=_fasta_file)
    chain_mapper = chain_pdb_combination_generator(_stoi, _chains, fasta_stoic_dict)

    return chain_mapper


# fasta_to_chain_mapper()
