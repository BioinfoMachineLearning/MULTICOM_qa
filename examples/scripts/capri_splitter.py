file_name = "/home/bdmlab/Downloads/Target192.pdb"
out_file = "/home/bdmlab/Downloads/capri_t192/"


def file_reader_remark(input_dir):
    contents = ""
    f = open(input_dir, "r")
    if f.mode == 'r':
        contents = f.read().splitlines()
        f.close()
    all_reamark = ""
    for remarks_line in contents:
        if "MD5" in remarks_line and "REMARK   9" in remarks_line and "MODEL" in remarks_line:
            all_reamark = all_reamark + remarks_line.replace("REMARK   9", "").replace("MODEL", "").replace("MD5",
                                                                                                            "").strip() + "\n"
        if "ATOM" in remarks_line:
            break
    model_dict = {}
    _temp_name = ""
    _temp_value = ""
    for line in all_reamark.splitlines():
        is_model_number = False
        _temp_name = ""
        _temp_value = ""
        _temp_line = line.split(" ")
        for value in _temp_line:
            if value != " " and is_model_number == False:
                _temp_name = value
                is_model_number = True
            elif value != " " and is_model_number == True:
                _temp_value = value

        model_dict[_temp_name] = _temp_value
    return model_dict


def write2File_list(_filename, _cont):
    print(_filename+"\n")
    _str_value = ""
    for values in _cont.splitlines():
        _str_value = _str_value + values + "\n"
    with open(_filename, "w") as f:
        f.writelines(_str_value.strip())
        f.close()


def file_reader_model(input_dir, _dict, _out):
    contents = ""
    f = open(input_dir, "r")
    if f.mode == 'r':
        contents = f.read().splitlines()
        f.close()
    all_atom = ""
    for _atom in contents:
        if "ATOM" in _atom or _atom.startswith("MODEL") or _atom.endswith("ENDMDL"):
            all_atom = all_atom + _atom + "\n"

    only_atoms = all_atom.split("ENDMDL")

    for value in only_atoms:
        name = ""
        full_model = ""
        for files in value.splitlines():
            if files.startswith("MODEL"):
                name = files.replace("MODEL", "").strip()
        temp_out_file = _out + name
        if len(name.strip()) > 0:
            write2File_list(temp_out_file, value)

    return None


remark_file = file_reader_remark(file_name)
# print(remark_file)
file_reader_model(file_name, remark_file, out_file)
