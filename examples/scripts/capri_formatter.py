import copy

result_file = "/home/bdmlab/T1110_casp_format.txt"
capri_dir = "/home/bdmlab/shawn_process/Target193/"

output_file = "/home/bdmlab/capri_target_193_raj.txt"


def pdb_file_reader(input_dir):
    # reads pdb
    _input_dir = copy.deepcopy(input_dir)
    contents = ""
    f = open(_input_dir, "r")
    if f.mode == 'r':
        contents = f.read().splitlines()
        f.close()
    file_str =""
    for _lines in contents[1:]:
        file_str = file_str+_lines+"\n"
    return file_str.strip()

def file_reader(input_dir):
    # reads pdb
    _input_dir = copy.deepcopy(input_dir)
    contents = ""
    f = open(_input_dir, "r")
    if f.mode == 'r':
        contents = f.read().splitlines()
        f.close()
    return contents


def write2File(_filename, _cont):
    with open(_filename, "w") as f:
        f.writelines(_cont)
        f.close()

def get_md5_code (_file):
    md5_file = file_reader(_file)
    md5_code = md5_file[0].split("MD5")[1].strip()
    return md5_code

def space_returner(_input):
    i = 0
    space = ""
    while i < _input:
        space = space + " "
        i = i + 1
    return space

def capri_formatter_header (_array,_capri_dir):
    model_file = capri_dir + "model/"
    md5_file = capri_dir + "md5/"
    hdr_str = "HEADER    TEMPLATE FOR CAPRI ROUND 50 TARGET 193"+"\n"
    hdr_str = hdr_str+ "COMPND    MOL_ID: 3;"+"\n"
    hdr_str = hdr_str+ "COMPND   2    MOLECULE: T1110;"+"\n"
    hdr_str = hdr_str+ "COMPND   3    CHAIN: A, B;"+"\n"
    hdr_str = hdr_str+"AUTHOR: Raj S.Roy , Jianlin Cheng"+"\n"

    for _models in _array:
        hdr_str = hdr_str+ "REMARK   9      MODEL"+space_returner(4-len(_models))+str(_models)+" MD5 "+get_md5_code(md5_file+str(_models)+".md5")+"\n"
    counter = 1
    for _pdbs in _array:
        hdr_str =hdr_str+"MODEL        "+str(counter)+"\n"
        hdr_str= hdr_str+pdb_file_reader(model_file+str(_pdbs)+".pdb")+"\n"
        counter = counter+1
    hdr_str =hdr_str+"END"+"\n"

    return hdr_str



res_str_file = file_reader(result_file)[6:16]
top_10_model = []
for value in res_str_file:
    top_10_model.append(value.split(" ")[0].strip())

final_capri_result = capri_formatter_header(top_10_model,capri_dir)
write2File(_filename=output_file,_cont=final_capri_result)
# print(top_10_model)

