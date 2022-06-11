"""
@ Description: Split CAPRI file, and md5 file
"""

import os
from argparse import ArgumentParser

parser = ArgumentParser(description='Generate label file')
parser.add_argument('--capri_file', '-c', type=str, required=True)
parser.add_argument('--output_folder', '-o', type=str, required=True)

args = parser.parse_args()
capri_file = args.capri_file
output_folder = args.output_folder

capri_file = os.path.abspath(capri_file)

if not os.path.isfile(capri_file):
    raise FileNotFoundError(capri_file)

target_id = capri_file.split('/')[-1].split('.')[0]

output_folder = os.path.abspath(output_folder)
output_folder = os.path.join(output_folder, target_id)
model_folder = os.path.join(output_folder, 'model')
md5_folder = os.path.join(output_folder, 'md5')

if not os.path.isdir(output_folder):
    os.makedirs(output_folder, exist_ok=True)

os.makedirs(model_folder, exist_ok=True)
os.makedirs(md5_folder, exist_ok=True)

with open(capri_file) as f:
    contents = f.readlines()
f.close()

MD5_record = []
MODEL_record =[]


for i in contents:
    if i.startswith('REMARK   9 MODEL'):
        MD5_record.append(i)
    if i.startswith('MODEL') or i.startswith('ATOM') or i.startswith('TER') or i.startswith('ENDMDL'):
        MODEL_record.append(i)

for i in MD5_record:
    model_id = i.split(' ')[-3].strip()
    with open(f'{md5_folder}/{model_id}.md5', 'w') as f:
        f.writelines(i.strip())

for i in MODEL_record:
    if i.startswith('MODEL'):
        model_id = i.split(' ')[-1].strip()
        print(model_id)
        new_file = f'{model_folder}/{model_id}.pdb'
        f = open(new_file, 'w')
    f.writelines(i.strip() + '\n')