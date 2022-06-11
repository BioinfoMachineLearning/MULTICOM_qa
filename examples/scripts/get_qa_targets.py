import os
from subprocess import call
from datetime import datetime
import logging
import wget
import tarfile
from typing import List
from pathlib import Path
from urllib.error import HTTPError



father_path = Path(__file__).resolve().parents[1]
# DPROQ = f'{father_path}/DPROQ/src/evaluate_complex.py'



def txt_to_list(txt_file: str, pattern='\n') -> List:
	def remove_n(lst: List, pattern='\n') -> List:
		return [i.strip(pattern) for i in lst]



	with open(txt_file, 'r') as f:
		tmp_list = f.readlines()
	tmp_list = remove_n(tmp_list, pattern=pattern)
	return tmp_list




# read downloaded records
base_output_folder = '/bml/bml_casp15/casp15_qa_multimer_eva/test/auto_downloads/'
downloaded_record = txt_to_list(f'{base_output_folder}/downloaded_record.txt')
new_target = []



# global setting
now = datetime.now()
CURRENT_TIME = now.strftime("%m_%d")
output_folder = os.path.join(base_output_folder, CURRENT_TIME)
qa_decoys_tar_folder = f'{base_output_folder}/QA_decoys_tar'
qa_decoys_folder_stage1 = f'{base_output_folder}/QA_decoys/stage1'
qa_decoys_folder_stage2 = f'{base_output_folder}/QA_decoys/stage2'
os.makedirs(output_folder, exist_ok=True)
os.makedirs(qa_decoys_tar_folder, exist_ok=True)
os.makedirs(qa_decoys_folder_stage1, exist_ok=True)
os.makedirs(qa_decoys_folder_stage2, exist_ok=True)



# define logger
logging.basicConfig(filename=os.path.join('/bml/bml_casp15/casp15_qa_multimer_eva/test/auto_downloads/LOG', f'mm_eva_{CURRENT_TIME}.log'),
filemode='a',
format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S',
level=os.environ.get("LOGLEVEL", "INFO"))



# get target info
cmd = 'wget -r -nH --cut-dirs=2 --no-parent --reject="index.html*" http://sysbio.rnet.missouri.edu/multicom_cluster/QA/ > /dev/null 2>&1'
os.chdir(output_folder)
call([cmd], shell=True)



for i in os.listdir(output_folder):
	if i.startswith('CASP'):
		with open(os.path.join(output_folder, i)) as f:
			content = f.readlines()
		content = [i.strip('\n') for i in content]
		target_name, sequence, email, decoy_url, Date, Stoichiom = content

		# download decoy if it is a new target
		if decoy_url not in downloaded_record:
			try:
				wget.download(url=decoy_url, out=qa_decoys_tar_folder)
				logging.info(f'Downloaded {decoy_url}')
				downloaded_record.append(decoy_url)
				new_target.append(decoy_url.split('/')[-1])
			except HTTPError: # some download link is wrong
				continue



for i in os.listdir(qa_decoys_tar_folder):
	if i in new_target:
		tar = tarfile.open(name=os.path.join(qa_decoys_tar_folder, i), mode='r:gz')
	if 'stage1' in i: # need to change when offical run.
		tar.extractall(path=qa_decoys_folder_stage1)
		tar.close()
	else:
		tar.extractall(path=qa_decoys_folder_stage2)
		tar.close()



# update downloaded record
print('\n')
print(f'New targets are {new_target}')
if new_target:
	with open(f'{base_output_folder}/downloaded_record.txt', 'a') as f:
		for i in new_target:
			f.writelines(f'https://predictioncenter.org/download_area/CASP15/predictions/oligo/{i}' + '\n')
else:
	print('No new targets')