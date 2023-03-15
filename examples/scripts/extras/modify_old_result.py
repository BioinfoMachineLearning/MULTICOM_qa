import csv

import pandas as pd


def file_reader(input_dir):
    contents = ""
    f = open(input_dir, "r")
    if f.mode == 'r':
        contents = f.read().splitlines()
        f.close()
    return contents
class eva_score:
    name = ''
    monomer_score = ''
    ds_score = ''
    icps_score = ''
    recall = ''
    mm_align_score = ''
    final_score = ''
    new_final_score = ''
    new_final_rank=''
    icps_rank = 0
    mm_rank = 0

    pass
multi_eva_dir = "/home/bdmlab/MM_Eva/DG/mm_eva_dockG_results/"
output_dir= "/home/bdmlab/MM_Eva/DG/modified_output/"

list_of_file = file_reader("/home/bdmlab/MM_Eva/DG/mm_eva_dockG_results/list_done.txt")

for pro_values in list_of_file:
    temp_result_file_mm_eva_file = multi_eva_dir+"/"+pro_values+".csv"
    data_eva = pd.read_csv(temp_result_file_mm_eva_file)

    eva_score_list = {}
    for row in data_eva.iterrows():
        temp = eva_score()
        temp.name = row[1].Name
        temp.monomer_score = row[1].average_MS
        temp.ds_score = row[1].average_DS
        temp.icps_score = row[1].average_ICPS
        temp.recall = row[1].average_R
        temp.mm_align_score = row[1].average_MMS
        temp.final_score = row[1].Final_score
        # temp.new_final_score = row[1].New_Score
        # temp.new_final_rank = row[1].W_Rank
        eva_score_list[temp.name] = temp

    rank_icps_1_name = sorted(eva_score_list.items(), key=lambda x: x[1].icps_score, reverse=True)
    rank_mm_1_name = sorted(eva_score_list.items(), key=lambda x: x[1].mm_align_score, reverse=True)
    counter =1
    for values in rank_icps_1_name:
        eva_score_list[values[1].name].icps_rank = counter
        # print(values[1].icps_score)
        counter = counter +1
        # ms_top_1_score = ms_top_1_name[1].true_mm_score
    counter =1
    for values in rank_mm_1_name:
        eva_score_list[values[1].name].mm_rank = counter
        # print(values[1].mm_align_score)
        counter = counter +1
        # ms_top_1_score = ms_top_1_name[1].tr


    for values in eva_score_list:
        temp = eva_score_list.get(values)
        temp.new_final_score = temp.icps_score  * 0.4 + 0.6 *temp.mm_align_score
        temp.new_final_rank = temp.icps_rank  * 0.4 + 0.6 *temp.mm_rank
        # print(values[1])
        # print("here")
    _header_row = ["Name",	"average_MS",	"average_DS",	"average_ICPS",	"average_R",	"average_MMS",	"Final_score", "New_Score",	"ICPS_Rank","MM_RANK","FINAL_RANK"]
    with open(output_dir+"/"+pro_values+".csv", 'w') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(_header_row)
        for data in eva_score_list:
            temp = eva_score_list.get(data)
            filewriter.writerow([data,temp.monomer_score,temp.ds_score,temp.icps_score,temp.recall,temp.mm_align_score,temp.final_score,temp.new_final_score,temp.icps_rank,temp.mm_rank,temp.new_final_rank])


