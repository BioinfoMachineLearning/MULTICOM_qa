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
    true_mm_score = ''

    pass


final_results = []
mm_align_dir = "/home/bdmlab/MM_Eva/casp14/casp_14mmalign_results/"
multi_eva_dir = "/home/bdmlab/CDPRED_csv/"
output_dir= "/home/bdmlab/"
list_of_file = file_reader("/home/bdmlab/CDPRED_csv/list.txt")
for _targets in list_of_file:
    Title = str(_targets)+".csv"

    temp_result_file_mm_align_file =mm_align_dir + str(Title)
    data_align = pd.read_csv(temp_result_file_mm_align_file)

    temp_result_file_mm_eva_file =multi_eva_dir + str(Title)
    data_eva = pd.read_csv(temp_result_file_mm_eva_file)

    eva_score_list = {}
    for row in data_eva.iterrows():
        temp = eva_score()
        temp.name = row[1].Name
        # temp.monomer_score = row[1].average_MS
        # temp.ds_score = row[1].average_DS
        temp.icps_score = row[1].average_ICPS
        # temp.recall = row[1].average_R
        temp.mm_align_score = row[1].average_MMS
        temp.final_score = row[1].Final_score
        # temp.new_final_score = row[1].New_Score
        temp.new_final_rank = row[1].Final_Rank
        eva_score_list[temp.name] = temp

    for values in data_align.iterrows():
        temp = values[1].Target
        eva_score_list.get(temp).true_mm_score = values[1].MMalign_Score
    # eva_score_list = sorted(eva_score_list.items(), key=lambda x: x[1].final,reverse=True)
    true_top_1_name = sorted(eva_score_list.items(), key=lambda x: x[1].true_mm_score, reverse=True)[0]
    true_top_1_score = true_top_1_name[1].true_mm_score
    # true_top_1_name = true_top_1_name[0]
    # print(true_top_1_name)
    # print(true_top_1_score)
    # print("here")

    #
    # ms_top_1_name = sorted(eva_score_list.items(), key=lambda x: x[1].monomer_score, reverse=True)[0]
    # ms_top_1_score = ms_top_1_name[1].true_mm_score
    # loss_ms_top_1_score = float(true_top_1_score)-float(ms_top_1_score)

    #
    # ds_top_1_name = sorted(eva_score_list.items(), key=lambda x: x[1].ds_score, reverse=True)[0]
    # ds_top_1_score = ds_top_1_name[1].true_mm_score
    # loss_ds_top_1_score = float(true_top_1_score)-float(ds_top_1_score)

    icps_top_1_name = sorted(eva_score_list.items(), key=lambda x: x[1].icps_score, reverse=True)[0]
    icps_1_score = icps_top_1_name[1].true_mm_score
    loss_icps_top_1_score = float(true_top_1_score)-float(icps_1_score)

    #
    # recall_top_1_name = sorted(eva_score_list.items(), key=lambda x: x[1].recall, reverse=True)[0]
    # recall_top_1_score = recall_top_1_name[1].true_mm_score
    # loss_recall_top_1_score = float(true_top_1_score)-float(recall_top_1_score)

    mm_align_top_1_name = sorted(eva_score_list.items(), key=lambda x: x[1].mm_align_score, reverse=True)[0]
    mm_align_top_1_score = mm_align_top_1_name[1].true_mm_score
    loss_mm_align_top_1_score = float(true_top_1_score)-float(mm_align_top_1_score)

    # new_score_top_1_name = sorted(eva_score_list.items(), key=lambda x: x[1].new_final_score, reverse=True)[0]
    # new_score_top_1_score = new_score_top_1_name[1].true_mm_score
    # loss_new_score_top_1_score = float(true_top_1_score) - float(new_score_top_1_score)

    new_rank_top_1_name = sorted(eva_score_list.items(), key=lambda x: x[1].new_final_rank)[0]
    new_rank_top_1_score = new_rank_top_1_name[1].true_mm_score
    loss_new_rank_top_1_score = float(true_top_1_score) - float(new_rank_top_1_score)


    final_top_1_name = sorted(eva_score_list.items(), key=lambda x: x[1].final_score, reverse=True)[0]
    final_top_1_score = final_top_1_name[1].true_mm_score
    loss_final_top_1_score = float(true_top_1_score)-float(final_top_1_score)
    final_results.append(
        [Title.replace(".csv",""), icps_1_score, mm_align_top_1_score, new_rank_top_1_score, final_top_1_score,new_rank_top_1_score,
         true_top_1_score,loss_icps_top_1_score,loss_mm_align_top_1_score,loss_final_top_1_score,loss_new_rank_top_1_score])

# print(final_results)

name_of_output_file = output_dir + str("casp_14_CDPRED_mm_eva_040722") + ".csv"
_header_row = ['Target',"icps_1_score","mm_align_top_1_score",'final_top_1_score','new_rank_top_1_score', "true_top_1_score",'icps_loss','mm_align_loss','final_loss','new_rank_top']
with open(name_of_output_file, 'w') as csvfile:
    filewriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    filewriter.writerow(_header_row)
    for data in final_results:
        filewriter.writerow(data)
