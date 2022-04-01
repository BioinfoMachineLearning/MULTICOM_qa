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
    true_mm_score = ''

    pass


final_results = []
mm_align_dir = "/home/bdmlab/label_csv/"
multi_eva_dir = "/home/bdmlab/csv/"
output_dir= "/home/bdmlab/"
list_of_file = file_reader("/home/bdmlab/dq_done.txt")
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
        temp.monomer_score = row[1].average_MS
        temp.ds_score = row[1].average_DS
        temp.icps_score = row[1].average_ICPS
        temp.recall = row[1].average_R
        temp.mm_align_score = row[1].average_MMS
        temp.final_score = row[1].Final_score
        eva_score_list[temp.name] = temp

    for values in data_align.iterrows():
        temp = "r-l_"+str(values[1].model)+"_tidy"
        eva_score_list.get(temp).true_mm_score = float(values[1].DcokQ)
        # print(float(values[1].DcokQ))
    print(_targets)
    # eva_score_list = sorted(eva_score_list.items(), key=lambda x: x[1].final,reverse=True)

    true_top_1_name = sorted(eva_score_list.items(), key=lambda x: x[1].true_mm_score, reverse=True)[0]
    true_top_1_score = true_top_1_name[1].true_mm_score
    # true_top_1_name = true_top_1_name[0]
    # print(true_top_1_name)
    # print(true_top_1_score)
    # print("here")


    ms_top_1_name = sorted(eva_score_list.items(), key=lambda x: x[1].monomer_score, reverse=True)[0]
    ms_top_1_score = ms_top_1_name[1].true_mm_score
    loss_ms_top_1_score = float(true_top_1_score)-float(ms_top_1_score)


    ds_top_1_name = sorted(eva_score_list.items(), key=lambda x: x[1].ds_score, reverse=True)[0]
    ds_top_1_score = ds_top_1_name[1].true_mm_score
    loss_ds_top_1_score = float(true_top_1_score)-float(ds_top_1_score)

    icps_top_1_name = sorted(eva_score_list.items(), key=lambda x: x[1].icps_score, reverse=True)[0]
    icps_1_score = icps_top_1_name[1].true_mm_score
    loss_icps_top_1_score = float(true_top_1_score)-float(icps_1_score)


    recall_top_1_name = sorted(eva_score_list.items(), key=lambda x: x[1].recall, reverse=True)[0]
    recall_top_1_score = recall_top_1_name[1].true_mm_score
    loss_recall_top_1_score = float(true_top_1_score)-float(recall_top_1_score)

    mm_align_top_1_name = sorted(eva_score_list.items(), key=lambda x: x[1].mm_align_score, reverse=True)[0]
    mm_align_top_1_score = mm_align_top_1_name[1].true_mm_score
    loss_mm_align_top_1_score = float(true_top_1_score)-float(mm_align_top_1_score)

    final_top_1_name = sorted(eva_score_list.items(), key=lambda x: x[1].final_score, reverse=True)[0]
    final_top_1_score = final_top_1_name[1].true_mm_score
    loss_final_top_1_score = float(true_top_1_score)-float(final_top_1_score)
    final_results.append(
        [Title.replace(".csv",""), ms_top_1_score, ds_top_1_score, icps_1_score, recall_top_1_score, mm_align_top_1_score, final_top_1_score,
         true_top_1_score, loss_ms_top_1_score,loss_ds_top_1_score,loss_icps_top_1_score,loss_recall_top_1_score,loss_mm_align_top_1_score,loss_final_top_1_score])

# print(final_results)

name_of_output_file = output_dir + str("DOCKQ_mm_eva") + ".csv"
_header_row = ['Target',"ms_top_1_score","ds_top_1_score","icps_1_score","recall_top_1_score","mm_align_top_1_score",'final_top_1_score', "true_top_1_score",'ms_loss','ds_loss','icps_loss','recall_loss','mm_align_loss','final_loss']
with open(name_of_output_file, 'w') as csvfile:
    filewriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    filewriter.writerow(_header_row)
    for data in final_results:
        filewriter.writerow(data)
