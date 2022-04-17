# perl pairwise_model_eva.pl /home/rajroy/multi_eva_test/monomer_chains/A/ /home/rajroy/multi_eva_test/Multimet_evatest_samples/casp_fasta/H1036A.fasta /home/rajroy/pairwiseQA/q_score /home/rajroy/Downloads/tools/TMscore A /home/rajroy/q_A/
class config:
    DEVICE = "Current"
    PARIWISE_QA_SCRIPT = ""
    Q_SCORE = ""
    TM_SCORE_PATH = ""
    if "lab_pc" == DEVICE:
        PARIWISE_QA_SCRIPT = "/home/rajroy/pairwise_model_eva.pl"
        Q_SCORE = "/home/rajroy/pairwiseQA/q_score"
        TM_SCORE_PATH = "/home/rajroy/Downloads/tools/TMscore"
        MM_ALIGN_PATH = "/home/bdmlab/Documents/tools/MMalign"
    elif "home_pc" == DEVICE:
        print("home_pc")
        PARIWISE_QA_SCRIPT = "/home/bdmlab/tools/pairwise_model_eva.pl"
        Q_SCORE = "/home/bdmlab/tools/pairwiseQA/q_score"
        TM_SCORE_PATH = "/home/bdmlab/tools/TMscore"
        MM_ALIGN_PATH = "/home/bdmlab/tools/MMalign"
        GLINTER_DIR = "/home/bdmlab/anaconda3/envs/multi_eva/"
        DOCK_Q_PATH = "/home/bdmlab/tools/DockQ/DockQ.py"
    elif "Current" == DEVICE:
        PARIWISE_QA_SCRIPT = "/home/multicom4s_tool/pairwise_model_eva.pl"
        Q_SCORE = "/home/multicom4s_tool/pairwiseQA/q_score"
        TM_SCORE_PATH = "/home/multicom4s_tool/TMscore"
        MM_ALIGN_PATH = "/home/multicom4s_tool/MMalign"
        GLINTER_DIR = "/home/rsr3gt/programs/glinter/"
        DOCK_Q_PATH = "/home/multicom4s_tool/DockQ/DockQ.py"
        CDPRED_PATH="/home/multicom4s_tool/CDpred_mm_eva/CDPred/run_CDFold.sh"
