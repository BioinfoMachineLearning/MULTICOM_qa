# perl pairwise_model_eva.pl /home/rajroy/multi_eva_test/monomer_chains/A/ /home/rajroy/multi_eva_test/Multimet_evatest_samples/casp_fasta/H1036A.fasta /home/rajroy/pairwiseQA/q_score /home/rajroy/Downloads/tools/TMscore A /home/rajroy/q_A/
class config:
    DEVICE = "home_pc"
    PARIWISE_QA_SCRIPT = ""
    Q_SCORE = ""
    TM_SCORE_PATH = ""
    if "lab_pc" == DEVICE:
        PARIWISE_QA_SCRIPT = "/home/rajroy/pairwise_model_eva.pl"
        Q_SCORE = "/home/rajroy/pairwiseQA/q_score"
        TM_SCORE_PATH = "/home/rajroy/Downloads/tools/TMscore"
        MM_ALIGN_PATH = "/home/bdmlab/Documents/tools/MMalign"
    elif "home_pc" == DEVICE:
        PARIWISE_QA_SCRIPT = "/home/bdmlab/pairwise_model_eva.pl"
        Q_SCORE = "/home/bdmlab/pairwiseQA/q_score"
        TM_SCORE_PATH = "/home/bdmlab/Documents/tools/TMscore"
        MM_ALIGN_PATH = "/home/bdmlab/Documents/tools/MMalign"

