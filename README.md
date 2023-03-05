# MULTICOM_qa 
Evaluate the quality of protein multimer structure models. It uses the pairwise similarity between structural models and a deep learning-based interface contact probability score to predict the quality of the models. It was ranked first in predicting the global accuracy of strucural models of protein multimers/assemblies/complexes in the 15th Critical Assessment of Techniques for Protein Structure Prediction (CASP15) in 2022 (see the SCORE ranking at https://predictioncenter.org/casp15/zscores_EMA.cgi).

# Usage:

```python multi_eva_modified_parallel_lite.py <fasta_file in casp format> <dir of the structure for evaluation > <stoichiometry> <dir of the tertiary structures (Alphafold 2  recommended)> <CPU_count > <out put directory>```

```e.g python multi_eva_modified_parallel_lite.py /home/rsr3gt/programs/Multi_Eva/data/fasta_casp14/casp_capri_fasta/H1060.fasta /home/rsr3gt/programs/Multi_Eva/data/predictions_cleaned/H1060/ A6B3C12D6 /home/rsr3gt/programs/Multi_Eva/data/pdbs_casp_alphafold/H1060/ 20 /home/rsr3gt/programs/Multi_Eva/output/example_H1060_beta/```

Final output would be a csv file with the final score of all the structures in the output directory, e.g "H1060.csv"

Note the fasta file is expected in CAPRI format check the dir examples for clarification.
# Requirements:
1. Glinter (https://github.com/zw2x/glinter)
2. Qscore 
3. TMscore (https://zhanggroup.org/TM-score/)
4. DockQ (https://github.com/bjornwallner/DockQ)
5. MMalign (https://zhanggroup.org/MM-align/)
6. Biopython 
7. Futures
  
## Installation
1. Install the required softwares.
2. Open the ./config file and edit the path of the following variablse under the "Current section" with the appropiate path:
    1. PARIWISE_QA_SCRIPT = "path of the pairwise qa script"
    2. Q_SCORE = "path of the qscore script"
    3. TM_SCORE_PATH = "path of the tmscore"
    4. MM_ALIGN_PATH = "path of the MMalign tools"
    5. GLINTER_DIR = "path of the installed glinter directory"
    6. DOCK_Q_PATH = "path of the DOCKQ.py file"
   

##Current Version
To run the current and faster version of MultimerEva run the "multi_eva_modified_parallel_lite.py"


## Caution
Sometimes Glinter's env doesn't activate properly. So it is better to activate it before runnning and also if you encounter any error regarding "MKL_SERVICE" then use the command from below before running.

export MKL_SERVICE_FORCE_INTEL=1
