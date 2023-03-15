# MULTICOM_qa 
Estimate (predict) the quality of protein multimer structure models. It uses the pairwise similarity between structural models and a deep learning-based interface contact probability score to predict the quality of the models. It was ranked first in predicting the global accuracy of strucural models of protein multimers/assemblies/complexes in the 15th Critical Assessment of Techniques for Protein Structure Prediction (CASP15) in 2022 (see the SCORE ranking at https://predictioncenter.org/casp15/zscores_EMA.cgi).



# Requirements:
1. CDPred (https://github.com/BioinfoMachineLearning/CDPred)
5. MMalign (https://zhanggroup.org/MM-align/)
6. Biopython (https://biopython.org/)
7. Futures (https://docs.python.org/3/library/concurrent.futures.html)
  
## Installation
1. Install the required softwares.
2. Open the ./config file and edit the path of the following variables with the appropiate path:
    1. MM_ALIGN_PATH = "path of the MMalign tools"
    2. CDPred = "path of the installed glinter directory"
   

# Usage:

```python multi_eva_modified_parallel_lite.py <fasta_file in casp format> <dir of the structure for evaluation > <stoichiometry> <dir of the tertiary structures (Alphafold 2  recommended)> <CPU_count > <out put directory>```

```e.g python multi_eva_modified_parallel_lite.py /home/rsr3gt/programs/Multi_Eva/data/fasta_casp14/casp_capri_fasta/H1060.fasta /home/rsr3gt/programs/Multi_Eva/data/predictions_cleaned/H1060/ A6B3C12D6 /home/rsr3gt/programs/Multi_Eva/data/pdbs_casp_alphafold/H1060/ 20 /home/rsr3gt/programs/Multi_Eva/output/example_H1060_beta/```

Final output would be a csv file with the final score of all the structures in the output directory, e.g "H1060.csv"

Note: the fasta file is expected in CAPRI format check the dir examples for clarification.

