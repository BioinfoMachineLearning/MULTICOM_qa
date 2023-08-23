# MULTICOM_qa 
Estimate (predict) the quality of protein multimer structure models. It uses the pairwise similarity between structural models and a deep learning-based interface contact probability score to predict the quality of the models. It was ranked first in predicting the global accuracy of strucural models of protein multimers/assemblies/complexes in the 15th Critical Assessment of Techniques for Protein Structure Prediction (CASP15) in 2022 (see the SCORE ranking at https://predictioncenter.org/casp15/zscores_EMA.cgi).



<p align="center">
  <img alt="not found"  width="60%" height="20%" src="/examples/asset/mqa.gif">
</p>


# Requirements:
1. [CDPred](https://github.com/BioinfoMachineLearning/CDPred)
2. [MMalign](https://zhanggroup.org/MM-align/)
3. [Biopython](https://biopython.org/)
4. [Futures](https://docs.python.org/3/library/concurrent.futures.html)
5. [Numpy](https://numpy.org/install/)
  
## Installation
1. Download Multicom_qa
    
    git clone https://github.com/BioinfoMachineLearning/MULTICOM_qa.git


2. Download and setup the CDPred.


3. Activate the  CDPred enviroment.


4. Pip install -r requirements.txt.


5. Copy the  file "run_CDFold.sh" inside the CDPred directory


6. Download the MMalign 


7. Open the "config" file and edit the path of the following variables with the appropiate path:

    1. MM_ALIGN_PATH = "path of the MMalign tools"
        
       e.g.   MM_ALIGN_PATH = "/home/multicom4s_tool/MMalign"
   
    2. CDPred = "path of the installed CDPred directory" 

       e.g.   CDPRED_PATH = "/home/multicom4s_tool/CDpred_mm_eva/CDPred/run_CDFold.sh
   

# Usage

```python multi_eva.py <fasta_file in casp format> <dir of the structure for evaluation > <stoichiometry> <dir of the tertiary structures (Alphafold 2  recommended)> <CPU_count > <out put directory>```

```e.g python multi_eva.py ./examples/H1143/H1143.fasta ./examples/H1143/decoys/ A1B1 ./examples/H1143/H1143_alphafold_monomer/ 4 ./examples/H1143/output/```

Final output would be a csv file with the final score of all the structures in the output directory, e.g "H1143.csv". The csv will contain "Final Rank" which is the weighted average of the Interchain Contact Probability Score (ICPS) rank and Pairwise Similarity Score (PSS) RANK.  And will also contain "Final_Score" which is the weighted average of the ICPS and PSS. The expected final output is provided in the expected_output directory.   



### Predicted score vs True score



The predicted result coincides with true one as displayed below. 

| Model  | Final Rank | Final Score | True TM-Score|
| ------------- | ---: | ---:  | ---:  |
|H1143TS035_5 | 1.4  |0.4712|0.9528 |
|H1143TS035_3 | 1.6  |0.4676|0.8534 |
|H1143TS035_4 | 3    |0.4227 |0.6771 |
|H1143TS399_4 | 4  |0.41186 |0.6509 |
|H1143TS037_3 | 5.4  |0.36549 |0.5926 |
|H1143TS493_2 | 5.6  |0.36231 |0.5919 |



# Reference


[Roy, R. S., Liu, J., Giri, N., Guo, Z., & Cheng, J. (2023). Combining pairwise structural similarity and deep learning interface contact prediction to estimate protein complex model accuracy in CASP15. bioRxiv.](https://doi.org/10.1101/2023.03.08.531814)