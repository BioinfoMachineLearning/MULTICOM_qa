# MultimerEva 
Evaluate the quality of protein multimer structure models

# Usage:

```python multi_eva_modified_parallel.py <fasta_file in casp format> <dir of the structure for evaluation > <stoichiometry> <dir of the tertiary structures (Alphafold 2  recommended)> <CPU_count > <out put directory>```

```e.g python multi_eva_modified_parallel.py /home/rsr3gt/programs/Multi_Eva/data/fasta_casp14/casp_capri_fasta/H1060.fasta /home/rsr3gt/programs/Multi_Eva/data/predictions_cleaned/H1060/ A6B3C12D6 /home/rsr3gt/programs/Multi_Eva/data/pdbs_casp_alphafold/H1060/ 20 /home/rsr3gt/programs/Multi_Eva/output/example_H1060_beta/```

Final output would be a csv file with the final score of all the structures in the output directory, e.g "H1060.csv"

# Requirements:
1. Glinter (https://github.com/zw2x/glinter)
2. Qscore 
3. TMscore (https://zhanggroup.org/TM-score/)
4. DockQ (https://github.com/bjornwallner/DockQ)
5. MMalign (https://zhanggroup.org/MM-align/)
6. Biopython 
7. Futures
  
## First Time Setup
Change the path of the tools/softwares accroding to your system in the config file.

## Caution
Sometimes Glinter's env doesn't activate properly. So it is better to activate it before runnning and also if you encounter any error regarding "MKL_SERVICE" then use the command from below before running.

export MKL_SERVICE_FORCE_INTEL=1
(READNE NEEDS UPDATES)
