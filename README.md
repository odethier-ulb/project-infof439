# INFOF439 - Caretta
Reimplementation of Caretta ([researchgate page](https://www.researchgate.net/publication/340479604_Caretta_-_A_Multiple_Protein_Structure_Alignment_and_Feature_Extraction_Suite)) 
for the project of the course INFOF439 by <i><b>ODT</b></i> and <i><b>TI</b></i>.

## Examples

N.B. By default, the alignment results will be written in 

> caretta/alignment_results

### To see the help:

```
cd caretta
python ./caretta.py --help
```

### To generate the msa in Fig. 3. A in the original paper

```
python ./caretta.py ../test_data/matt_caretta --extension=.atm
```

The normal execution should output something like:
```
2022-05-11 18:45:33,050 | INFO: Found 5 proteins to align
2022-05-11 18:45:33,055 | INFO: Guide tree construction...
2022-05-11 18:45:34,552 | INFO: Protein alignment...
2022-05-11 18:45:35,277 | INFO: Done, the alignment results have been written to alignment_results
Computed RMSD : 3.835492856349312 A
```

### To generate the msa in Fig. 3. B in the original paper

```
python ./caretta.py ../test_data/mtm_caretta --extension=.atm
```
