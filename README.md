GATE
====

### Genomic Annotation from Time-course Epigenomic data 

#### 1. Preprocessing of time-series data

Use the script "Preprocessing\_selectRegion.py"

  help:
```
usage: Preprocessing_selectRegion.py [-h] [-i INFO_BED] [-g GENOME]
                                     [-p PARAMETERS [PARAMETERS ...]]
                                     [-n NORMALIZE] [-o OUTPUT]

Data preprocessing and region selction for GATE

optional arguments:
  -h, --help            show this help message and exit
  -i INFO_BED, --info_bed INFO_BED
                        information of aligned bed files to be inputed,must be
                        stored in the folder "bed_file"
  -g GENOME, --genome GENOME
                        genome information for input data, can be
                        mm9/hg19/hg18/
  -p PARAMETERS [PARAMETERS ...], --parameters PARAMETERS [PARAMETERS ...]
                        parameters [m,n] to select regions. those regions with
                        normalized signal no less than 'm' in at least one
                        time point for at least 'n' epi-marks were selected,
                        default:[5,5]
  -n NORMALIZE, --normalize NORMALIZE
                        normalize total read counts into this number,
                        default:1E8
  -o OUTPUT, --output OUTPUT
                        output txt file
None
```
The files in folder "/bed_file/" should be the ones with names in the third column of "bed_info.txt".
The output file of this script can be the [input for GATE software](http://systemsbio.ucsd.edu/GATE/Input.htm).

#### 2. Run your data using GATE model
Use the output of step 1 as input file and run the program following [the steps here](http://systemsbio.ucsd.edu/GATE/Usage.htm).
