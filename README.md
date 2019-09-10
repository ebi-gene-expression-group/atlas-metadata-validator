# atlas-metadata-validator

## Single-cell MAGE-TAB validator

A MAGE-TAB pre-validation module for running checks that guarantee the experiment can be processed for [Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/home). The checks are mainly covering the pre-validation by https://github.com/ebi-gene-expression-group/scxa-control-workflow/blob/master/bin/sdrfToNfConf.R in order to guarantee correct processing. It reads metadata directly from the MAGE-TAB and generates a log file in the directory of the IDF file.

The checks can be invoked using the atlas_validation script with an IDF file path as input:
```
python atlas_validation.py path/to/test.idf.txt 
```

### Options
- The SDRF file is expected in the same directory as the IDF file. If this is not the case, the location of the SDRF and other data files can be specified with `-d PATH_TO_DATA` option. Provide the file path after the 
- The script guesses the experiment type (sequencing, microarray or single-cell) from the MAGE-TAB. If this was unsuccessful the experiment type can be set by specifying the respective argument `-seq`, `-ma` or `-sc`. 
- The data file and URI checks may take long time. Hence there is an option to skip these checks with `-x`.
- Verbose logging can be activated with `-v`.