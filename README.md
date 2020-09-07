# Expression Atlas metadata validator

This is a Python module to parse a set of MAGE-TAB files and check for compatibility with the Expression Atlas and Single Cell Expression Atlas analysis pipelines. The main validation script automatically detects the experiment type from the MAGE-TAB and runs the respective tests. Currently general checks (for bulk and single-cell experiment) as well as Single Cell Expression Atlas specific checks are supported. The controlled vocabulary for the checks are stored in the [metadata-validation-config](https://github.com/ebi-gene-expression-group/metadata-validation-config) repository. 

## Requirements

Python 3.6, requests (tested with version 2.20.1), and GitPython (tested with version 3.1.7)


## Install

Install with pip:
```
pip install atlas-metadata-validator
```

Set environment variables to fetch the validation configuration for controlled vocabulary terms:

```bash
export VALIDATION_CONFIG_REPO='https://github.com/ebi-gene-expression-group/metadata-validation-config'
export VALIDATION_CONFIG_FILE='atlas_validation_config.json'
```
(if not specified the above values are used as default)

Set local path where to store the validation config files, default is working dir when calling the script
```bash
export VALIDATION_CONFIG_LOCAL_PATH='/path/to/your/local/copy/'
```


## Single-cell MAGE-TAB validator

A MAGE-TAB pre-validation module for running checks that guarantee the experiment can be processed for [Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/home). The checks are mainly covering the pre-validation by https://github.com/ebi-gene-expression-group/scxa-control-workflow/blob/master/bin/sdrfToNfConf.R in order to guarantee correct processing. It reads metadata directly from the MAGE-TAB and generates a log file in the directory of the IDF file.

The checks can be invoked using the atlas_validation script with an IDF file path as input:
```
atlas_validation.py path/to/test.idf.txt 
```

### Options
- The SDRF file is expected in the same directory as the IDF file. If this is not the case, the location of the SDRF and other data files can be specified with `-d PATH_TO_DATA` option.
- The script guesses the experiment type (sequencing, microarray or single-cell) from the MAGE-TAB. If this was unsuccessful the experiment type can be set by specifying the respective argument `-seq`, `-ma` or `-sc`. 
- The data file and URI checks may take long time. Hence there is an option to skip these checks with `-x`.
- Verbose logging can be activated with `-v`.
- Special validation rules for HCA-imported experiments can be invoked with `-hca` option. The validator will otherwise guess if the experiment is an HCA import based on the HCAD accession code in the ExpressionAtlasAccession field. 
