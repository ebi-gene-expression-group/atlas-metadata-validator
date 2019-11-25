#!/usr/bin/env python3

"""
This script takes an IDF file as input and runs validation of the metadata in the common datamodel.
"""

import argparse
import os
import sys

from atlas_metadata_validator.file import create_logger, file_exists, is_utf8
from atlas_metadata_validator.parser import get_sdrf_path, guess_submission_type_from_sdrf, guess_submission_type_from_idf, \
    read_sdrf_file, simple_idf_parser

from atlas_metadata_validator.atlas_checker import AtlasMAGETABChecker


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('idf',
                        help="Path to the MAGE-TAB IDF file")
    parser.add_argument('-d', "--data_dir", default="",
                        help="Path to the directory with SDRF and data files")
    parser.add_argument('-v', '--verbose', action='store_const', const=10, default=20,
                        help="Option to output detailed logging (debug level).")
    parser.add_argument('-hca', action='store_true', default=False,
                        help="Mark experiment as HCA import")
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('-sc', '--singlecell', action='store_const', const="singlecell", dest='submission_type',
                       help="Force submission type to be 'singlecell'")
    group.add_argument('-seq', '--sequencing', action='store_const', const="sequencing", dest='submission_type',
                       help="Force submission type to be 'sequencing'")
    group.add_argument('-ma', '--microarray', action='store_const', const="microarray", dest='submission_type',
                       help="Force submission type to be 'microarray'")
    parser.add_argument('-x', '--skip-file-checks', action='store_true',
                        help="Skip file and URI checks")
    args = parser.parse_args()

    return args


def main():
    process_name = "atlas_validation"

    args = parse_args()
    idf_file, data_dir, logging_level = args.idf, args.data_dir, args.verbose
    submission_type = args.submission_type
    
    # Exit if IDF file doesn't exist or isn't in UTF-8 encoding
    file_exists(idf_file)
    is_utf8(idf_file)

    # Create logger
    current_dir, idf_file_name = os.path.split(idf_file)
    logger = create_logger(current_dir, process_name, idf_file_name, logger_name="ATLAS", log_level=logging_level)

    # Print input file name for clarity when doing multiple validations
    logger.info("Validating {} and associated SDRF file".format(idf_file))

    # Get path to SDRF file
    sdrf_file_path = get_sdrf_path(idf_file, logger, data_dir)

    # Read IDF/SDRF and get submission type
    idf_dict = simple_idf_parser(idf_file)
    sdrf_data, header, header_dict = read_sdrf_file(sdrf_file_path)

    # Set submission type
    if not submission_type:
        submission_type = guess_submission_type_from_sdrf(sdrf_data, header, header_dict)
        if not submission_type:
            submission_type = guess_submission_type_from_idf(idf_dict)
        logger.info("Detected submission type: {}".format(submission_type))
    else:
        logger.info("Setting submission type to \"{}\"".format(submission_type))

    atlas_checker = AtlasMAGETABChecker(idf_file, sdrf_file_path, submission_type,
                                        skip_file_checks=args.skip_file_checks,
                                        is_hca=args.hca)
    atlas_checker.check_all(logger)

    # Collect error codes
    error_codes = atlas_checker.errors

    if error_codes:
        logger.info("Validation finished with the following error codes: {}".format(", ".join(error_codes)))
        sys.exit(1)
    else:
        logger.info("Validation was successful! 🍰")


if __name__ == '__main__':
    main()
