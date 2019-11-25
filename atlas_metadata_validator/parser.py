
import codecs
import os
import re

from collections import OrderedDict, defaultdict

from atlas_metadata_validator.fetch import get_controlled_vocabulary
from atlas_metadata_validator.file import file_exists, is_utf8

SDRF_FILE_NAME_REGEX = r"^\s*SDRF\s*File"
DEFAULT_DATA_DIRECTORY = "unpacked"


def simple_idf_parser(idf_file):
    """Return an OrderedDict with IDF headers as keys and list of all values per key"""

    with codecs.open(idf_file, encoding='utf-8') as fi:
        idf_raw = fi.readlines()
        idf_dict = OrderedDict()
        for row in idf_raw:
            idf_row = row.rstrip('\n').split('\t')
            # Skip empty lines and comments
            characters = ''.join(idf_row).strip()
            if not (len(characters) == 0 or characters.startswith('#')):
                row_label = idf_row.pop(0)
                if row_label in idf_dict:
                    idf_dict[row_label].extend(idf_row)
                    # The single cell pipeline does not handle duplicated IDF fields, need to collect here
                    idf_dict["duplicates"] = row_label
                else:
                    idf_dict[row_label] = idf_row
    return idf_dict


def read_sdrf_file(sdrf_file):
    """
    Read SDRF file and return the table content as nested list,
    the header row as list, and a dictionary of the fields and their indexes
    :param sdrf_file: string, path to SDRF file
    """

    with codecs.open(sdrf_file, encoding='utf-8') as sf:
        header = sf.readline().rstrip().split('\t')
        sdrf_raw = sf.readlines()

    sdrf_list = [x.rstrip('\n').split('\t') for x in sdrf_raw]
    header_dict = defaultdict(list)
    for i, field in enumerate(header):
        short_name = get_name(field)
        header_dict[short_name].append(i)

    return sdrf_list, header, header_dict


def get_name(header_string):
    """Return the first part of an SDRF header in lower case and without spaces."""
    no_spaces = header_string.replace(' ', '')
    field_name = no_spaces.split('[')[0]
    return field_name.lower()


def get_value(header_string):
    """Return the value within square brackets of an SDRF header."""
    field_value = header_string.split('[')[-1]
    return field_value.strip(']')


def get_sdrf_path(idf_file_path, logger, data_dir=None):
    """Read IDF and get the SDRF file name, look for the SDRF in the data directory (i.e. "unpacked")
    or in the same directory as the IDF.

    :param idf_file_path: full or relative path to IDF file
    :param logger: log handler
    :param data_dir: path to folder with SDRF
    :return: path to SDRF file
    """

    current_dir = os.path.dirname(idf_file_path)
    sdrf_file_path = ""
    if not data_dir:
        data_dir = DEFAULT_DATA_DIRECTORY
    # Figure out the name and location of sdrf files
    with codecs.open(idf_file_path, 'rU', encoding='utf-8') as f:
        # U flag makes it portable across in unix and windows (\n and \r\n are treated the same)
        for line in f:
            if re.search(SDRF_FILE_NAME_REGEX, line):
                sdrf_file_name = line.split('\t')[1].strip()
                data_path = os.path.join(current_dir, data_dir)
                if os.path.exists(data_path):
                    sdrf_file_path = os.path.join(data_path, sdrf_file_name)
                else:
                    sdrf_file_path = os.path.join(current_dir, sdrf_file_name)
    logger.debug("Generated SDRF file path: {}".format(sdrf_file_path))
    # Check if file exists and is in UTF-8 encoding and exit if not
    file_exists(sdrf_file_path)
    is_utf8(sdrf_file_path)

    return sdrf_file_path


def guess_submission_type_from_sdrf(sdrf_data, header, header_dict):
    """ Guess the basic experiment type (microarray or sequencing) from SDRF"""

    if 'arraydesignref' in header_dict or 'labeledextractname' in header_dict:
        return "microarray"
    elif "comment" in header_dict:
        for comment_index in header_dict.get("comment"):
            if get_value(header[comment_index]) in ("library construction", "single cell isolation"):
                return "singlecell"
    if "technologytype" in header_dict:
        index = header_dict.get("technologytype")
        if len(index) > 0:
            index = index[0]
            if sdrf_data[0][index] == "array assay":
                return "microarray"
            elif sdrf_data[0][index] == "sequencing assay":
                return "sequencing"


def guess_submission_type_from_idf(idf_dict):
    """Based on the experiment type, we can try to infer the basic experiment type
    This returns the type of the first experiment type found. We cannot handle mixed type experiments.
    """

    if "AEExperimentType" in idf_dict:
        all_types = get_controlled_vocabulary("ae_experiment_types", "atlas")

        for exptype in idf_dict["AEExperimentType"]:
            if exptype in all_types["sequencing"]:
                return "sequencing"
            elif exptype in all_types["microarray"]:
                return "microarray"
            elif exptype in all_types["singlecell"]:
                return "singlecell"