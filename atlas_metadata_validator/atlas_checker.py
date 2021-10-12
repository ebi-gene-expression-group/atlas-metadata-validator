"""Module with validation checks to test for eligibility for Expression Atlas and Single Cell Expression Atlas"""

import logging
import re
from collections import defaultdict
from itertools import filterfalse

from atlas_metadata_validator.parser import simple_idf_parser, get_name, get_value, read_sdrf_file
from atlas_metadata_validator.fetch import get_taxon, check_urls, get_controlled_vocabulary


class AtlasMAGETABChecker:
    def __init__(self, idf_file, sdrf_file, submission_type, skip_file_checks=False, is_hca=False):
        self.idf_file = idf_file
        self.sdrf_file = sdrf_file
        self.submission_type = submission_type
        self.skip_file_checks = skip_file_checks
        self.is_hca = is_hca
        self.errors = set()

        # Read in IDF/SDRF files
        try:
            self.sdrf, self.sdrf_header, self.header_dict = read_sdrf_file(sdrf_file)
            self.idf = simple_idf_parser(idf_file)
        except Exception as e:
            raise Exception("Failed to open MAGE-TAB files: {}".format(e))

        # Get metadata fields
        self.sdrf_comment_values = [self.normalise_header(self.sdrf_header[i])
                                    for i in self.header_dict.get('comment', [])]
        self.sdrf_charactistics_values = [self.normalise_header(self.sdrf_header[i])
                                          for i in self.header_dict.get('characteristics', [])]
        self.sdrf_factor_values = [self.normalise_header(self.sdrf_header[i])
                                   for i in self.header_dict.get('factorvalue', [])]
        self.sdrf_values = self.sdrf_comment_values + self.sdrf_charactistics_values + self.sdrf_factor_values
        self.idf_values = [self.normalise_header(field) for field in self.idf]

        self.sample2datafile = self.map_samples_to_files()

    def run_general_checks(self, logger):
        """Checks applicable to all experiments"""

        # Warn about technical replicates
        if "technical replicate group" not in self.sdrf_comment_values:
            for s in self.sample2datafile:
                if s.endswith("PAIRED") and len(self.sample2datafile[s]) > 2:
                    logger.warn("Experiment is likely to contain technical replicates. "
                                "Please add \"Comment[technical replicate group]\".")
                    break
                elif not s.endswith("PAIRED") and len(self.sample2datafile[s]) > 1:
                        logger.warn("Experiment is likely to contain technical replicates. "
                                    "Please add \"Comment[technical replicate group]\".")
                        break

        # Required IDF fields
        required_idf_fields = get_controlled_vocabulary("required_idf_fields", "atlas", logger)
        for field in required_idf_fields:
            if field.lower() not in self.idf_values:
                logger.error("No \"{}\" found in IDF.".format(field))
                self.errors.add("GEN-E01")

        # No duplications between comments and characteristics
        duplicates = set(self.sdrf_comment_values).intersection(self.sdrf_charactistics_values)
        duplicates.update(set(self.sdrf_comment_values).intersection(self.sdrf_factor_values))
        for c in duplicates:
            logger.error("Column name \"{}\" appears as Comment and Characteristics/Factor Value.".format(c))
            self.errors.add("GEN-E02")

        # Species checks
        for i, field in enumerate(self.sdrf_header):
            if re.search("organism", get_value(field), re.IGNORECASE):
                organisms = {row[i] for row in self.sdrf}
                # Only 1 species allowed
                if len(organisms) > 1:
                    logger.error("Experiment contains more than 1 organism.")
                    self.errors.add("GEN-E03")
                # Organism must be found in NCBI taxonomy
                for o in organisms:
                    taxon_id = get_taxon(o)
                    if not taxon_id:
                        logger.error("Organism \"{}\" is not found in NCBI taxonomy.".format(o))
                        self.errors.add("GEN-E04")
                break

        # Sequencing experiments must have RUN or ENA_RUN comment column
        if self.submission_type in ("sequencing", "singlecell"):
            if not ("run" in self.sdrf_values or "ena_run" in self.sdrf_values):
                logger.error("No ENA_RUN or RUN column found in SDRF.")
                self.errors.add("GEN-E05")

        # FASTQ_URIs checks
        # We allow two types of values here: web URIs ("uris") or submitted file names ("non_uris")
        uri_pattern = re.compile(r"^(https?|ftp)://", re.IGNORECASE)  # Identify based on start of string
        non_uri_pattern = re.compile(r"^[^/\s]+$")  # Not allowing slashes and spaces if we don't have a web address

        uri_dict = defaultdict(list)
        uris_to_check = defaultdict(set)

        # Fetch the allowed field names for raw data URIs and get values from SDRF (store in uri_dict)
        uri_fields = get_controlled_vocabulary("raw_data_download_sdrf_fields")
        for uri_field in uri_fields:
            uri_index = [i for i, c in enumerate(self.sdrf_header)
                         if re.search(uri_field, c, flags=re.IGNORECASE)]
            for row in self.sdrf:
                for i in uri_index:
                    uri_dict[uri_field].append(row[i])

        # Check format of the strings and collect the web URIs for look up
        for uri_field, uri_list in uri_dict.items():
            for uri in uri_list:
                if re.match(uri_pattern, str(uri)):
                    uris_to_check[uri_field].add(uri)
                elif not re.match(non_uri_pattern, uri):
                    logger.error("{} {} is not in the right format.".format(uri_field.upper(), uri))
                    self.errors.add("GEN-E15")

        # Check for duplicated URI values
        all_uris = [uri for uri_field, uri_list in uri_dict.items() for uri in uri_list]
        duplicates = {uri for uri in all_uris if all_uris.count(uri) > 1}
        if duplicates:
            for uri in duplicates:
                logger.error(f"{uri} is duplicated.")
                self.errors.add("GEN-E16")

        # Looking up web URIs if they are valid web addresses
        if not self.skip_file_checks:
            logger.info("Checking for valid web URIs. This may take a while... (Skip this check with -x option)")
            invalid_uris = check_urls(logger, uris_to_check)
            for uri_field in invalid_uris:
                for uri in invalid_uris[uri_field]:
                    logger.error("{} {} is not valid.".format(uri_field.upper(), uri))
                    self.errors.add("GEN-E06")

    def run_singlecell_checks(self, logger):
        """Check requirements for loading an experiment into Single Cell Expression Atlas"""

        # Single cell IDF checks
        required_comments = get_controlled_vocabulary("required_singlecell_idf_comments", "atlas", logger)
        for comment in required_comments:
            if comment.lower() not in self.idf_values:
                logger.error("Comment \"{}\" not found in IDF. Required for Single Cell Atlas.".format(comment))
                self.errors.add("SC-E01")
        optional_comments = get_controlled_vocabulary("optional_singlecell_idf_comments", "atlas", logger)
        for comment in optional_comments:
            if comment.lower() not in self.idf_values:
                logger.warn("Comment \"{}\" not found in IDF. This is optional.".format(comment))
        if self.idf.get("duplicates"):
            logger.error("Detected duplicated IDF field \"{}\".".format(self.idf.get("duplicates")))
            self.errors.add("SC-E11")

        # Atlas IDF comment value checks
        for k, attribs in self.idf.items():
            if re.search("EAAdditionalAttributes", k, flags=re.IGNORECASE):
                for attrib in attribs:
                    if attrib and attrib.strip() not in self.sdrf_values:
                        logger.error("Additional attribute \"{}\" not found in SDRF.".format(attrib))
                        self.errors.add("SC-E02")
            elif re.search("EAExpectedClusters", k, flags=re.IGNORECASE):
                for number in attribs:
                    if number and not re.match(r"^\d+$", number):
                        logger.error("Expected clusters value \"{}\" is not numerical.".format(number))
                        self.errors.add("SC-E03")
            elif re.search("EAExperimentType", k, flags=re.IGNORECASE):
                for attrib in attribs:
                    if attrib and attrib.strip() not in get_controlled_vocabulary("singlecell_experiment_type", "atlas", logger):
                        logger.error("Unknown EAExperimentType: \"{}\"".format(attrib))
                        self.errors.add("SC-E04")

        # Required SDRF fields
        required_sdrf_names = get_controlled_vocabulary("required_singlecell_sdrf_fields", "atlas", logger)
        # Slightly different SDRF comments required for HCA experiments
        if self.is_hca:
            logger.debug("Found HCA imported experiment, checking presence of special fields.")
            required_sdrf_names = get_controlled_vocabulary("required_hca_sdrf_fields", "atlas", logger)
        for field in required_sdrf_names:
            if not (field.lower() in self.sdrf_values or get_name(field) in self.header_dict):
                logger.error("Required SDRF field \"{}\" not found.".format(field))
                self.errors.add("SC-E05")

        # Require at least one of those fields to be present
        required_download_fields = get_controlled_vocabulary("raw_data_download_sdrf_fields", "atlas", logger)
        found_download_field = False
        for field in required_download_fields:
            if field.lower() in self.sdrf_values or get_name(field) in self.header_dict:
                found_download_field = True
                break
        if not found_download_field:
            logger.error("Required data download SDRF field \"{}\" not found.".format(" or ".join(required_download_fields)))
            self.errors.add("SC-E12")

        # Valid SDRF values
        library_construction_terms = get_controlled_vocabulary("singlecell_library_construction", "atlas", logger)
        sc_protocol_values = {}
        for i, c in enumerate(self.sdrf_header):
            # Check for supported library construction terms
            if re.search(r"library construction", self.normalise_header(c), flags=re.IGNORECASE):
                sc_protocol_values = {row[i] for row in self.sdrf}
                logger.debug("Found library construction: \"{}\"".format(", ".join(sc_protocol_values)))
                if len(sc_protocol_values) > 1:
                    logger.warn("Experiment contains more than 1 single cell library construction protocol.")
                for protocol in sc_protocol_values:
                    if protocol.lower() not in library_construction_terms.get("all", []):
                        logger.error("Library construction protocol \"{}\" is not supported for Expression Atlas."
                                     .format(protocol))
                        self.errors.add("SC-E06")
            # Not all rows should be "not OK"
            elif re.search(r"single cell (well)? quality", self.normalise_header(c), flags=re.IGNORECASE):
                well_quality_values = {row[i] for row in self.sdrf}
                if len(well_quality_values) == 1 and "not OK" in well_quality_values:
                    logger.error("Single cell quality values are all \"not OK\".")
                    self.errors.add("SC-E07")
            # Technical replicate group checks
            elif re.search(r"technical replicate group", self.normalise_header(c), flags=re.IGNORECASE):
                tech_rep_values = {row[i] for row in self.sdrf}
                wrong_values = {v for v in tech_rep_values if not re.match(r"^[A-Za-z0-9]+$", v)}
                # Technical replicate group column should not be empty
                if len(wrong_values) == 1 and "" in wrong_values:
                    logger.error("Technical replicate group values are all empty.")
                    self.errors.add("SC-E09")
                # Technical replicate group values must only contain letters and numbers
                elif len(wrong_values) > 0:
                    logger.error(
                        "Technical replicate group values can only contain letters and numbers and cannot be empty: {}"
                        .format(", ".join([str(v) for v in wrong_values])))
                    self.errors.add("SC-E08")

        # SDRF terms required for droplet experiments
        for protocol in sc_protocol_values:
            if protocol.lower() in library_construction_terms.get("droplet"):
                droplet_terms = get_controlled_vocabulary("required_droplet_sdrf_fields", "atlas", logger)
                allowed_read_values = get_controlled_vocabulary("allowed_read_values", "atlas", logger)
                for dt in droplet_terms:
                    if dt not in self.sdrf_values:
                        logger.error("Required SDRF droplet field \"{}\" not found.".format(dt))
                        self.errors.add("SC-E10")
                    else:
                        # Check that the values match the allowed
                        droplet_term_values = self.get_sdrf_values_for_field(dt, unique_only=True)
                        if dt.endswith("read"):
                            non_match = droplet_term_values.difference(set(allowed_read_values))
                            for droplet_value in non_match:
                                logger.error(f"Read value \"{droplet_value}\" for \"{dt}\" is not allowed.")
                                self.errors.add("SC-E13")

                # check for numerical values in the optional comments
                droplet_numerical_terms = get_controlled_vocabulary("optional_droplet_numerical_fields", "atlas", logger)
                numerical_value_pattern = re.compile(r"^\d+$")
                for dnt in droplet_numerical_terms:
                    droplet_term_values = self.get_sdrf_values_for_field(dnt, unique_only=True)
                    non_match = filterfalse(numerical_value_pattern.match, droplet_term_values)
                    for droplet_value in non_match:
                        logger.error(f"Value \"{droplet_value}\" for \"{dnt}\" is not a numerical value.")
                        self.errors.add("SC-E14")

                break

    def check_all(self, logger=logging.getLogger()):
        """Trigger all applicable checks"""

        self.run_general_checks(logger)
        if self.submission_type == "singlecell":
            self.run_singlecell_checks(logger)

    @staticmethod
    def normalise_header(field_name):
        """Strip field names such as Comment and make everything lowercase without spaces"""
        return get_value(field_name).lower().strip()

    def map_samples_to_files(self):
        """Make a dictionary with the mapping of samples to their associated raw data files"""
        # Compare source name values against files
        sample2datafile = defaultdict(list)
        # Find source name and data file index
        sample_node = "sourcename"
        if self.submission_type == "microarray":
            data_node = "arraydatafile"
        else:
            data_node = "scanname"
        layout_comment = "library_layout"
        # The fist columns for each node (expecting only one raw data column)
        sample_index = next(iter(self.header_dict.get(sample_node, [])), None)
        data_index = next(iter(self.header_dict.get(data_node, [])), None)
        layout_index = next(iter([i for i, x in enumerate(self.sdrf_header)
                                  if re.search(layout_comment, x, flags=re.IGNORECASE)]), None)
        # Collect all data files for each sample
        if sample_index is not None and data_index is not None:
            for row in self.sdrf:
                if layout_index is not None:
                    name = str(row[sample_index])+"~"+str(row[layout_index])
                    sample2datafile[name].append(row[data_index])
                else:
                    sample2datafile[row[sample_index]].append(row[data_index])

        return sample2datafile

    def get_sdrf_values_for_field(self, field_name, unique_only=False):
        """Look up the values for a given SDRF field name (term inside square brackets for Comments/Characteristics)
        and return them as a list. If the unique_only option is set to True, return a set"""
        if unique_only:
            return set(self.get_sdrf_values_for_field(field_name))
        return [row[i]
                for row in self.sdrf
                for i, c in enumerate(self.sdrf_header)
                if self.normalise_header(field_name) == self.normalise_header(c)]
