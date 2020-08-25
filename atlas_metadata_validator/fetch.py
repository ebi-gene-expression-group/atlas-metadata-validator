
import json
import logging
import pkg_resources
import re
import requests
import time
import urllib
import socket


# To store organisms that we have already looked-up in the taxonomy (this is slow...)
organism_lookup = {}

# To store the Atlas validation config with controlled vocabulary terms
atlas_config = None


def get_taxon(organism, logger=logging.getLogger()):
    """Return the NCBI taxonomy ID for a given species name."""

    if organism and organism not in organism_lookup:
        # If we have more than one organism mixed in one sample - in the case assign the 'mixed
        # sample' taxon_id (c.f. https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1427524)
        if re.search(r" and | \+ ", organism):
            return 1427524
        logger.info("Looking up species in NCBI taxonomy. Please wait...")
        db = 'taxonomy'
        term = organism.replace('(', ' ').replace(')', ' ')
        url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
        data = {'db': db, 'term': term, 'retmode': 'json'}
        r = requests.get(url, params=data)
        try:
            a = json.loads(r.text)
            taxon_id = int(a['esearchresult']['idlist'][0])
            organism_lookup[organism] = taxon_id
            return taxon_id
        except Exception as e:
            logger.error("Failed to retrieve organism data from ENA taxonomy service for {} due to {}".format(organism, str(e)))
    else:
        return organism_lookup.get(organism)


def is_valid_url(url, logger=None, retry=10):
    """Check if a given URL exists without downloading the page/file

    For HTTP and HTTPS URLs, urllib.requests returns a http.client.HTTPResponse object,
    for FTP URLs it returns a urllib.response.addinfourl object
    """

    # The global timeout for waiting for the response from the server before giving up
    timeout = 2
    socket.setdefaulttimeout(timeout)

    try:
        r = urllib.request.urlopen(url)
        logger.debug("Checking {}... Done.".format(url))
        if r:
            return True
    except urllib.error.URLError:
        if retry > 0:
            logger.debug("URI check failed for {}. Retrying {} more time(s).".format(url, str(retry)))
            time.sleep(60/retry)
            return is_valid_url(url, logger, retry-1)
        return False


def get_controlled_vocabulary(category, resource="atlas", logger=None):
    """Read the json with controlled vocab and return the dict for the given category.
    The config is fetched from the GitHub online copy as the primary source,
    if this isn't possible ("offline mode") the config packaged with the validator is used.
    """

    # Using global variable to only parse the file the first time it is used
    global atlas_config

    if resource == "atlas":
        if not atlas_config:
            resource_path = "atlas_validation_config.json"
            online_path = "https://raw.githubusercontent.com/ebi-gene-expression-group/metadata-validation-config/master/atlas_validation_config.json"

            try:
                if logger:
                    logger.debug("Getting online configuration from {}".format(online_path))

                raw_repsonse = requests.get(online_path)
                atlas_config = json.loads(raw_repsonse.text)

            except Exception:
                # Deliberately keeping this broad, any failure will lead to the back-up being used
                if logger:
                    logger.debug("Beware! Using local configuration file that might be out of date.")

                resource_package = "atlas_metadata_validator"
                atlas_config = json.loads(pkg_resources.resource_string(resource_package, resource_path))

        return atlas_config[category]

