
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


def get_controlled_vocabulary(category, resource="translations"):
    """Read the json with controlled vocab and return the dict for the given category.
    The resource parameter specifies which file to read."""
    resource_package = "atlas_metadata_validator"

    if resource == "ontology":
        resource_path = "ontology_terms.json"
    elif resource == "atlas":
        resource_path = "atlas_validation_config.json"
    else:
        resource_path = "term_translations.json"
    all_terms = json.loads(pkg_resources.resource_string(resource_package, resource_path))

    return all_terms[category]

