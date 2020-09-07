
import json
import logging
import re
import requests
import os
import time
import urllib
import socket
import git


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
    The config is fetched from the GitHub online copy as the primary source, which is cloned
    on the first encounter and updated upon subsequent runs of the script.
    """

    # Using global variable to only parse the file the first time it is used
    global atlas_config

    if resource == "atlas":
        if not atlas_config:
            # Retrieve config file repo location from environment or use default
            config_repo = os.environ.get("VALIDATION_CONFIG_REPO") \
                          or "https://github.com/ebi-gene-expression-group/metadata-validation-config"
            config_file_name = os.environ.get("VALIDATION_CONFIG_FILE") or "atlas_validation_config.json"
            config_repo_name = os.path.basename(config_repo)

            local_config_path = os.environ.get("VALIDATION_CONFIG_LOCAL_PATH") or os.getcwd()
            local_repo = os.path.join(local_config_path, config_repo_name)
            local_config = os.path.join(local_repo, config_file_name)

            if not os.path.exists(local_repo):
                logger.debug("No local config found, cloning from remote repo {} to {}".format(config_repo, local_repo))
                git.Repo.clone_from(config_repo, local_repo, branch="master")

            try:
                logger.debug("Trying to update local config.")
                cloned_repo = git.Repo(local_repo)
                origin = cloned_repo.remotes.origin
                origin.pull()

            except Exception:
                logger.debug("Failed to update local config.")

            if os.path.exists(local_config):
                with open(local_config) as lc:
                    atlas_config = json.load(lc)

        return atlas_config[category]
