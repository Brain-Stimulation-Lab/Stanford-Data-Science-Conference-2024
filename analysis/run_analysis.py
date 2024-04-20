""" A script to scape databases like PubMed.

    Usage: python3 run_analysis.py -n 100
"""
import argparse
import pandas as pd
import re
import urllib.request
import ssl
from Bio import Entrez
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

PUBLICATION_YEAR = "2023"

def search(query, count):
    Entrez.email = "example@email.com"
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retmax=count,
                            retmode='xml',
                            year=PUBLICATION_YEAR,
                            term=query)
    
    results = Entrez.read(handle)
    return results

def fetch_article_text(pubmed_id):
    """
    Fetch article using BioC API.
    https://www.ncbi.nlm.nih.gov/research/bionlp/APIs/BioC-PMC/
    
    """
    print(pubmed_id)
    raw = urllib.request.urlopen(f"https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/{pubmed_id}/unicode")
    contents = raw.read()
    return contents.decode("utf-8")

def get_does_text_have_code(text):

    # Regex for 'code availible'
    has_code_available = re.search("code availab", text, re.IGNORECASE) is not None
    
    # Regex for github, bitbucket, or gitlab
    has_repository_link = re.search("bitbucket|github|gitlab", text, re.IGNORECASE) is not None

    return has_code_available, has_repository_link

def ssl_patch():
    try:
        _create_unverified_https_context = ssl._create_unverified_context
    except AttributeError:
        # Legacy Python that doesn't verify HTTPS certificates by default
        pass
    else:
        # Handle target environment that doesn't support HTTPS verification
        ssl._create_default_https_context = _create_unverified_https_context


def calculate_statistics(results, count):


    pubmed_ids = results['IdList']

    articles_data = {
        "pubmed_id": [],
        "fake": [],
        "repository_link": []
    }

    for pubmed_id in pubmed_ids:
        print("Searching details")

        contents = fetch_article_text(pubmed_id)
        result = get_does_text_have_code(contents)

        # articles_data[]
        # articles_data[]
        # articles_data[]


    df = pd.DataFrame.from_dict(articles_data)

    return df



def main():
    ssl_patch()

    parser = argparse.ArgumentParser(description='Search Pubmed and other databases.')
    parser.add_argument("-n", type=int, default=20, dest='count', help="The number of results to fetch (default = 20)")
    args = parser.parse_args()

    count = args.count

    query = "task fmri"

    print(f"Searching PubMed for '{query}'")
    results = search(query, count)
    print(f"Found {results['Count']} results searching for '{query}'")

    print(results)


    df = calculate_statistics(results, count)

    # df.to_csv("results.csv", index=False)


if __name__ == "__main__":
    main()