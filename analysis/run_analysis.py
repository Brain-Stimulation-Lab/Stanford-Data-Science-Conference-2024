""" A script to scape databases like PubMed.

    Usage: python3 run_analysis.py -n 100
"""
import argparse
import pandas as pd
import re
import urllib.request
import ssl
from Bio import Entrez
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

PUBLICATION_YEAR = "2023"

def search(query):
    Entrez.email = "example@email.com"
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            
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
    raw = urllib.request.urlopen(f"https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/{pubmed_id}/unicode")
    contents = raw.read()
    return contents.decode("utf-8")

def parse_article_text(text):

    # An arbitrary threshold of 20 characters chosen, BioC often returns articles that don't have full text
    has_full_text_available = len(text) > 20

    # Regex for 'code available', the standard
    has_code_available = re.search("code availab", text, re.IGNORECASE) is not None
    
    # Regex for github, bitbucket, or gitlab
    repository_match = re.search("bitbucket|github|gitlab", text, re.IGNORECASE)

    repository_link = repository_match.group(0) if repository_match is not None else False

    return has_full_text_available, has_code_available, repository_link

def ssl_patch():
    try:
        _create_unverified_https_context = ssl._create_unverified_context
    except AttributeError:
        # Legacy Python that doesn't verify HTTPS certificates by default
        pass
    else:
        # Handle target environment that doesn't support HTTPS verification
        ssl._create_default_https_context = _create_unverified_https_context


def calculate_statistics(results, max_results):

    pubmed_ids = results['IdList']

    articles_data = []

    print(len(pubmed_ids))

    for pubmed_id in pubmed_ids:
        print(f"Searching details for {pubmed_id}")

        contents = fetch_article_text(pubmed_id)

        # TODO filter original research

        has_full_text_available, has_code_available, repository_link = parse_article_text(contents)

        if has_full_text_available:
            articles_data += [[pubmed_id, has_full_text_available, has_code_available, repository_link]]
        
    df = pd.DataFrame(articles_data, columns=["PubMedID", "FullText", "CodeAvailable", "RepositoryLink"])

    return df



def main():
    ssl_patch()

    parser = argparse.ArgumentParser(description='Search Pubmed and other databases.')
    parser.add_argument("-n", type=int, default=20, dest='max_results', help="The number of results to fetch (default = 20)")
    args = parser.parse_args()

    max_results = args.max_results

    query = "task fmri"

    print(f"Searching PubMed for '{query}'")
    results = search(query)
    print(f"Found {results['Count']} results searching for '{query}'")

    print("Calculating statistics")
    df = calculate_statistics(results, max_results)

    print("DataFrame created")
    print(df)
    df.to_csv("results.csv", index=False)


if __name__ == "__main__":
    main()