# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
import os
import json
from mbio.packages.ref_rna_v2.functions import tryforgood
from biocluster.config import Config
import time


class SearchPubmed(object):
    """
    This script is used to search pubmed database, use pmid or keyword.
    publication_fields = array('year', 'journal', 'title', 'authors', 'pmid', 'doi').
    """

    def __init__(self, pmid=None, query=None, start=0, stop=20):
        perl_path = Config().SOFTWARE_DIR + '/miniconda2/bin'
        os.environ['PATH'] = perl_path + ":" + os.environ['PATH']
        self.pmid = pmid
        self.query = query
        self.start = start
        self.stop = stop
        self.program = {
            'esearch': Config().SOFTWARE_DIR + '/bioinfo/biodb/edirect/esearch',
            'efetch': Config().SOFTWARE_DIR + '/bioinfo/biodb/edirect/efetch',
            'esummary': Config().SOFTWARE_DIR + '/bioinfo/biodb/edirect/esummary'
        }

    def run(self):
        if self.pmid:
            result = self.search_pmid(self.pmid)
            self.parse_result(result)
        if self.query:
            result = self.search_keyword(self.query, self.start, self.stop)
            self.parse_result(result)

    @tryforgood
    def search_pmid(self, pmid):
        cmd = "{} -db pubmed -query {} | {}  -format docsum -json".format(self.program['esearch'], pmid,
                                                                          self.program['efetch'])
        print(cmd)
        try:
            search_json = os.popen(cmd)
        except:
            print("search pubmed online failed: {}".format(search_json.readlines()))
            return False
        try:
            results = json.load(search_json)
            with open('pubmed_info.json', 'w') as f:
                f.write(json.dumps(results, indent=4))
            return results
        except:
            print("search pubmed online failed")
            return False

    @tryforgood
    def search_keyword(self, query, start, stop):
        cmd = "{} -db pubmed -query \"{}\" | {}  -format docsum -json -start {} -stop {}".format(self.program['esearch'],
                                                                                            query,
                                                                                            self.program['efetch'],
                                                                                            start, stop)
        print(cmd)
        try:
            search_json = os.popen(cmd)
        except:
            print("search pubmed online failed: {}".format(search_json.readlines()))
            time.sleep(15)
            return False
        try:
            results = json.load(search_json)
            with open('pubmed_info.json', 'w') as f:
                f.write(json.dumps(results, indent=4))
            return results
        except:
            print("search pubmed online failed")
            time.sleep(15)
            return False

    def parse_result(self, result):
        result_list = list()
        if isinstance(result['DocumentSummarySet']['DocumentSummary'], dict):
            summary = result['DocumentSummarySet']['DocumentSummary']
            authors = list()
            if isinstance(summary['Authors']['Author'], list):
                for each in summary['Authors']['Author']:
                    if isinstance(each, dict) and 'Name' in each:
                        authors.append(each['Name'])
                    else:
                        print(each)
            else:
                if 'Name' in summary['Authors']['Author']:
                    authors.append(summary['Authors']['Author']['Name'])
            data = {
                "pubdate": summary['PubDate'],
                "authors": ",".join(authors),
                "pmid": summary['Id'],
                "doi": summary['ELocationID'],
                "title": summary['Title'],
                "journal": summary['FullJournalName'],
            }
            result_list.append(data)
        elif isinstance(result['DocumentSummarySet']['DocumentSummary'], list):
            for summary in result['DocumentSummarySet']['DocumentSummary']:
                authors = list()
                if isinstance(summary['Authors']['Author'], list):
                    for each in summary['Authors']['Author']:
                        if isinstance(each, dict) and 'Name' in each:
                            authors.append(each['Name'])
                        else:
                            print(each)
                else:
                    if 'Name' in summary['Authors']['Author']:
                        authors.append(summary['Authors']['Author']['Name'])
                data = {
                    "pubdate": summary['PubDate'],
                    "authors": ",".join(authors),
                    "pmid": summary['Id'],
                    "title": summary['Title'],
                    "journal": summary['FullJournalName'],
                }
                result_list.append(data)
        with open('metadata.txt', 'w') as w:
            header = ["Pmid", "Title", "PubDate", "Journal", "Authors"]
            w.write("\t".join(header) + "\n")
            for doc in result_list:
                content = [doc["pmid"], doc["title"], doc["pubdate"], doc["journal"], doc["authors"]]
                w.write("\t".join(content) + "\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='This script is used to search pubmed database, use pmid or keyword.')
    parser.add_argument('-pmid', type=int, default=None, help='Pubmed id.')
    parser.add_argument('-keyword', type=str, default=None, help='Query string used to search pubmed.')
    parser.add_argument('-start', type=int, default=0, help='First record to fetch, used when -keyword is set.')
    parser.add_argument('-stop', type=int, default=20, help='Last record to fetch, used when -keyword is set.')

    args = parser.parse_args()

    if args.pmid and args.keyword:
        raise Exception("pmid and keyword cannot be searched at the same time.")
    if not (args.pmid or args.keyword):
        raise Exception("pmid and keyword cannot be empty at the same time.")
    search = SearchPubmed(args.pmid, args.keyword, args.start, args.stop)
    search.run()