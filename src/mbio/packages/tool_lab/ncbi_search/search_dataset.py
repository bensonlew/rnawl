# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
import os
import json
from mbio.packages.ref_rna_v2.functions import tryforgood
from biocluster.config import Config
import requests
import re
import time


class SearchDataset(object):
    """
    This script is used to search dataset, use accession or keyword.
    dataset_fields = ["Accession", "TotalReads", "TotalBases", "FastqURL", "SraURL", "ScientificName", "Taxid",
                      "Platform", "LibraryName", "LibraryStrategy", "LibrarySource", "LibraryLayout",
                      "LibrarySelection", "StudyAcc", "StudyName", "Bioproject", "Biosample", "Submitter"]
    accession examples: SAMN00015175, SRX1067547, SRS084172, SRR1448793, SRA020176

    NCBI databases
    gds: GEO database, processed sequence data files
    sra: SRA database, raw sequence data files

    GDS accession
    Sample acc: GSM758572
    Series acc: GSE30567
    DataSets acc: GDS30567
    Platform acc: GPL10999

    SRA accession
    Bioproject: PRJNA128733 or PRJEB8073
    Study: SRP002628 or ERP002628
    Biosample: SAMN00015175 or SAMEA3178895
    Sample acc: SRS084172 or ERS631799
    Experiment acc: SRX022089 or ERX987758
    Submitter acc: SRA020176 or ERA4447
    Run acc: SRR1448793 or ERR908503 ---------target download data
    """

    def __init__(self, accession=None, query=None, start=0, stop=20):
        perl_path = Config().SOFTWARE_DIR + '/miniconda2/bin'
        os.environ['PATH'] = perl_path + ":" + os.environ['PATH']
        self.accession = accession
        self.query = query
        self.start = start
        self.stop = stop
        self.program = {
            'esearch': Config().SOFTWARE_DIR + '/bioinfo/biodb/edirect/esearch',
            'efetch': Config().SOFTWARE_DIR + '/bioinfo/biodb/edirect/efetch',
            'esummary': Config().SOFTWARE_DIR + '/bioinfo/biodb/edirect/esummary'
        }

    def run(self):
        if self.accession:
            if self.accession[0:3] in ["SRR", "SRX", "SRS", "ERR", "DRR", "SRP", "ERP", "DRP", "ERX", "ERS", "DRS",
                                       "DRX", "GSM", "SRA"] \
                    or self.accession.startswith("SAM") or self.accession.startswith("PRJEB") or \
                    self.accession.startswith("PRJNA"):
                result = self.search_sra(self.accession)
                if result:
                    self.parse_sra_result(result, self.accession)
                else:
                    print("search sra online failed")
            elif self.accession[0:3] in ["GSE"]:
                bioprojects = self.search_gds(self.accession)
                if isinstance(bioprojects, list):
                    for bioproject in bioprojects:
                        i += 1
                        result = self.search_sra(bioproject)
                        if result:
                            self.parse_sra_result(result, self.accession + "_" + str(i))
                        else:
                            print("search sra online failed")
                else:
                    result = self.search_sra(bioprojects)
                    if result:
                        self.parse_sra_result(result, self.accession)
                    else:
                        print("search sra online failed")
            else:
                raise Exception("{} not supported right now!".format(self.accession))
        if self.query:
            result = self.search_keyword(self.query, self.start, self.stop)
            if result:
                self.parse_sra_result(result)
            else:
                print("search sra online failed")

    @tryforgood
    def search_sra(self, accession):
        cmd = "{} -db sra -query {} | {}  -format docsum -json".format(self.program['esearch'], accession,
                                                                       self.program['efetch'])
        print(cmd)
        try:
            search_json = os.popen(cmd)
        except:
            print("search sra online failed: {}".format(search_json.readlines()))
            time.sleep(15)
            return False
        try:
            results = json.load(search_json)
            with open('sra_info.json', 'w') as f:
                f.write(json.dumps(results, indent=4))
            return results
        except:
            print("search sra online failed")
            time.sleep(15)
            return False

    @tryforgood
    def search_gds(self, accession):
        cmd = "{} -db gds -query {} | {}  -format docsum -json".format(self.program['esearch'], accession,
                                                                       self.program['efetch'])
        print(cmd)
        try:
            search_json = os.popen(cmd)
            result = json.load(search_json)
            if isinstance(result['DocumentSummarySet']['DocumentSummary'], dict):
                entry_type = result['DocumentSummarySet']['DocumentSummary']["entryType"]
                if entry_type == accession[0:3]:
                    bioproject = result['DocumentSummarySet']['DocumentSummary']["BioProject"]
                    if bioproject:
                        return bioproject
                    else:
                        print("No datasets found!")
                else:
                    print("No datasets found!")
            elif isinstance(result['DocumentSummarySet']['DocumentSummary'], list):
                bioprojects = list()
                for summary in result['DocumentSummarySet']['DocumentSummary']:
                    entry_type = summary["entryType"]
                    if entry_type == accession[0:3]:
                        bioproject = summary["BioProject"]
                        if bioproject and bioproject not in bioprojects:
                            bioprojects.append(bioproject)
                if bioprojects:
                    if len(bioprojects) == 1:
                        return "".join(bioprojects)
                    return bioprojects
                else:
                    print("No datasets found!")
        except:
            print("search gds online failed")
            time.sleep(15)
            return False

    @tryforgood
    def search_keyword(self, query, start, stop):
        cmd = "{} -db sra -query \"{}\" | {}  -format docsum -json -start {} -stop {}".format(self.program['esearch'],
                                                                                              query,
                                                                                              self.program['efetch'],
                                                                                              start, stop)
        print(cmd)
        try:
            search_json = os.popen(cmd)
        except:
            print("search sra online failed: {}".format(search_json.readlines()))
            time.sleep(15)
            return False
        try:
            results = json.load(search_json)
            with open('sra_info.json', 'w') as f:
                f.write(json.dumps(results, indent=4))
            return results
        except:
            print("search sra online failed")
            time.sleep(15)
            return False

    def parse_sra_result(self, result, accession=None):
        result_list = list()
        if isinstance(result['DocumentSummarySet']['DocumentSummary'], dict):
            expxml = result['DocumentSummarySet']['DocumentSummary']['ExpXml']
            runs = result['DocumentSummarySet']['DocumentSummary']['Runs']['Run']
            total_runs = expxml['Summary']['Statistics']['total_runs']
            if 'LIBRARY_STRATEGY' in expxml['Library_descriptor']:
                if isinstance(expxml['Library_descriptor']['LIBRARY_STRATEGY'], dict):
                    if 'content' in expxml['Library_descriptor']['LIBRARY_STRATEGY']:
                        library_strategy = expxml['Library_descriptor']['LIBRARY_STRATEGY']['content']
                    else:
                        library_strategy = ''
                else:
                    library_strategy = expxml['Library_descriptor']['LIBRARY_STRATEGY']
            else:
                library_strategy = ''
            if 'LIBRARY_SOURCE' in expxml['Library_descriptor']:
                if isinstance(expxml['Library_descriptor']['LIBRARY_SOURCE'], dict):
                    if 'content' in expxml['Library_descriptor']['LIBRARY_SOURCE']:
                        library_source = expxml['Library_descriptor']['LIBRARY_SOURCE']['content']
                    else:
                        library_source = ''
                else:
                    library_source = expxml['Library_descriptor']['LIBRARY_SOURCE']
            else:
                library_source = ''
            if 'LIBRARY_LAYOUT' in expxml['Library_descriptor']:
                if isinstance(expxml['Library_descriptor']['LIBRARY_LAYOUT'], dict):
                    if 'content' in expxml['Library_descriptor']['LIBRARY_LAYOUT']:
                        library_layout = expxml['Library_descriptor']['LIBRARY_LAYOUT']['content']
                    elif 'PAIRED' in expxml['Library_descriptor']['LIBRARY_LAYOUT']:
                        library_layout = 'PAIRED'
                    elif 'SINGLE' in expxml['Library_descriptor']['LIBRARY_LAYOUT']:
                        library_layout = 'SINGLE'
                    else:
                        library_layout = ''
                else:
                    library_layout = expxml['Library_descriptor']['LIBRARY_LAYOUT']
            else:
                library_layout = ''
            if 'LIBRARY_SELECTION' in expxml['Library_descriptor']:
                if isinstance(expxml['Library_descriptor']['LIBRARY_SELECTION'], dict):
                    if 'content' in expxml['Library_descriptor']['LIBRARY_SELECTION']:
                        library_selection = expxml['Library_descriptor']['LIBRARY_SELECTION']['content']
                    else:
                        library_selection = ''
                else:
                    library_selection = expxml['Library_descriptor']['LIBRARY_SELECTION']
            else:
                library_selection = ''
            if 'LIBRARY_NAME' in expxml['Library_descriptor']:
                if isinstance(expxml['Library_descriptor']['LIBRARY_NAME'], dict):
                    if 'content' in expxml['Library_descriptor']['LIBRARY_NAME']:
                        library_name = expxml['Library_descriptor']['LIBRARY_NAME']['content']
                    else:
                        library_name = ''
                else:
                    library_name = expxml['Library_descriptor']['LIBRARY_NAME']
            else:
                library_name = ''
            data = {
                "scientific_name": expxml['Organism']['ScientificName'] if 'ScientificName' in expxml[
                    'Organism'] else "",
                "taxonomy_id": expxml['Organism']['taxid'] if 'taxid' in expxml['Organism'] else "",
                "platform": expxml['Summary']['Platform']['instrument_model'] if 'instrument_model' in
                                                                                 expxml['Summary']['Platform'] else "",
                "study_acc": expxml['Study']['acc'] if 'acc' in expxml['Study'] else '',
                "study_name": expxml['Study']['name'] if 'name' in expxml['Study'] else '',
                "library_name": library_name,
                "library_strategy": library_strategy,
                "library_source": library_source,
                "library_layout": library_layout,
                "library_selection": library_selection,
                "biosample": expxml['Biosample'] if 'Biosample' in expxml else '',
                "bioproject": expxml['Bioproject'] if 'Bioproject' in expxml else '',
                "submitter": expxml['Submitter']['acc'] if 'Submitter' in expxml else '',
            }
            if total_runs and int(total_runs) >= 2:
                for eachrun in runs:
                    acc = eachrun['acc']
                    if self.accession:
                        if self.accession[0:3] in ["SRR", "DRR", "ERR"]:
                            if self.accession != acc:
                                continue
                    sra_ftp = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/" + acc[
                                                                                                    0:3] + "/" + acc[
                                                                                                                 0:6] + "/" + acc[
                                                                                                                              0:9] + "/" + acc + ".sra"
                    url = "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=" + acc + "&result=read_run&fields=run_accession,fastq_ftp"
                    try:
                        r = requests.get(url)
                        fastq_ftp = r.text.encode('utf-8').split("\n")[1].split("\t")[1]
                    except:
                        fastq_ftp = ""
                    data.update({
                        "acc": acc,
                        "total_bases": eachrun['total_bases'],
                        "total_reads": eachrun['total_spots'],
                        "fastq_ftp": fastq_ftp,
                        "sra_ftp": sra_ftp,
                    })
                    run_data = data.copy()
                    result_list.append(run_data)
            elif total_runs and int(total_runs) == 1:
                acc = runs['acc']
                if self.accession:
                    if self.accession[0:3] in ["SRR", "DRR", "ERR"]:
                        if self.accession != acc:
                            return
                sra_ftp = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/" + acc[
                                                                                                0:3] + "/" + acc[
                                                                                                             0:6] + "/" + acc[
                                                                                                                          0:9] + "/" + acc + ".sra"
                url = "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=" + acc + "&result=read_run&fields=run_accession,fastq_ftp"
                try:
                    r = requests.get(url)
                    fastq_ftp = r.text.encode('utf-8').split("\n")[1].split("\t")[1]
                except:
                    fastq_ftp = ""
                data.update({
                    "acc": acc,
                    "total_bases": runs['total_bases'],
                    "total_reads": runs['total_spots'],
                    "fastq_ftp": fastq_ftp,
                    "sra_ftp": sra_ftp,
                })
                run_data = data.copy()
                result_list.append(run_data)
            else:
                return
        elif isinstance(result['DocumentSummarySet']['DocumentSummary'], list):
            for summary in result['DocumentSummarySet']['DocumentSummary']:
                if 'ExpXml' in summary:
                    expxml = summary['ExpXml']
                    if 'LIBRARY_STRATEGY' in expxml['Library_descriptor']:
                        if isinstance(expxml['Library_descriptor']['LIBRARY_STRATEGY'], dict):
                            if 'content' in expxml['Library_descriptor']['LIBRARY_STRATEGY']:
                                library_strategy = expxml['Library_descriptor']['LIBRARY_STRATEGY']['content']
                            else:
                                library_strategy = ''
                        else:
                            library_strategy = expxml['Library_descriptor']['LIBRARY_STRATEGY']
                    else:
                        library_strategy = ''
                    if 'LIBRARY_SOURCE' in expxml['Library_descriptor']:
                        if isinstance(expxml['Library_descriptor']['LIBRARY_SOURCE'], dict):
                            if 'content' in expxml['Library_descriptor']['LIBRARY_SOURCE']:
                                library_source = expxml['Library_descriptor']['LIBRARY_SOURCE']['content']
                            else:
                                library_source = ''
                        else:
                            library_source = expxml['Library_descriptor']['LIBRARY_SOURCE']
                    else:
                        library_source = ''
                    if 'LIBRARY_LAYOUT' in expxml['Library_descriptor']:
                        if isinstance(expxml['Library_descriptor']['LIBRARY_LAYOUT'], dict):
                            if 'content' in expxml['Library_descriptor']['LIBRARY_LAYOUT']:
                                library_layout = expxml['Library_descriptor']['LIBRARY_LAYOUT']['content']
                            elif 'PAIRED' in expxml['Library_descriptor']['LIBRARY_LAYOUT']:
                                library_layout = 'PAIRED'
                            elif 'SINGLE' in expxml['Library_descriptor']['LIBRARY_LAYOUT']:
                                library_layout = 'SINGLE'
                            else:
                                library_layout = ''
                        else:
                            library_layout = expxml['Library_descriptor']['LIBRARY_LAYOUT']
                    else:
                        library_layout = ''
                    if 'LIBRARY_SELECTION' in expxml['Library_descriptor']:
                        if isinstance(expxml['Library_descriptor']['LIBRARY_SELECTION'], dict):
                            if 'content' in expxml['Library_descriptor']['LIBRARY_SELECTION']:
                                library_selection = expxml['Library_descriptor']['LIBRARY_SELECTION']['content']
                            else:
                                library_selection = ''
                        else:
                            library_selection = expxml['Library_descriptor']['LIBRARY_SELECTION']
                    else:
                        library_selection = ''
                    if 'LIBRARY_NAME' in expxml['Library_descriptor']:
                        if isinstance(expxml['Library_descriptor']['LIBRARY_NAME'], dict):
                            if 'content' in expxml['Library_descriptor']['LIBRARY_NAME']:
                                library_name = expxml['Library_descriptor']['LIBRARY_NAME']['content']
                            else:
                                library_name = ''
                        else:
                            library_name = expxml['Library_descriptor']['LIBRARY_NAME']
                    else:
                        library_name = ''
                    data = {
                        "scientific_name": expxml['Organism']['ScientificName'] if 'ScientificName' in expxml[
                            'Organism'] else '',
                        "taxonomy_id": expxml['Organism']['taxid'] if 'taxid' in expxml['Organism'] else '',
                        "platform": expxml['Summary']['Platform']['instrument_model'] if 'instrument_model' in
                                                                                         expxml['Summary'][
                                                                                             'Platform'] else '',
                        "study_acc": expxml['Study']['acc'] if 'acc' in expxml['Study'] else '',
                        "study_name": expxml['Study']['name'] if 'name' in expxml['Study'] else '',
                        "library_name": library_name,
                        "library_strategy": library_strategy,
                        "library_source": library_source,
                        "library_layout": library_layout,
                        "library_selection": library_selection,
                        "biosample": expxml['Biosample'],
                        "bioproject": expxml['Bioproject'],
                        "submitter": expxml['Submitter']['acc'],
                    }
                    total_runs = expxml['Summary']['Statistics']['total_runs']
                    if total_runs and int(total_runs) >= 2:
                        for eachrun in summary['Runs']['Run']:
                            acc = eachrun['acc']
                            if self.accession:
                                if self.accession[0:3] in ["SRR", "DRR", "ERR"]:
                                    if self.accession != acc:
                                        continue
                            sra_ftp = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/" + acc[
                                                                                                            0:3] + "/" + acc[
                                                                                                                         0:6] + "/" + acc[
                                                                                                                                      0:9] + "/" + acc + ".sra"
                            url = "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=" + acc + "&result=read_run&fields=run_accession,fastq_ftp"
                            try:
                                r = requests.get(url)
                                fastq_ftp = r.text.encode('utf-8').split("\n")[1].split("\t")[1]
                            except:
                                fastq_ftp = ""
                            data.update({
                                "acc": acc,
                                "total_bases": eachrun['total_bases'],
                                "total_reads": eachrun['total_spots'],
                                "fastq_ftp": fastq_ftp,
                                "sra_ftp": sra_ftp,
                            })
                            run_data = data.copy()
                            result_list.append(run_data)
                    elif total_runs and int(total_runs) == 1:
                        run = summary['Runs']['Run']
                        acc = run['acc']
                        if self.accession:
                            if self.accession[0:3] in ["SRR", "DRR", "ERR"]:
                                if self.accession != acc:
                                    continue
                        sra_ftp = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/" + acc[
                                                                                                        0:3] + "/" + acc[
                                                                                                                     0:6] + "/" + acc[
                                                                                                                                  0:9] + "/" + acc + ".sra"
                        url = "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=" + acc + "&result=read_run&fields=run_accession,fastq_ftp"
                        try:
                            r = requests.get(url)
                            fastq_ftp = r.text.encode('utf-8').split("\n")[1].split("\t")[1]
                        except:
                            fastq_ftp = ""
                        data.update({
                            "acc": acc,
                            "total_bases": run['total_bases'],
                            "total_reads": run['total_spots'],
                            "fastq_ftp": fastq_ftp,
                            "sra_ftp": sra_ftp,
                        })
                        run_data = data.copy()
                        result_list.append(run_data)
                    else:
                        continue
        if accession:
            output = '{}_metadata.xls'.format(accession)
        else:
            output = 'metadata.xls'
        with open(output, 'w') as w:
            header = ["Accession", "TotalReads", "TotalBases", "FastqURL", "SraURL", "ScientificName", "Taxid",
                      "Platform", "LibraryName", "LibraryStrategy", "LibrarySource", "LibraryLayout",
                      "LibrarySelection", "StudyAcc", "StudyName", "Bioproject", "Biosample", "Submitter"]
            w.write("\t".join(header) + "\n")
            for sample in result_list:
                content = [sample["acc"], sample["total_reads"], sample["total_bases"], sample["fastq_ftp"],
                           sample["sra_ftp"], sample["scientific_name"], sample["taxonomy_id"], sample["platform"],
                           sample["library_name"], sample["library_strategy"], sample["library_source"],
                           sample["library_layout"], sample["library_selection"], sample["study_acc"],
                           sample["study_name"], sample["bioproject"], sample["biosample"], sample["submitter"]]
                w.write("\t".join(content) + "\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='This script is used to search datasets, use accession or keyword.')
    parser.add_argument('-accession', type=str, default=None,
                        help='accession examples: SAMN00015175, SRX1067547, SRS084172, SRR1448793, SRA354008, GSE102990, PRJNA312176, PRJEB8073')
    parser.add_argument('-keyword', type=str, default=None, help='Query string used to search.')
    parser.add_argument('-start', type=int, default=0, help='First record to fetch, used when -keyword is set.')
    parser.add_argument('-stop', type=int, default=20, help='Last record to fetch, used when -keyword is set.')

    args = parser.parse_args()

    if args.accession and args.keyword:
        raise Exception("accession and keyword cannot be searched at the same time.")
    if not (args.accession or args.keyword):
        raise Exception("accession and keyword cannot be empty at the same time.")
    search = SearchDataset(args.accession, args.keyword, args.start, args.stop)
    search.run()
