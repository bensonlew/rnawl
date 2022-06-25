# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
import os
import json
import unittest
import requests
from biocluster.workflow import Workflow
from mbio.packages.ref_rna_v2.functions import tryforgood


class GetSraDatasetWorkflow(Workflow):
    """
    Script to handle NCBI SRA lookups using different accessions
    SRR to NCBI, ERR to EBI and DRR to DDBJ

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
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GetSraDatasetWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "acc", "type": "string", "default": ""},  # SRP accession38
            {"name": "retmax", "type": "int", "default": 20},
            # Total number of UIDs from the retrieved set to be shown in the XML output (default=20).
            {"name": "retstart", "type": "int", "default": 0},
            # Sequential index of the first UID in the retrieved set to be shown in the XML output (default=0, corresponding
            # to the first record of the entire set). This parameter can be used in conjunction with retmax to download an
            # arbitrary subset of UIDs retrieved from a search.
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        perl_path = self.config.SOFTWARE_DIR + '/miniconda2/bin'
        os.environ['PATH'] = perl_path + ":" + os.environ['PATH']
        self.program = {
            'esearch': self.config.SOFTWARE_DIR + '/bioinfo/biodb/edirect/esearch',
            'efetch': self.config.SOFTWARE_DIR + '/bioinfo/biodb/edirect/efetch',
            'esummary': self.config.SOFTWARE_DIR + '/bioinfo/biodb/edirect/esummary'
        }

    def run(self):
        self.start_listener()
        self.fire("start")
        if self.option("acc")[0:3] in ["GSM", "SRR", "ERR", "DRR", "SRP", "ERP", "DRP"]:
            self.search_sra()
        self.end()

    @tryforgood
    def search_sample(self):
        results = list()
        cmd = "{} -db sra -query {} | {}  -format docsum -json".format(self.program['esearch'], self.option("acc"),
                                                                       self.program['efetch'])
        self.logger.info(cmd)
        try:
            search_json = os.popen(cmd)
        except:
            self.logger.info("search sra online failed: {}".format(search_json.readlines()))
            return False
        try:
            summarys = json.load(search_json)
            with open(self.work_dir + '/sra_info.json', 'w') as f:
                f.write(json.dumps(summarys, indent=4))
        except:
            self.logger.info("search sra online failed")
            return False
        if 'ExpXml' in summarys['DocumentSummarySet']['DocumentSummary']:
            dict = {
                "scientific_name": summarys['DocumentSummarySet']['DocumentSummary']['ExpXml']['Organism'][
                    'ScientificName'],
                "taxonomy_id": summarys['DocumentSummarySet']['DocumentSummary']['ExpXml']['Organism'][
                    'taxid'],
                "platform": summarys['DocumentSummarySet']['DocumentSummary']['ExpXml']['Summary'][
                    'Platform']['instrument_model'],
                "study_acc": summarys['DocumentSummarySet']['DocumentSummary']['ExpXml']['Study']['acc'],
                "study_name": summarys['DocumentSummarySet']['DocumentSummary']['ExpXml']['Study']['name'],
                "library_name": summarys['DocumentSummarySet']['DocumentSummary']['ExpXml']['Library_descriptor'][
                    'LIBRARY_NAME'],
                "library_strategy": summarys['DocumentSummarySet']['DocumentSummary']['ExpXml']['Library_descriptor'][
                    'LIBRARY_STRATEGY'],
                "library_source": summarys['DocumentSummarySet']['DocumentSummary']['ExpXml']['Library_descriptor'][
                    'LIBRARY_SOURCE'],
                "library_layout": summarys['DocumentSummarySet']['DocumentSummary']['ExpXml']['Library_descriptor'][
                    'LIBRARY_LAYOUT'].keys(),
                "library_selection": summarys['DocumentSummarySet']['DocumentSummary']['ExpXml']['Library_descriptor'][
                    'LIBRARY_SELECTION'],
            }
            total_runs = summarys['DocumentSummarySet']['DocumentSummary']['ExpXml']['Summary']['Statistics'][
                             'total_runs'],
            if total_runs >= 1:
                for run in summarys['DocumentSummarySet']['DocumentSummary']['Runs']:
                    acc = summarys['DocumentSummarySet']['DocumentSummary']['Runs'][run]['acc']
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
                    dict.update({
                        "acc": acc,
                        "total_bases": summarys['DocumentSummarySet']['DocumentSummary']['Runs'][run]['total_bases'],
                        "total_reads": summarys['DocumentSummarySet']['DocumentSummary']['Runs'][run]['total_spots'],
                        "fastq_ftp": fastq_ftp,
                        "sra_ftp": sra_ftp,
                    })
                    results.append(dict)
            else:
                results.append(dict)

        with open(self.work_dir + '/metadata.txt', 'w') as w:
            header = ["Accession", "TotalReads", "TotalBases", "FastqURL", "SraURL", "ScientificName", "Taxid",
                      "Platform", "LibraryName", "LibraryStrategy", "LibrarySource", "LibraryLayout",
                      "LibrarySelection", "StudyAcc", "StudyName"]
            w.write("\t".join(header) + "\n")
            for sample in results:
                content = [sample["acc"], sample["total_reads"], sample["total_bases"], sample["fastq_ftp"],
                           sample["sra_ftp"], sample["scientific_name"], sample["taxonomy_id"], sample["platform"],
                           sample["library_name"], sample["library_strategy"], sample["library_source"],
                           ",".join(sample["library_layout"]), sample["library_selection"], sample["study_acc"],
                           sample["study_name"]]
                w.write("\t".join(content) + "\n")

    @tryforgood
    def search_sra(self):
        results = list()
        cmd = "{} -db sra -query {} | {}  -format docsum -json".format(self.program['esearch'], self.option("acc"),
                                                                       self.program['efetch'])
        self.logger.info(cmd)
        try:
            search_json = os.popen(cmd)
        except:
            self.logger.info("search sra online failed: {}".format(search_json.readlines()))
            return False
        try:
            summarys = json.load(search_json)
            with open(self.work_dir + '/sra_info.json', 'w') as f:
                f.write(json.dumps(summarys, indent=4))
        except:
            self.logger.info("search sra online failed")
            return False
        for summary in summarys['DocumentSummarySet']['DocumentSummary']:
            if 'ExpXml' in summary:
                expxml = summary['ExpXml']
                dict = {
                    "scientific_name": expxml['Organism']['ScientificName'],
                    "taxonomy_id": expxml['Organism']['taxid'],
                    "platform": expxml['Summary']['Platform']['instrument_model'],
                    "study_acc": expxml['Study']['acc'],
                    "study_name": expxml['Study']['name'],
                    "library_name": expxml['Library_descriptor']['LIBRARY_NAME'],
                    "library_strategy": expxml['Library_descriptor']['LIBRARY_STRATEGY'],
                    "library_source": expxml['Library_descriptor']['LIBRARY_SOURCE'],
                    "library_layout": expxml['Library_descriptor']['LIBRARY_LAYOUT'].keys(),
                    "library_selection": expxml['Library_descriptor']['LIBRARY_SELECTION'],
                }
                total_runs = expxml['Summary']['Statistics']['total_runs'],
                if total_runs >= 1:
                    for eachrun in summary['Runs']:
                        run = summary['Runs'][eachrun]
                        acc = run['acc']
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
                        dict.update({
                            "acc": acc,
                            "total_bases": run['total_bases'],
                            "total_reads": run['total_spots'],
                            "fastq_ftp": fastq_ftp,
                            "sra_ftp": sra_ftp,
                        })
                        results.append(dict)
                else:
                    results.append(dict)
        with open(self.work_dir + '/metadata.txt', 'w') as w:
            header = ["Accession", "TotalReads", "TotalBases", "FastqURL", "SraURL", "ScientificName", "Taxid",
                      "Platform", "LibraryName", "LibraryStrategy", "LibrarySource", "LibraryLayout",
                      "LibrarySelection", "StudyAcc", "StudyName"]
            w.write("\t".join(header) + "\n")
            for sample in results:
                content = [sample["acc"], sample["total_reads"], sample["total_bases"], sample["fastq_ftp"],
                           sample["sra_ftp"], sample["scientific_name"], sample["taxonomy_id"], sample["platform"],
                           sample["library_name"], sample["library_strategy"], sample["library_source"],
                           ",".join(sample["library_layout"]), sample["library_selection"], sample["study_acc"],
                           sample["study_name"]]
                w.write("\t".join(content) + "\n")


    def end(self):
        super(GetSraDatasetWorkflow, self).end()

class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from biocluster.wpm.client import worker_client
        import datetime
        worker = worker_client()
        data = {
            "id": "SraExplorer_" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.get_sra_dataset",
            "options": dict(
                acc="SRP002628"
            )
        }
        info = worker.add_task(data)
        print info


if __name__ == '__main__':
    unittest.main()