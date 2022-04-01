# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify: 20190124

import os
import re
import json
import types
import datetime
from bson.objectid import ObjectId
from biocluster.api.database.base import Base, report_check

class Datasplit(Base):
    def __init__(self, bind_object):
        super(Datasplit, self).__init__(bind_object)
        self._project_type = 'datasplit'

    def insert_sg_enzyme(self):
        """
        插入dna barcode信息
        """
        enzyme_dict = {"10": "GGCTAC", "16": "ACCTCT", "18": "ATAATC", "46": "TGCTTA", "52": "TTCCAC", "53": "CTCC",
                       "61": "AGGC", "62": "GATC", "63": "TCAC", "66": "TCACC", "78": "CATCT", "85": "TCGTT",
                       "96": "ACGTGTT", "100": "CGCTGAT", "101": "CGGTAGA", "104": "TAGCGGA", "108": "ACGACTAC",
                       "111": "TGCAAGGA", "114": "CCGGATAT", "119": "CCATGGGT", "121": "CGTGTGGT", "122": "GCTGTGGA",
                       "123": "GGATTGGT", "124": "GTGAGGGT", "M1": "AACG", "M3": "TTACA", "M6": "CCTCAG",
                       "M9": "GAAGATCC", "T1": "AGAACATA", "T5": "TTCAATG", "T7": "GAACTA", "T12": "CCTA",
                       "PstI124": "GTGAGGGT", "PstI123": "GGATTGGT", "PstI122": "GCTGTGGA", "PstI121": "CGTGTGGT",
                       "PstI119": "CCATGGGT", "PstI114": "CCGGATAT", "PstI111": "TGCAAGGA", "PstI108": "ACGACTAC",
                       "PstI104": "TAGCGGA", "PstI101": "CGGTAGA", "PstI100": "CGCTGAT", "PstI96": "ACGTGTT",
                       "PstI85": "TCGTT", "PstI78": "CATCT", "PstI62": "GATC", "PstI61": "AGGC", "PstI53": "CTCC",
                       "PstI52": "TTCCAC", "PstI46": "TGCTTA", "PstI18": "ATAATC", "PstI16": "ACCTCT", "PstI10": "GGCTAC",
                       "EcoRI124": "GTGAGGGT", "EcoRI123": "GGATTGGT", "EcoRI122": "GCTGTGGA", "EcoRI121": "CGTGTGGT",
                       "EcoRI119": "CCATGGGT", "EcoRI114": "CCGGATAT", "EcoRI111": "TGCAAGGA", "EcoRI108": "ACGACTAC",
                       "EcoRI104": "TAGCGGA", "EcoRI101": "CGGTAGA", "EcoRI100": "CGCTGAT", "EcoRI96": "ACGTGTT",
                       "EcoRI85": "TCGTT", "EcoRI78": "CATCT", "EcoRI62": "GATC", "EcoRI61": "AGGC", "EcoRI53": "CTCC",
                       "EcoRI52": "TTCCAC", "EcoRI46": "TGCTTA", "EcoRI18": "ATAATC", "EcoRI16": "ACCTCT", "EcoRI10": "GGCTAC",
                       "TaqI124": "GTGAGGGT", "TaqI123": "GGATTGGT", "TaqI122": "GCTGTGGA", "TaqI121": "CGTGTGGT",
                       "TaqI119": "CCATGGGT", "TaqI114": "CCGGATAT", "TaqI111": "TGCAAGGA", "TaqI108": "ACGACTAC",
                       "TaqI104": "TAGCGGA", "TaqI101": "CGGTAGA", "TaqI100": "CGCTGAT", "TaqI96": "ACGTGTT", "TaqI85": "TCGTT",
                       "TaqI78": "CATCT", "TaqI62": "GATC", "TaqI61": "AGGC", "TaqI53": "CTCC", "TaqI52": "TTCCAC", "TaqI46": "TGCTTA",
                       "TaqI18": "ATAATC", "TaqI16": "ACCTCT", "TaqI10": "GGCTAC", "MseI124": "GTGAGGGT", "MseI123": "GGATTGGT",
                       "MseI122": "GCTGTGGA", "MseI121": "CGTGTGGT", "MseI119": "CCATGGGT", "MseI114": "CCGGATAT", "MseI111": "TGCAAGGA",
                       "MseI108": "ACGACTAC", "MseI104": "TAGCGGA", "MseI101": "CGGTAGA", "MseI100": "CGCTGAT", "MseI96": "ACGTGTT",
                       "MseI85": "TCGTT", "MseI78": "CATCT", "MseI62": "GATC", "MseI61": "AGGC", "MseI53": "CTCC", "MseI52": "TTCCAC",
                       "MseI46": "TGCTTA", "MseI18": "ATAATC", "MseI16": "ACCTCT", "MseI10": "GGCTAC"}
        data_list = []
        for enzyme in enzyme_dict.keys():
            insert_data = {
                "enzyme": enzyme,
                "seq": enzyme_dict[enzyme]
            }
            data_list.append(insert_data)
        self.db["sg_enzyme"].insert_many(data_list)


if __name__ == "__main__":
    a = Datasplit(None)
    a.insert_sg_enzyme()
