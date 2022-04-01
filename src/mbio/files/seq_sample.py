# -*- coding: utf-8 -*-
# __author__ = 'guoquan'

import os
from biocluster.iofile import Directory
from biocluster.core.exceptions import FileError
import os
import re


class SeqSampleFile(Directory):
    def __init__(self):
        super(SeqSampleFile, self).__init__()
        self._config = None
        self._r1_file = None
        self._r2_file = None

    def set_path(self, path):
        self._config = path
        self.set_property("sample_name", path["sample_name"])
        self.set_property("library", path["library"])
        for key in ["genome_size", "library_name", "insert_size", "read_length"]:
            if key in path.keys():
                self.set_property(key, path[key])
        current_path = os.path.join(self.parent.path, self.name)
        super(SeqSampleFile, self).set_path(current_path)
        self.__set_files()

    def __set_files(self):
        if self.library.upper() in ["PE", "MP"]:
            if "sequence_r1_file_ids" not in self._config.keys() or "sequence_r2_file_ids" not in self._config.keys():
                raise FileError("样品%s 文库 %s 必须含有sequence_r1_file_ids和sequence_r2_file_ids" %
                                (self.name, self.library))
            # r1_file = os.path.join(self.parent.path, self.name, self._config["sequence_r1_file_ids"]["alias"])
            # r2_file = os.path.join(self.parent.path, self.name, self._config["sequence_r2_file_ids"]["alias"])
            self.add_files("sequence.fastq", self._config["sequence_r1_file_ids"]["alias"])
            self._r1_file = self._file_list[0]
            self.add_files("sequence.fastq", self._config["sequence_r2_file_ids"]["alias"])
            self._r2_file = self._file_list[1]
        elif self.library.lower() == "nanopore":
            if "sequence_file_ids" not in self._config.keys():
                raise FileError("样品%s 文库 %s 必须含有sequence_file_ids" % (self.name, self.library))
            for f in self._config["sequence_file_ids"]:
                # file_path = os.path.join(self.parent.path, self.name, f["alias"])
                self.add_files("sequence.fastq", f["alias"])
        elif self.library.lower() == "pacbio":
            if "sequence_file_ids" not in self._config.keys():
                raise FileError("样品%s 文库 %s 必须含有sequence_file_ids" % (self.name, self.library))
            for f in self._config["sequence_file_ids"]:
                # file_path = os.path.join(self.parent.path, self.name, f["alias"])
                if re.search('\.fastq$', f["alias"]) or re.search('\.fq$', f["alias"]) \
                    or re.search('\.fq\.gz$', f["alias"]) or re.search('\.fq\.gzip$', f["alias"]) \
                        or re.search('\.fastq\.gz$', f["alias"]) or re.search('\.fastq\.gzip$', f["alias"]):
                    self.add_files("sequence.fastq", f["alias"])
                elif re.search('\.bax\.h5$', f["alias"]) or re.search('\.bam$', f["alias"]) \
                        or re.search('\.bam\.pbi$', f["alias"]):
                    self.add_files("bacgenome.simple_file", f["alias"])

    @property
    def name(self):
        return self.prop["sample_name"]

    @property
    def genome_size(self):
        i = self.prop["genome_size"]
        if i:
            float(i)
        return i

    @property
    def library(self):
        return self.prop["library"]

    @property
    def library_name(self):
        return self.prop["library_name"]

    @property
    def insert_size(self):
        i = self.prop["insert_size"]
        if i:
            int(i)
        return i

    @property
    def read_length(self):
        i = self.prop["read_length"]
        if i:
            int(i)
        return i

    @property
    def r1_file(self):
        return self._r1_file

    @property
    def r2_file(self):
        return self._r2_file

