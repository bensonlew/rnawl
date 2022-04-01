# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/3/10 18:07

from __future__ import print_function
import sys, re
import subprocess
from biocluster.iofile import File
from biocluster.core.exceptions import FileError
from biocluster.config import Config

''' https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/master/releases/so-xp.owl/so.obo'''
'''http://song.cvs.sourceforge.net/viewvc/song/ontology/so.obo?revision=1.263'''


class SequenceOntologyFile(File):
    def __init__(self):
        # self._path = os.path.join(os.getcwd(), 'temp_sequence_ontology.obo')
        # self._download(site)
        self._items_info = dict()
        self._candidate_items = set()
        self._properties = {}

    def check(self):
        super(SequenceOntologyFile, self).check()

    def _download(self, url):
        subprocess.call('wget -c -O {} {}'.format(self._path, url), shell=True)

    def parse(self):
        self._candidate_items = re.split(r'\[Term\]', open(self.path).read())
        self._candidate_items = self._candidate_items[1:len(self._candidate_items)]
        temp = re.split(r'\n+\[[^:\s]+?\]\n+', self._candidate_items[len(self._candidate_items) - 1])
        self._candidate_items.append(temp[0])

    def findAll(self, target_term, relation, term_id_set=None, term_name_set=None):
        if term_id_set == None:
            term_id_set = set()
        if term_name_set == None:
            term_name_set = set()
        for record in self._candidate_items:
            item_id_m = re.search(
                r'id:\s+(\S+)\nname:\s+(\S+)\n(.*?\n)*' + relation + ':\s+' + str(target_term) + '\s+!\s+(\S+)',
                record)
            if not item_id_m:
                continue
            else:
                son_so_id = item_id_m.group(1).strip()
                term_id_set.add(son_so_id)
                son_so_name = item_id_m.group(2).strip()
                term_name_set.add(son_so_name)
                print(son_so_id + ": " + son_so_name)
                term_id_set.union(self.findAll(son_so_id, term_id_set, term_name_set)[0])
                term_name_set.union(self.findAll(son_so_id, term_id_set, term_name_set)[1])
        return (term_id_set, term_name_set)


if __name__ == '__main__':
    # so_file_path = 'F:\\temp\\SO-Ontologies-master\\so.obo'
    so_file_path_shell = '/mnt/ilustre/users/sanger-dev/sg-users/jinlinfang/tmp/so.obo'
    so_file = SequenceOntologyFile()
    # so_file.set_path(so_file_path)
    so_file.set_path(so_file_path_shell)

    so_file.parse()
    (id_set, name_set) = so_file.findAll('SO:0000704')
    print(len(id_set))
    print(len(name_set))
    print(name_set)
