# -*- coding: utf-8 -*-
from biocluster.config import Config
from biocluster.api.database.base import Base
from mainapp.controllers.project.prok_rna_controller import ProkRNAController
from collections import defaultdict
from mbio.api.database.ref_rna_v2.api_base import ApiBase
from bson.objectid import ObjectId
import types


class GetRelationList(ProkRNAController):
    def __init__(self, project_type='prok_rna'):
        super(GetRelationList, self).__init__()
        self._project_type = project_type
        self.collections = self.db.collection_names()
        self.relation_list = list()

    def get_list(self):
        for col in self.collections:
            col_conn = self.db[col]
            if 'task_id' not in col_conn.find_one({}):
                continue
            detail = list()
            for l_col in self.collections:
                if col in l_col and col != l_col:
                    detail.append(l_col)
            if not detail:
                tur_res = (col, None, None)
                self.relation_list.append(tur_res)
            else:
                rela_id = ''
                if len(detail) == 1:
                    conn = self.db[detail[0]]
                    dict_one = conn.find_one({})
                    for key in dict_one:
                        if key != '_id' and key != 'main_id':
                            if isinstance(dict_one[key], ObjectId):
                                rela_id = key
                                break
                    tur_res = (col, detail[0], rela_id)
                    self.relation_list.append(tur_res)

                else:
                    rela_id = set()
                    for de in detail:
                        conn = self.db[de]
                        dict_one = conn.find_one({})
                        for key in dict_one:
                            if key != '_id' and key != 'main_id':
                                if isinstance(dict_one[key], ObjectId):
                                    rela_id.add(key)
                                    break
                    tur_res = (col, detail, ';'.join(rela_id))
                    self.relation_list.append(tur_res)

    def run(self):
        self.get_list()
        return self.relation_list


if __name__ == '__main__':
    # with open('relation.list', 'w') as rel_r:
    #     rel_r.write(GetRelationList().run())
    print(GetRelationList().run())