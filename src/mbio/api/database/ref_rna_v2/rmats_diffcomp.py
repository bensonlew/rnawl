from mbio.api.database.ref_rna_v2.api_base import ApiBase
import os
import pandas as pd
from bson.objectid import ObjectId

class RmatsDiffcomp(ApiBase):
    def __init__(self, bind_object):
        super(RmatsDiffcomp, self).__init__(bind_object)
        self._project_type = 'ref_rna_v2'

    def add_rmats_diffcomp_detail(self,result_dir, main_id):
        result_file = os.path.join(result_dir, 'diffcomp_stat.txt')
        result_pd = pd.read_table(result_file)
        result_pd = result_pd.fillna('')
        result_pd['rmats_diffcomp_id'] = ObjectId(main_id)
        result_pd = result_pd.drop(columns=['Gene_name'])
        result_dict_list = result_pd.to_dict("records")
        not_float = ['New_ID', 'Type', 'Gene_ID']
        for dict in result_dict_list:
            for key in dict:
                if key not in not_float:
                    try:
                        dict[key] = float(dict[key])
                    except:
                        pass
        self.create_db_table("sg_splicing_rmats_diffcomp_detail",result_dict_list)