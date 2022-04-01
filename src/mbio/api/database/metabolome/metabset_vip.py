# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:20180523
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
# from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import glob
from mbio.packages.metabolome.common import check_metab
import pandas as pd


class MetabsetVip(Base):
    def __init__(self, bind_object):
        super(MetabsetVip, self).__init__(bind_object)
        self._project_type = "metabolome"

    @report_check
    def add_metabset_vip(self, metab_table_id, name=None, diff_id =None, main_id=None, params =None, diff_detail=None,table_type=None,metab_set=None):
        metab_table_id = self.check_id(metab_table_id, 'metab_table_id')
        if not main_id:
            task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else 'MetabVip_Origin',
                'params': params if params else '',
                'status': 'end',
                #'metab_table_id': metab_table_id,
                'main_id': '',
                'diff_detail': diff_detail if diff_detail else '',
                'diff_id': diff_id if diff_id else ''
            }
            if table_type:
                insert_data['table_type'] = table_type
            if metab_set:
                insert_data['metab_set'] = metab_set
            try:
                collection = self.db['metabset_vip']
                metabset_vip_id = collection.insert_one(insert_data).inserted_id
                self.update_table("main_id", metabset_vip_id, metabset_vip_id)
            except Exception, e:
                self.bind_object.set_error('导入metabset_vip主表异常:%s', variables=(e), code="54702101")
        else:
            self.update_table("main_id", main_id, main_id)
            metabset_vip_id = main_id
        return metabset_vip_id

    @report_check
    def add_metabset_vip_detail(self, metabset_vip_id, vip_dir, vip_type, metab_desc, scale=None,out_dir=None):
        metabset_vip_id = self.check_id(metabset_vip_id, "metabset_vip_id")
        if not os.path.exists(vip_dir):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(vip_dir), code="54702102")
        data_list = []
        result = self.db['metabset_vip'].find_one({'_id': metabset_vip_id})
        if not result:
            self.bind_object.set_error('找不到%s对应的主表id', variables=(metabset_vip_id), code="54702103")
        else:
            task_id = result['task_id']
            metab_table_id = eval(result['params'])["metab_table"]
            #samples_dic = name2id(task_id, type="task")
        groups = os.listdir(vip_dir)
        #mydict = self.get_some_info(metab_table_id, table_type)
        mydict = self.get_some_info_from_file(metab_desc)
        if not groups:
            self.bind_object.set_error('%s为空目录', variables=(vip_dir), code="54702104")
        for eachgroup in groups:
            if scale == "scale":
                vip_file = glob.glob(vip_dir + "/" + eachgroup + '/Vip_scale_exp.xls')
            else:
                vip_file = glob.glob(vip_dir + "/" + eachgroup + '/Vip_exp.xls')
            #tree_file = glob.glob(vip_dir + "/" + eachgroup + '/*tree*')
            ####vip_workdir = vip_dir.replace("/output", "")
            if out_dir:
                vip_workdir = out_dir
                tree_file = glob.glob(vip_workdir + "/" + eachgroup + '/metab_id.cluster_tree.xls')  #
            else:
                tree_file = glob.glob(vip_dir + "/" + eachgroup + '/metab_id.cluster_tree.xls')
            mdata = pd.read_table(vip_file[0],sep='\t',header=0)
            if 'ID' in mdata.columns:
                sams = mdata.columns[6:]
            else:
                sams = mdata.columns[5:]
            #with open(vip_file[0], 'rb') as f:
                # head = f.next()
                # heads = head.strip().split("\t")
                # if 'o_id' in heads:
                #     sams = heads[6:len(heads)]
                # else:
                #     sams = heads[5:len(heads)]
            for i in range(len(mdata)):
                #for line in f:
                #    line = line.strip().split('\t')
                #metab = mdata['Metabolite'][i].decode('utf-8')
                metab = mdata['Metabolite'][i]
                metab_search = metab
                if  metab_search in mydict:
                    metab_id = mydict[metab_search]
                else:
                    metab_id = ''
                if vip_type == "oplsda":
                    vip = ("%.4f" % float(mdata['Vip_oplsda'][i]))
                elif vip_type == "plsda":
                    vip = ("%.4f" % float(mdata['Vip_plsda'][i]))
                pvalue = float(mdata['P_value'][i])
                insert_data = {
                    'vip_id': metabset_vip_id,
                    'metab': metab,
                    'metab_id': metab_id,
                    'vip': float(vip),
                    'pvalue': pvalue,
                    'diff_group': eachgroup,
                    'type': "vip",
                }
                if 'ID' in mdata.columns:
                    insert_data['o_id'] = mdata['ID'][i]
                #for i in range(0, len(sams)):
                for s in sams:
                    sample_id = s
                    insert_data[sample_id] = float(("%.4f" % float(mdata[s][i])))
                data_list.append(insert_data)

            if tree_file:

                with open(tree_file[0], 'rb') as f:
                    metab_tree = f.readline().strip()
                    insert_data = {
                        'vip_id': metabset_vip_id,
                        'diff_group': eachgroup,
                        'type': "tree",
                        'metab_tree': metab_tree
                    }
                    data_list.append(insert_data)
        try:
            collection = self.db['metabset_vip_detail']
            if data_list:
                collection.insert_many(data_list)
            if 'ID' in mdata.columns:
                main_collection = self.db['metabset_vip']
                main_collection.update({"_id":metabset_vip_id},{"$set":{"var_head":'o_id'}})
        except Exception as e:
            self.bind_object.set_error("导入vip表格信息出错:%s", variables=(e), code="54702105")
        else:
            self.bind_object.logger.info("导入表格vip_detail信息成功!")

    @report_check
    def update_table(self, str, name, main_table_id):
        try:
            self.db['metabset_vip'].update_one({'_id': ObjectId(main_table_id)}, {'$set': {str: name}})
        except Exception as e:
            self.bind_object.set_error('更新metabset_vip主表%s字段出错:%s', variables=(str, e), code="54702106")

    @report_check
    def check_id(self, object_id, type):
        if not isinstance(object_id, ObjectId):
            if isinstance(object_id, types.StringTypes):
                object_id = ObjectId(object_id)
            else:
                self.bind_object.set_error('%s必须为ObjectId对象或其对应的字符串！', variables=(type), code="54702107")
        return object_id

    @report_check
    def get_some_info(self, metab_table_id, table_type):
        self.mydict = {}
        coll = "exp_" + table_type
        metab_table_id = self.check_id(metab_table_id, "metab_table_id")
        main_id = self.db["exp"].find_one({"_id": metab_table_id})["main_id"]
        main_id = self.check_id(main_id, "main_id")
        self.bind_object.logger.info(main_id)
        result = self.db[coll].find({'exp_id': main_id})
        self.bind_object.logger.info(result.count())
        for each in result:
            each_list = []
            metab_id = each['metab_id']
            metab = each['metab']
            self.mydict[metab] = metab_id
        return self.mydict

    def get_some_info_from_file(self,desc_file):
        self.mydict = {}
        with open(desc_file) as f:
            f.readline()
            for line in f:
                spline = line.split('\t')
                metab = spline[1]
                metab_id = spline[0]
                self.mydict[metab] = metab_id
        return self.mydict
