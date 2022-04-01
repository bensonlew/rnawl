# -*- coding: utf-8 -*-

from bson.objectid import ObjectId
from api_base import ApiBase
import json
import re

class PlsdaMixomics(ApiBase):
    def __init__(self, bind_object):
        super(PlsdaMixomics, self).__init__(bind_object)

    def add_plsda_result(self, dir_path, main_id=None, group_file=None):
        self.s_group = {}
        with open(group_file) as fr:
            fr.readline()
            for line in fr:
                spline = line.strip('').split()
                self.s_group[spline[0]] = spline[1]
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        self.main_id = main_id

        site_path = dir_path.rstrip('/') + '/plsda_sites.xls'
        self.insert_plsda_meta_comp(site_path)

        rotation_path = dir_path.rstrip('/') + '/plsda_rotation.xls'
        self.insert_plsda_meta_vars(rotation_path)

        importance_path = dir_path.rstrip('/') + '/plsda_importance.xls'
        self.insert_plsda_meta_group(importance_path)

        ellipse_path =  dir_path.rstrip('/') + '/ellipse.xls'
        self.insert_plsda_meta_comp_group_ellipse(ellipse_path)

        importance_pre_path = dir_path.rstrip('/') + '/plsda_importancepre.xls'
        self.insert_explain(importance_pre_path)
        self.update_main()

    def insert_plsda_meta_comp(self,infile):
        insert_data = []
        with open(infile) as fr:
            fr.readline()
            for line in fr:
                spline = line.strip().split('\t')
                tmp = {
                    'name' : spline[0],
                    'group' : self.s_group[spline[0]],
                    'pc1' : float(spline[1]),
                    'pc2': float(spline[2]),
                    #'pc3' : float(spline[3]),
                    'type' :'scatter',
                    'plsda_id' : self.main_id
                }
                if len(spline) > 3:
                    tmp['pc3'] = float(spline[3])
                insert_data.append(tmp)
        self.db['plsda_meta_comp'].insert_many(insert_data)

    def insert_plsda_meta_comp_group_ellipse(self,infile):
        insert_data = []
        with open(infile) as fr:
            groups = fr.readline().strip().split('\t')[1:]
            for line in fr:
                spline = line.strip().split('\t')
                for  id,group in enumerate(groups,1):
                    spl = spline[id].split(',')
                    pc_x = 'pc'+ re.search('(\d+)\D+(\d+)',spline[0]).group(1)
                    pc_y = 'pc'+ re.search('(\d+)\D+(\d+)',spline[0]).group(2)
                    tmp = {
                        'name' : pc_x+'|'+pc_y,
                        'group' : group,
                        #'data': {'cx':float(spl[0]),'c11': float(spl[2]),'c12': float(spl[3]) ,'cy':float(spl[1]) ,'c21': float(spl[4]), 'c22':float(spl[5])},
                        'x': pc_x,
                        'y': pc_y,
                        #'data':   ','.join([spl[5],spl[4],spl[1],spl[0],spl[3],spl[2]]),
                        "data" :  spline[id], #",".join([spl[0], spl[2],spl[3],spl[1],spl[4],spl[5]]),  #cx,c11,c12,cy,c21,c22
                        'type' : 'ellipse',
                        'plsda_id' : self.main_id,
                        'method' : 'trend',
                        'method_type' : 'group'
                    }
                    insert_data.append(tmp)
        self.db['plsda_meta_comp'].insert_many(insert_data)


    def insert_plsda_meta_comp_ellipse(self,infile, main_id, all_sample_ellipse):
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        insert_data = []
        with open(infile) as fr, open(all_sample_ellipse) as fr2:
            fr.readline()
            fr2.readline()
            for id, file_o in enumerate([fr,fr2]):
                for line in file_o:
                    spline = line.strip().split('\t')
                    pc_x = 'pc' + re.search('(\d+)\D+(\d+)',spline[0]).group(1)
                    pc_y= 'pc' + re.search('(\d+)\D+(\d+)',spline[0]).group(2)
                    tmp = {
                        'name' : pc_x+'|'+pc_y,
                        'group' : spline[1],
                        #'data' : ','.join([spline[7],spline[6],spline[5],spline[2],spline[4],spline[3]]),
                        #'data': {'cx':float(spline[2]),'c11': float(spline[3]),'c12': float(spline[4]) ,'cy':float(spline[5]) ,'c21': float(spline[6]), 'c22':float(spline[7])},
                        'data' : ','.join([spline[2],spline[3],spline[4], spline[5],spline[6],spline[7]]) , # cx,c11,c12,cy,c21,c22
                        'x': pc_x,
                        'y': pc_y,
                        'type' : 'ellipse',
                        'plsda_id' : main_id,
                        'method' : 'ci'
                    }
                    if  id == 0:
                        tmp['method_type'] = 'group'
                    else:
                        tmp['method_type'] = 'sample'
                    insert_data.append(tmp)
        self.db['plsda_meta_comp'].insert_many(insert_data)

    def insert_plsda_meta_vars(self,infile):
        insert_data = []
        with open(infile) as fr:
            fr.readline()
            for line in fr:
                spline = line.strip().split('\t')
                tmp = {
                    'name' : spline[0],
                    'pc1' : float(spline[1]),
                    'pc2': float(spline[2]),
                    #'pc3': float(spline[3]),
                    'plsda_id' : self.main_id
                }
                if len(spline) > 3:
                    tmp['pc3'] = float(spline[3])
                insert_data.append(tmp)
        self.db['plsda_meta_vars'].insert_many(insert_data)

    def insert_plsda_meta_group(self,infile):
        insert_data = []
        with open(infile) as fr:
            fr.readline()
            for line in fr:
                spline = line.strip().split('\t')
                tmp = {
                    'group' : spline[0],
                    'pc1' : float(spline[1]),
                    'pc2': float(spline[2]),
                    #'pc3' : float(spline[3]),
                    'plsda_id' : self.main_id
                }
                if len(spline) > 3:
                    tmp['pc3'] = float(spline[3])
                insert_data.append(tmp)
        self.db['plsda_meta_group'].insert_many(insert_data)

    def insert_explain(self,infile):
        comp_data = {'plsda_id' : self.main_id, 'type': "explain"}
        with open(infile) as fr:
            lines = fr.readlines()
            for id,line in enumerate(lines[1:],1):
                nid = 'pc'+str(id)
                num = line.strip().split('\t')[1]
                comp_data[nid] = float(num)

        self.db['plsda_meta_comp'].insert_one(comp_data)

    def update_main(self):
        table_data ={
            "condition":{"type": "explain"},
            "column" :[
                {"field": "pc1", "type": "string", "sort": "false", "title": "comp1"},
                {"field": "pc2", "type": "string", "sort": "false", "title": "comp2"},
                #{"field": "pc3", "type": "string", "sort": "false", "title": "comp3"}
                 ]
            }

        table_data2 = {
            "condition":{"type":"scatter"},
            "column":[
                {"field": "name", "type": "string", "sort": "false", "title": "sample"},
                {"field": "pc1", "type": "string", "sort": "false", "title": "comp1"},
                {"field": "pc2", "type": "string", "sort": "false", "title": "comp2"},
               #{"field": "pc3", "type": "string", "sort": "false", "title": "comp3"}
            ]
            }

        table_data3 = {
            "condition":{},
            "column":[
                {"field": "group", "type": "string", "sort": "false", "title": "group"},
                {"field": "pc1", "type": "string", "sort": "false", "title": "comp1"},
                {"field": "pc2", "type": "string", "sort": "false", "title": "comp2"}
            ]
            }

        table_data4 = {
            "condition":{},
            "column":[
                {"field": "name", "type": "string", "sort": "false", "title": "var"},
                {"field": "pc1", "type": "string", "sort": "false", "title": "comp1"},
                {"field": "pc2", "type": "string", "sort": "false", "title": "comp2"},
                #{"field": "pc3", "type": "string", "sort": "false", "title": "comp3"}
                 ]
            }

        scatter_data = {"name":"var","data":["pc1","pc2","pc3"],"category":"group","condition":{"type":"scatter"}}
        ellipse_data = {"name":"name","data_option":["ci","trend"],"data":["data","data"],"condition":{"type":"ellipse"}}

        update_data  = {
            'main_id' : self.main_id,
            'pc_table_data' : json.dumps(table_data),
            'sample_table_data': json.dumps(table_data2),
            'group_table_data': json.dumps(table_data3),
            'varload_table_data': json.dumps(table_data4),
            'scatter_data': json.dumps(scatter_data),
            'ellipse_data' : json.dumps(ellipse_data)
        }

        self.db['plsda_meta'].update({"_id": self.main_id}, {"$set": update_data})






