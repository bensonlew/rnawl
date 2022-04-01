# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
# __last_modify__ = '2019/4'
from biocluster.api.database.base import Base, report_check
import os,re
import datetime
import types
from bson.son import SON
from bson.objectid import ObjectId
import json
import copy


class DiffKegg(Base):
    def __init__(self, bind_object):
        super(DiffKegg, self).__init__(bind_object)
        self._project_type = "bacgenome"


    @report_check
    def add_diff_kegg(self, name, task_id=None, samples=None, graph=None, params=None):
        """
        没有实际使用，只在交互分析里面做
        :param name:
        :param task_id:
        :param group_name: 分组名称，逗号分隔
        :param graph: 生成图片的对象存储路径
        :param params:
        :return:
        """
        task_id = task_id if task_id else self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        collection = self.db['diff_kegg']
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'desc': '',
            'created_ts': created_ts,
            'name': name if name else "DIFFKEGG_origin",
            'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
            'status': 'end',
            'samples': samples,
            'graph_dir': graph
        }
        main_id = collection.insert_one(insert_data).inserted_id
        return main_id

    @report_check
    def add_diff_kegg_detail(self, main_id,samples, mark_dir,pathway=None):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="52804701")

        diff_kegg_table = self.db['diff_kegg'].find_one({'_id':main_id})
        task_id = diff_kegg_table['task_id']
        r = self.db['anno_kegg'].find_one({"task_id" :task_id})
        #kegg_id = r['_id']
        kegg_id = r['main_id']
        kegg_detail = self.db['anno_kegg_detail']
        sample_list = samples.split(',')
        k_pat = re.compile('(K\d+)')

        if not os.path.isdir(mark_dir):
            self.bind_object.set_error('所指定的文件夹不存在，请检查: %s' , variables=(mark_dir), code="52804702")
        data_list = []
        files = os.listdir(mark_dir)
        for f in files:
            if f.endswith('.mark'):
                mark_path = os.path.join(mark_dir, f)
            else:
                continue
            with open(mark_path, 'rb') as f:
                line = f.readline()
                tmp_pathway = line.strip().split("\t")[0]
            with open(mark_path, 'rb') as f:
                lines = f.readlines()
                for line in lines:
                    line = line.strip().split("\t")
                    ori_data = [
                        ("diff_kegg_id", main_id),
                        ("pathway_id", line[0]),
                        ("shape", line[1]),
                        ("coords", line[4]),
                        ("title", line[5]),
                        ("href", line[7])
                    ]
                    mat = k_pat.findall(line[5])

                    if mat:
                        #k = mat.group(0)
                        has_one_least = False
                        for k in mat:
                            has_ = False
                            data = copy.deepcopy(ori_data)
                            self.bind_object.logger.info('match:'+ k)
                            for sample in sample_list:
                                self.bind_object.logger.info('specimen_id:'+sample+',ko:'+k+ ',kegg_id:'+str(kegg_id))
                                res = kegg_detail.find({'specimen_id':sample, 'ko':k, 'kegg_id':kegg_id})
                                if res >0:
                                    tmp_genes = [r['gene_id'] for r in res]
                                    if len(tmp_genes) >0:
                                        has_ = True
                                        data.append((sample , ';'.join(tmp_genes)))
                                    else:
                                        data.append((sample ,''))

                            if has_:
                                has_one_least = True
                                data.append(("query",k))
                                data = SON(data)
                                data_list.append(data)

                        if not has_one_least:
                            data = SON(ori_data)
                            data_list.append(data)

                    else:
                        #self.bind_object.logger.info(data)
                        data = SON(ori_data)
                        data_list.append(data)

        try:
            collection = self.db["diff_kegg_detail"]
            collection.insert_many(data_list)
            if not pathway:
                pathway = tmp_pathway
            graph_dir = os.path.join(self.bind_object.sheet.output, "pathway_img/{}.png".format(pathway))
            collection_main = self.db["diff_kegg"]
            collection_main.update_one({"_id": main_id}, {"$set": {"graph_dir": graph_dir, "main_id":main_id}}, upsert=True)
        except Exception,e:
            self.bind_object.set_error("diff_kegg_detail error: %s" , variables=( e), code="52804703")
        else:
            self.bind_object.logger.info("导入diff_kegg_detail成功")





    # @report_check
    # def add_diff_kegg_diff(self, main_id, path):
    #     if not isinstance(main_id, ObjectId):
    #         if isinstance(main_id, types.StringTypes):
    #             main_id = ObjectId(main_id)
    #         else:
    #             self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="52804704")
    #     if not os.path.isfile(path):
    #         self.bind_object.set_error('所指定的文件不存在，请检查: %s' , variables=( path), code="52804705")
    #     data_list = []
    #     with open(path, 'rb') as f:
    #         lines = f.readlines()
    #         head = lines[0].strip().split('\t')[1:]
    #         for i in head:
    #             if i[-3:] == "_sd" and not self.have_export_pic:
    #                 self.profile[i[:-3]] = {}  # 常见self.profile,以分组名称为key
    #         for line in lines[1:]:
    #             line = line.strip().split('\t')
    #             data = [
    #                 ("query", line[0]),
    #                 ("diff_kegg_id", main_id)
    #             ]
    #             for index,i in enumerate(head):
    #                 if line[index + 1] == "NaN":
    #                     data.append((i, 1))  # 防止pvalue qvalue 为NaN，影响作图
    #                 else:
    #                     data.append((i, float(line[index + 1])))
    #                 if i[-5:] != "value" and i[-3:] != "_sd":
    #                     self.profile[i][line[0]] = line[index + 1]
    #                     if float(line[index + 1]) > float(self.max_profile):
    #                         self.max_profile = line[index + 1]
    #             data = SON(data)
    #             data_list.append(data)
    #     try:
    #         collection = self.db["diff_kegg_diff"]
    #         collection.insert_many(data_list)
    #         if self.have_export_pic:
    #             self.update_kegg_detail(main_id)
    #     except Exception, e:
    #         self.bind_object.logger.error('导入%s信息出错：%s' % (path, e))
    #     else:
    #         self.bind_object.logger.info('导入%s信息成功！' % path)

    # def update_main(self, main_id, path):
    #     if not isinstance(main_id, ObjectId):
    #         if isinstance(main_id, types.StringTypes):
    #             main_id = ObjectId(main_id)
    #         else:
    #             self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="52804706")
    #     if not os.path.isfile(path):
    #         self.bind_object.set_error('所指定的文件不存在，请检查: %s' , variables=( path), code="52804707")
    #     graph_dir = os.path.join(self.bind_object.sheet.output, "pathway_img/")
    #     legend_color = [rotate(self.base_color, -i) for i in range(0, 81, 20)]  # 0% - 80%的颜色图例
    #     if float(self.max_profile) < 0.1:
    #         max_scale = float(re.match(r"0.0*\d?\d?",self.max_profile).group())
    #     else:
    #         max_scale = round(float(self.max_profile), 2)
    #     legend_label = [i * max_scale / 100 for i in range(0, 81, 20)]
    #     legend_update = zip(legend_label, legend_color)
    #     with open(path, "r") as f:
    #         line = f.readlines()[-1].strip().split("\t")
    #         href = line[-1]
    #         pathway_id = line[0]
    #     collection = self.db["diff_kegg"]
    #     collection.update_one({"_id": main_id}, {"$set": {"href": href, "graph_dir": graph_dir, "legend_color": legend_update}}, upsert=True)
    #     return pathway_id

    # def update_kegg_detail(self, main_id):
    #     collection = self.db["diff_kegg_detail"]
    #     for query_id in self.query:
    #         data = {}
    #         for group_name,group_profile in self.profile.items():
    #             profiles = group_profile.get(query_id, 0.0)
    #             if float(self.max_profile) < 0.1:
    #                 max_scale = float(re.match(r"0.0*\d?\d?",self.max_profile).group())
    #             else:
    #                 max_scale = round(float(self.max_profile), 2)
    #             profile_color = rotate(self.base_color, -float(profiles) / max_scale * 100)
    #             data[group_name] = profile_color
    #         collection.update_one({"_id": main_id, "query": query_id}, {"$set": data}, upsert=True)