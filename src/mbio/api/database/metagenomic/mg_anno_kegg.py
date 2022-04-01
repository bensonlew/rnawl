# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:20171114
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
from mbio.packages.metagenomic.id_convert import name2id
from mbio.packages.metagenomic.id_convert import id2name
import json


class MgAnnoKegg(Base):
    def __init__(self, bind_object):
        super(MgAnnoKegg, self).__init__(bind_object)
        self._project_type = "metagenomic"
        self.sanger_prefix = Config().get_project_region_bucket(project_type="metagenomic")
        # self._db_name = Config().MONGODB + '_metagenomic'

    @report_check
    def add_anno_kegg(self, geneset_id, specimen, anno_file_path, main=True, main_table_id=None, name=None,
                      params=None, group_id=None, group_detail=None, software_ver={}):  # 删除xml_file
        if not isinstance(geneset_id, ObjectId):  # 检查传入的anno_kegg_id是否符合ObjectId类型
            if isinstance(geneset_id, types.StringTypes):  # 如果是string类型，则转化为ObjectId
                geneset_id_str = geneset_id
                geneset_id = ObjectId(geneset_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('geneset_id必须为ObjectId对象或其对应的字符串！', code="52802101")
        else:
            geneset_id_str = str(geneset_id)
            #if not os.path.exists(anno_file_path): # 调整为上传的永久路径，由于先导表，暂不检查
            #raise Exception('anno_file_path所指定的路径不存在，请检查！')
        if main:
            task_id = self.bind_object.sheet.id
            project_sn = self.bind_object.sheet.project_sn
            created_ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            # tmp = anno_file_path.split("/")
            # self.xml_path_prefix = "/".join(tmp[0:5])
            # if xml_file != None:
            #     if not os.path.exists(xml_file):
            #         self.bind_object.set_error('kegg_merge_xml所指定的路径%s不存在，请检查！', variables=(xml_file),
            #                                    code="52802102")
            #     else:
            #         # if self.bind_object.config.RGW_ENABLE:
            #         self.bind_object.logger.info("我用新的文件上传方式")
            #         xml_file = "/".join(xml_file.split('/')[8:])
            #         self.bind_object.upload_to_s3(xml_file, self.sanger_prefix + "metag/KEGG_Pathway/" + task_id + '/kegg_merge.xml')  # 改用相对路径
            #         # link_xml = self.xml_path_prefix + '/metag/KEGG_Pathway/' + task_id + "/kegg_merge.xml"  # 沿用旧存储格式
            #         link_xml = self.sanger_prefix + 'metag/KEGG_Pathway/' + task_id + '/kegg_merge.xml' # 换为相对路径
            #         # else:
            #         #     link_xml = self.xml_path_prefix + "/metag/KEGG_Pathway/" + task_id + "/kegg_merge.xml"
            #         #     if not os.path.exists(self.xml_path_prefix + "/metag/KEGG_Pathway/" + task_id):
            #         #         os.mkdir(self.xml_path_prefix + "/metag/KEGG_Pathway/" + task_id)
            #         #     if os.path.exists(link_xml):
            #         #         os.remove(link_xml)
            #         #     os.link(xml_file, link_xml)
            # else:
            #     link_xml = ""  # 所有if判断暂时去除，不在导表函数中上传结果
            if params == None:
                if group_id != None:
                    if not isinstance(group_id, types.StringTypes):
                        if isinstance(group_id, ObjectId):
                            group_id = str(group_id)
                        else:
                            self.bind_object.set_error('geneset_id必须为字符串类型！', code="52802103")
                    if group_detail == None:
                        self.bind_object.set_error('传入group_id时必须输入group_detail！', code="52802104")
                params = {
                    "database": "kegg",
                    "group_detail": group_detail,
                    "group_id": group_id,
                    "geneset_id": geneset_id_str,
                    "identity": 0,
                    "align_length": 0,
                    "submit_location": "annokegg",
                    "task_type": 2
                }
                params = json.dumps(params, sort_keys=True, separators=(",", ":"))
            insert_data = {
                'project_sn': project_sn,
                'task_id': task_id,
                'desc': '',
                'created_ts': created_ts,
                'name': name if name else "KEGG_Origin",
                'params': params,
                'status': 'end',
                'geneset_id': geneset_id,
                'specimen': specimen,
                'anno_file': anno_file_path,
                'lowest_level': 'Pathway',
                'version': "2",
                # "xml_file": link_xml
                # "xml_file": xml_file, # 取消了link_xml变量,
                "settled_params": json.dumps({"version": "kegg_v94.2"})
            }
            if software_ver:
                insert_data.update(software_ver)
            try:
                collection = self.db['anno_kegg']
                anno_kegg_id = collection.insert_one(insert_data).inserted_id
            except Exception, e:
                self.bind_object.logger.error('导入anno_kegg主表异常:{}'.format(e))
                self.bind_object.set_error("导入anno_kegg主表异常", code="52802105")
        else:
            if main_table_id is None:
                self.bind_object.set_error("main为False时需提供main_table_id!", code="52802106")
            if not isinstance(main_table_id, ObjectId):
                if isinstance(main_table_id, types.StringTypes):
                    main_table_id = ObjectId(main_table_id)
                else:
                    self.bind_object.set_error('main_table_id必须为ObjectId对象或其对应的字符串！', code="52802107")
            try:
                self.update_anno_file(main_table_id, anno_file_path)
                anno_kegg_id = main_table_id
            except Exception as e:
                self.bind_object.logger.error('更新anno_kegg主表anno_file_path出错:{}'.format(e))
                self.bind_object.set_error("更新anno_kegg主表anno_file_path出错", code="52802108")
        return anno_kegg_id

    @report_check
    def add_anno_kegg_gene(self, anno_kegg_id, kegg_profile_dir, update_main=True):
        kegg_profile = kegg_profile_dir + "/kegg_gene_profile.xls"
        if not isinstance(anno_kegg_id, ObjectId):  # 检查传入的anno_kegg_id是否符合ObjectId类型
            if isinstance(anno_kegg_id, types.StringTypes):  # 如果是string类型，则转化为ObjectId
                anno_kegg_id = ObjectId(anno_kegg_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('anno_kegg_id必须为ObjectId对象或其对应的字符串！', code="52802109")
        if not os.path.exists(kegg_profile):  # 检查要上传的数据表路径是否存在
            self.bind_object.set_error('kegg_gene_profile所指定的路径不存在，请检查！', code="52802110")
        data_list = []
        result = self.db['anno_kegg'].find_one({'_id': anno_kegg_id})
        if not result:
            self.bind_object.set_error('找不到anno_kegg_gene对应的主表id', code="52802111")
        else:
            task_id = result['task_id']
            # task_id = "metagenome3"
            samples_dic = name2id(task_id, type="task")
        with open(kegg_profile, 'rb') as f:
            head = f.next()
            if "#Gene" in head:
                heads = head.strip().split("\t")
                sams = heads[1:len(heads) - 1]
            else:
                self.bind_object.set_error('kegg_gene_profile.xls文件错误！', code="52802128")
            for line in f:
                line = line.strip().split('\t')
                gene = line[0]
                KO = line[-1]
                insert_data = {
                    'kegg_id': anno_kegg_id,
                    'gene': gene,
                    'orthology': KO
                }
                for i in range(0, len(sams)):
                    if not sams[i] == "Total":
                        sample_id = samples_dic[sams[i]]
                    else:
                        sample_id = sams[i]
                    insert_data[sample_id] = float(line[i + 1])
                data_list.append(insert_data)
            try:
                collection = self.db['anno_kegg_gene']
                # collection.insert_many(data_list)
                self.insert_detail(collection, data_list)
                main_table = self.db['anno_kegg']
                main_table.update_one({'_id': ObjectId(anno_kegg_id)}, {'$set': {"settled_params": json.dumps({"version": "kegg_v94.2"})}})
            except Exception as e:
                self.bind_object.logger.error("导入表格%s信息出错:%s" % (kegg_profile, e))
                self.bind_object.set_error("导入kegg_profile信息出错", code="52802113")
            else:
                self.bind_object.logger.info("导入表格%s信息成功!" % kegg_profile)
            if update_main:
                main_table = self.db['anno_kegg']
                if "Total" in sams:
                    sams.remove("Total")
                specimen = ",".join([samples_dic[i] for i in sams])
                # specimen = ",".join(sams)
                main_table.update_one({'_id': ObjectId(anno_kegg_id)}, {'$set': {'specimen': specimen,
                                                                                 "settled_params": json.dumps({"version": "kegg_v94.2"})}})

    @report_check
    def add_anno_kegg_orthology(self, anno_kegg_id, kegg_profile_dir):
        kegg_profile = kegg_profile_dir + "/kegg_KO_profile.xls"
        if not isinstance(anno_kegg_id, ObjectId):
            if isinstance(anno_kegg_id, types.StringTypes):
                anno_kegg_id = ObjectId(anno_kegg_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('anno_kegg_id必须为ObjectId对象或其对应的字符串！', code="52802109")
        if not os.path.exists(kegg_profile):
            self.bind_object.set_error('kegg_profile所指定的路径不存在，请检查！', code="52802114")
        data_list = []
        result = self.db['anno_kegg'].find_one({'_id': anno_kegg_id})
        if not result:
            self.bind_object.set_error('找不到anno_kegg_orthology对应的主表id', code="52802115")
        else:
            task_id = result['task_id']
            # task_id = "metagenome3"
            samples_dic = name2id(task_id, type="task")
        with open(kegg_profile, 'rb') as f:
            head = f.next()
            if "KO" in head:
                heads = head.strip().split("\t")
                sams = heads[1:len(heads) - 2]
            else:
                self.bind_object.set_error('kegg_KO_profile.xls文件错误！', code="52802116")
            for line in f:
                line = line.strip().split('\t')
                KO = line[0]
                des = line[-2]
                hyperlink = line[-1]
                insert_data = {
                    'kegg_id': anno_kegg_id,
                    'orthology': KO,
                    'description': des,
                    'hyperlink': hyperlink
                }
                for i in range(0, len(sams)):
                    if not sams[i] == "Total":
                        sample_id = samples_dic[sams[i]]
                    else:
                        sample_id = sams[i]
                    insert_data[sample_id] = float(line[i + 1])
                data_list.append(insert_data)
            try:
                collection = self.db['anno_kegg_orthology']
                # collection.insert_many(data_list)
                self.insert_detail(collection, data_list)
            except Exception as e:
                self.bind_object.logger.error("导入表格%s信息出错:%s" % (kegg_profile, e))
                self.bind_object.set_error("导入kegg_profile出错", code="52802117")
            else:
                self.bind_object.logger.info("导入表格%s信息成功!" % kegg_profile)

    @report_check
    def add_anno_kegg_module(self, anno_kegg_id, kegg_profile_dir):
        kegg_profile = kegg_profile_dir + "/kegg_module_profile.xls"
        if not isinstance(anno_kegg_id, ObjectId):
            if isinstance(anno_kegg_id, types.StringTypes):
                anno_kegg_id = ObjectId(anno_kegg_id)
            else:  # 如果是其他类型，则报错
                self.bind_object.set_error('anno_kegg_id必须为ObjectId对象或其对应的字符串！', code="52802109")
        if not os.path.exists(kegg_profile):
            self.bind_object.set_error('kegg_profile所指定的路径不存在，请检查！', code="52802114")
        data_list = []
        result = self.db['anno_kegg'].find_one({'_id': anno_kegg_id})
        if not result:
            self.bind_object.set_error('找不到kegg_module_profile对应的主表id', code="52802118")
        else:
            task_id = result['task_id']
            # task_id = "metagenome3"
            samples_dic = name2id(task_id, type="task")
        with open(kegg_profile, 'rb') as f:
            head = f.next()
            if "#Module" in head:
                heads = head.strip().split("\t")
                sams = heads[1:len(heads) - 1]
            else:
                self.bind_object.set_error('kegg_module_profile.xls文件错误！', code="52802119")
            for line in f:
                line = line.strip().split('\t')
                module = line[0]
                des = line[-1]
                insert_data = {
                    'kegg_id': anno_kegg_id,
                    'module': module,
                    'description': des
                }
                for i in range(0, len(sams)):
                    if not sams[i] == "Total":
                        sample_id = samples_dic[sams[i]]
                    else:
                        sample_id = sams[i]
                    insert_data[sample_id] = float(line[i + 1])
                data_list.append(insert_data)
            try:
                collection = self.db['anno_kegg_module']
                # collection.insert_many(data_list)
                self.insert_detail(collection, data_list)
            except Exception as e:
                self.bind_object.logger.error("导入表格%s信息出错:%s" % (kegg_profile, e))
                self.bind_object.set_error("导入kegg_profile信息出错", code="52802120")
            else:
                self.bind_object.logger.info("导入表格%s信息成功!" % kegg_profile)

    @report_check
    def add_anno_kegg_enzyme(self, anno_kegg_id, kegg_profile_dir):
        kegg_profile = kegg_profile_dir + "/kegg_enzyme_profile.xls"
        if not isinstance(anno_kegg_id, ObjectId):
            if isinstance(anno_kegg_id, types.StringTypes):
                anno_kegg_id = ObjectId(anno_kegg_id)
            else:
                self.bind_object.set_error('anno_kegg_id必须为ObjectId对象或其对应的字符串！', code="52802109")
        if not os.path.exists(kegg_profile):
            self.bind_object.set_error('kegg_profile所指定的路径不存在，请检查！', code="52802114")
        data_list = []
        result = self.db['anno_kegg'].find_one({'_id': anno_kegg_id})
        if not result:
            self.bind_object.set_error('找不到kegg_enzyme_profile对应的主表id', code="52802121")
        else:
            task_id = result['task_id']
            # task_id = "metagenome3"
            samples_dic = name2id(task_id, type="task")
        with open(kegg_profile, 'rb') as f:
            head = f.next()
            if "#Enzyme" in head:
                heads = head.strip().split("\t")
                sams = heads[1:len(heads) - 1]
            else:
                self.bind_object.set_error('kegg_enzyme_profile.xls文件错误！', code="52802122")
            for line in f:
                line = line.strip().split('\t')
                enzyme = line[0]
                des = line[-1].strip()
                insert_data = {
                    'kegg_id': anno_kegg_id,
                    'enzyme': enzyme,
                    'description': des
                }
                for i in range(0, len(sams)):
                    if not sams[i] == "Total":
                        sample_id = samples_dic[sams[i]]
                    else:
                        sample_id = sams[i]
                    insert_data[sample_id] = float(line[i + 1])
                data_list.append(insert_data)
            try:
                collection = self.db['anno_kegg_enzyme']
                # collection.insert_many(data_list)
                self.insert_detail(collection, data_list)
            except Exception as e:
                self.bind_object.logger.error("导入表格%s信息出错:%s" % (kegg_profile, e))
                self.bind_object.set_error("导入kegg_profile出错", code="52802123")
            else:
                self.bind_object.logger.info("导入表格%s信息成功!" % kegg_profile)

    @report_check
    def add_anno_kegg_pathway(self, anno_kegg_id, kegg_profile_dir, img_pathway):
        kegg_profile = kegg_profile_dir + "/kegg_pathway_profile.xls"
        if not isinstance(anno_kegg_id, ObjectId):
            if isinstance(anno_kegg_id, types.StringTypes):
                anno_kegg_id = ObjectId(anno_kegg_id)
            else:
                self.bind_object.set_error('anno_kegg_id必须为ObjectId对象或其对应的字符串！', code="52802109")
        if not os.path.exists(kegg_profile):
            self.bind_obect.set_error('kegg_profile所指定的路径不存在，请检查！', code="52802114")
        data_list = []
        result = self.db['anno_kegg'].find_one({'_id': anno_kegg_id})
        if not result:
            self.bind_object.set_error('找不到kegg_pathway_profile对应的主表id', code="52802124")
        else:
            task_id = result['task_id']
            # task_id = "metagenome3"
            samples_dic = name2id(task_id, type="task")
        self.bind_object.logger.info("start link pathway img")
        # pathway_img_dir = kegg_profile_dir + "/pathway_img"  # 取消了
        #if self.bind_object.config.RGW_ENABLE:
        # self.bind_object.logger.info("我用新的文件上传方式上传图片")
        # self.link_s3_img(pathway_img_dir, img_pathway)  # 不进行上传了
        #else:
        #self.link_img(pathway_img_dir, img_pathway)
        # self.bind_object.logger.info("链接pathway img成功")
        with open(kegg_profile, 'rb') as f:
            head = f.next()
            if "#Pathway" in head:
                heads = head.strip().split("\t")
                sams = heads[1:len(heads) - 2]
                self.bind_object.logger.info(sams)
            for line in f:
                line = line.strip().split('\t')
                pathway = line[0]
                self.add_anno_kegg_pic(anno_kegg_id, os.path.join(kegg_profile_dir, "pathway_img"), pathway)
                des = line[-2].strip()
                map = line[-1]
                # img_path = "/metag/KEGG_Pathway/" + str(task_id) + "/" + str(anno_kegg_id) + "/" + pathway
                img_path = os.path.join(img_pathway, pathway)  # img_pathway由将上传的路径，改为已上传的位置
                insert_data = {
                    'kegg_id': anno_kegg_id,
                    'pathway': pathway,
                    'description': des,
                    'pathwaymap': map,
                    'img_path': img_path
                }
                for i in range(0, len(sams)):
                    if not sams[i] == "Total":
                        sample_id = samples_dic[sams[i]]
                    else:
                        sample_id = sams[i]
                    insert_data[sample_id] = float(line[i + 1])
                data_list.append(insert_data)
            try:
                collection = self.db['anno_kegg_pathway']
                # collection.insert_many(data_list)
                self.insert_detail(collection, data_list)
            except Exception as e:
                self.bind_object.logger.error("导入表格%s信息出错:%s" % (kegg_profile, e))
                self.bind_object.set_error("导入kegg_profile出错", code="52802125")
            else:
                self.bind_object.logger.info("导入表格%s信息成功!" % kegg_profile)
        with open(kegg_profile_dir + "/kegg_level3_profile.xls", 'rb') as f3:
            head = f3.next()
            for line in f3:
                line = line.strip().split('\t')
                l1 = line[0]
                l2 = line[1]
                l3 = line[2].strip()
                try:
                    collection.update({"kegg_id": anno_kegg_id, 'description': l3},
                                      {'$set': {'level1': l1, "level2": l2}})
                except Exception as e:
                    self.bind_object.logger.error("更新%s信息出错level:%s" % (kegg_profile, e))
                    self.bind_object.set_error("更新kegg_profile信息出错", code="52802126")
            self.bind_object.logger.info("更新%s信息level成功!" % kegg_profile)
        with open(kegg_profile_dir + "/kegg_pathway_eachmap.xls", 'rb') as f4:
            head = f4.next()
            for line in f4:
                line = line.strip().split('\t')
                Pathway = line[0]
                ko_list = line[2]
                abundance = line[3]
                query_count = line[4]
                try:
                    collection.update({"kegg_id": anno_kegg_id, 'pathway': Pathway},
                                      {'$set': {'orthology_list': ko_list, 'orthology_abu': abundance,
                                                "query_count": query_count}})
                except Exception as e:
                    self.bind_object.logger.error("更新%s信息orthology_list、orthology_abu出错:%s" % (kegg_profile, e))
                    self.bind_object.set_error("更新orthology_list/orthology_abu出错", code="52802127")
            self.bind_object.logger.info("更新%s信息orthology_list、orthology_abu成功!" % kegg_profile)
        with open(kegg_profile_dir + "/kegg_level2_profile.xls", 'rb') as f5:
            head = f5.next()
            sams = head.strip().split("\t")[2:]
            data_list = []
            for line in f5:
                line = line.strip().split('\t')
                level1 = line[0]
                level2 = line[1]
                insert_data = {
                    'kegg_id': anno_kegg_id,
                    'level1': level1,
                    'level2': level2
                }
                for index,level2_prof in enumerate(line[2:]):
                    if sams[index] == "Total":
                        insert_data[sams[index]] = int(level2_prof)
                    else:
                        insert_data[samples_dic[sams[index]]] = int(level2_prof)
                data_list.append(insert_data)
            try:
                keggg_stat_collection = self.db['anno_kegg_stat']
                # keggg_stat_collection.insert_many(data_list)
                self.insert_detail(keggg_stat_collection, data_list)
            except Exception as e:
                self.bind_object.logger.error("import table kegg_level2_profile.xls error: %s" % e )
                self.bind_object.set_error("import anno_kegg_stat error", code="52802129")
            else:
                self.bind_object.logger.info("import mongo table kegg_level2_profile.xls success!")

    @report_check
    def add_anno_kegg_pic(self, anno_kegg_id, pic_dir, pathway):
        pic_path = os.path.join(pic_dir, pathway + ".html.mark")
        if not os.path.exists(pic_path):
            self.bind_object.logger.error("pic_path所指定的路径不存在，请检查！")
            return
        data_list = []
        with open(pic_path, "r") as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split("\t")
                data = [
                    ("kegg_id", anno_kegg_id),
                    ("pathway_id", line[0]),
                    ("shape", line[1]),
                    ("color", line[2]),
                    ("coords", line[4]),
                    ("title", line[5]),
                    ("query", line[6]),
                    ("href", line[7])
                ]
                data = SON(data)
                data_list.append(data)
        try:
            collection = self.db["anno_kegg_pic"]
            # collection.insert_many(data_list)
            self.insert_detail(collection, data_list)
        except Exception,e:
            self.bind_object.set_error("anno_kegg_pic error: %s" , variables=( e), code="52802130")
        else:
            self.bind_object.logger.info("导入anno_kegg_pic成功")

    @report_check
    def update_anno_file(self, main_table_id, anno_file):
        self.db['anno_kegg'].update_one({'_id': ObjectId(main_table_id)}, {'$set': {'anno_file': anno_file}})

    def link_s3_img(self, inputdir, outputdir):
        inputdir = "/".join(inputdir.split('/')[8:]) + "/"
        # outputdir = "/".join(outputdir.split('/')[5:]) + "/"
        # self.bind_object.upload_to_s3(inputdir, "s3://rerewrweset/" + outputdir)
        self.bind_object.upload_to_s3(inputdir, self.sanger_prefix + outputdir)

    @report_check
    def link_img(self, inputdir, outputdir):
        if not os.path.exists(outputdir):
            os.makedirs(outputdir)
        if os.path.exists(inputdir):
            files = os.listdir(inputdir)
            for i in files:
                eachfile = os.path.join(inputdir, i)
                newfile = os.path.join(outputdir, i)
                if os.path.exists(newfile):
                    os.remove(newfile)
                os.link(eachfile, newfile)

    def insert_detail(self, collection, datalist):
        if datalist:
            collection.insert_many(datalist)
        else:
            self.bind_object.logger.info("collection is empty")
