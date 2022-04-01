# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# last_modify:20170926
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import re
# from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import json
from mbio.packages.metagenomic.id_convert import name2id


class BetaDiversity(Base):
    def __init__(self, bind_object):
        super(BetaDiversity, self).__init__(bind_object)
        self._project_type = "metagenomic"
        self.main_col = None
        self.detail_col = None
        self.main_id_name = "beta_diversity_id"

    # self._db_name = Config().MONGODB + '_metagenomic'

    @staticmethod
    def get_main_table_name(analysis_type):
        if analysis_type == 'pca':
            return 'PCA'
        elif analysis_type == 'pcoa':
            return 'PCoA'
        elif analysis_type == 'nmds':
            return 'NMDS'
        elif analysis_type == 'dbrda':
            return 'db-RDA'
        elif analysis_type == 'rda_cca':
            return 'RDA/CCA'
        else:
            self.bind_object.set_error('错误的分析类型', code="52800301")

    @report_check
    def add_beta_diversity(self, dir_path, analysis, web_path=None, main_id=None, main=None, env_id=None, group_id=None,
                           task_id=None, remove=None, params=None, anno_type=None, name=None, main_col="beta_diversity",
                           diff_check=""):
        self._tables = []  # 记录存入了哪些表格
        if task_id is None:
            task_id = self.bind_object.sheet.id
        if env_id and not isinstance(env_id, ObjectId):
            env_id = ObjectId(env_id)
        else:
            if 'env_id' in self.bind_object.sheet.data['options']:
                env_id = ObjectId(self.bind_object.option('env_id'))  # 仅仅即时计算直接绑定workflow对象
        _main_collection = self.db[main_col]
        self.main_col = main_col
        self.detail_col = main_col + "_detail"
        self.main_id_name = main_col + "_id"
        if diff_check:
            self.update_check(dir_path, diff_check, main_id)
        if main:
            insert_mongo_json = {
                'project_sn': self.bind_object.sheet.project_sn,
                'task_id': task_id,
                'name': name if name else BetaDiversity.get_main_table_name(analysis) + '_Origin',
                'table_type': analysis,
                'env_id': env_id,
                'params': json.dumps(params, sort_keys=True, separators=(',', ':')),
                'status': 'end',
                'desc': '',
                'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                'anno_type': anno_type
            }
            if analysis == 'rda_cca':  # 在主表中添加必要的rda或者是cca分类信息
                rda_files = os.listdir(dir_path.rstrip('/') + '/')
                if 'cca_sites.xls' in rda_files:
                    insert_mongo_json['rda_cca'] = 'cca'
                elif 'rda_sites.xls' in rda_files:
                    insert_mongo_json['rda_cca'] = 'rda'
                else:
                    self.bind_object.set_error('RDA/CCA分析没有生成正确的结果数据', code="52800302")
            try:
                main_id = _main_collection.insert_one(insert_mongo_json).inserted_id
            except:
                self.bind_object.set_error('RDA/CCA导表失败', code="52800303")
            else:
                self.bind_object.logger.info("RDA/CCA导表成功")
        else:
            if not main_id:
                self.bind_object.set_error('不写入主表时，需要提供主表ID', code="52800304")
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
            if analysis == 'rda_cca':  # 在主表中添加必要的rda或者是cca分类信息
                rda_files = os.listdir(dir_path.rstrip('/') + '/')
                if 'cca_sites.xls' in rda_files:
                    _main_collection.update_one({'_id': main_id}, {'$set': {'rda_cca': 'cca'}})
                elif 'rda_sites.xls' in rda_files:
                    _main_collection.update_one({'_id': main_id}, {'$set': {'rda_cca': 'rda'}})
                else:
                    self.bind_object.set_error('RDA/CCA分析没有生成正确的结果数据', code="52800302")
        result = _main_collection.find_one({'_id': main_id})
        if not result:
            self.bind_object.set_error('找不到主表id对应的表', code="52800305")
        if analysis == 'pca':
            site_path = dir_path.rstrip('/') + '/pca_sites.xls'
            rotation_path = dir_path.rstrip('/') + '/pca_rotation.xls'
            all_rotation_path = web_path.rstrip('/') + '/pca_rotation_all.xls'
            importance_path = dir_path.rstrip('/') + '/pca_importance.xls'
            _main_collection.update_one({'_id': main_id}, {'$set': {'download_file': all_rotation_path}})
            self.insert_table_detail(site_path, 'specimen', update_id=main_id)
            self.insert_table_detail(rotation_path, 'species', update_id=main_id, split_fullname=False)
            self.insert_table_detail(importance_path, 'importance', update_id=main_id, colls=['proportion_variance'])
            # by zhigang.zhao 20200907 分组椭圆需要 >>>
            if group_id not in ['all', 'All', 'ALL']:
                circle_path = dir_path.rstrip('/') + '/ellipse.xls'
                if os.path.exists(circle_path):
                    self.insert_group_ellipse(circle_path, main_id=main_id)
                    os.remove(circle_path)
            # <<<
            if result['env_id']:
                # filelist = os.listdir(dir_path.rstrip('/') + '/Pca')
                filelist = os.listdir(dir_path)
                if 'pca_envfit_factor_scores.xls' in filelist:
                    env_fac_path = dir_path.rstrip('/') + '/pca_envfit_factor_scores.xls'
                    env_fac_pr_path = dir_path.rstrip('/') + 'pca_envfit_factor.xls'
                    self.insert_table_detail(env_fac_path, 'factor', update_id=main_id)
                    self.insert_table_detail(env_fac_pr_path, 'factor_stat', update_id=main_id)
                if 'pca_envfit_vector_scores.xls' in filelist:
                    env_vec_path = dir_path.rstrip('/') + '/pca_envfit_vector_scores.xls'
                    env_vec_pr_path = dir_path.rstrip('/') + '/pca_envfit_vector.xls'
                    self.insert_table_detail(env_vec_path, 'vector', update_id=main_id)
                    self.insert_table_detail(env_vec_pr_path, 'vector_stat', update_id=main_id)
            else:
                pass
            self.bind_object.logger.info('beta_diversity:PCA分析结果导入数据库完成.')
        elif analysis == 'pcoa':
            site_path = dir_path.rstrip('/') + '/pcoa_sites.xls'
            self.insert_table_detail(site_path, 'specimen', update_id=main_id)
            importance_path = dir_path.rstrip('/') + '/pcoa_eigenvalues.xls'
            self.insert_table_detail(importance_path, 'eigenvalues', update_id=main_id)
            importance_pre_path = dir_path.rstrip('/') + '/pcoa_eigenvaluespre.xls'
            self.insert_table_detail(importance_pre_path, 'importance', update_id=main_id)
            # by zhigang.zhao 20200907 分组椭圆需要 >>>
            if group_id not in ['all', 'All', 'ALL']:
                circle_path = dir_path.rstrip('/') + '/ellipse.xls'
                if os.path.exists(circle_path):
                    self.insert_group_ellipse(circle_path, main_id=main_id)
                    os.remove(circle_path)
            # <<<
            self.bind_object.logger.info('beta_diversity:PCoA分析结果导入数据库完成.')
        elif analysis == 'nmds':
            site_path = dir_path.rstrip('/') + '/nmds_sites.xls'
            self.insert_table_detail(site_path, 'specimen', update_id=main_id,
                                     colls=['NMDS1', 'NMDS2'])  # 这里暂且只有两个轴，故只写两个，后面如果修改进行多维NMDS，结果应该有多列，此种方式应该进行修改
            # by zhigang.zhao 20200907 分组椭圆需要 >>>
            if group_id not in ['all', 'All', 'ALL']:
                circle_path = dir_path.rstrip('/') + '/ellipse.xls'
                if os.path.exists(circle_path):
                    self.insert_group_ellipse(circle_path, main_id=main_id,type="nmds")
                    os.remove(circle_path)
            # <<<
            nmds_stress = float(open(dir_path.rstrip('/') + '/nmds_stress.xls').readlines()[1])
            _main_collection.update_one({'_id': main_id}, {'$set': {'nmds_stress': nmds_stress}})
            self.bind_object.logger.info('beta_diversity:NMDS分析结果导入数据库完成.')
        elif analysis == 'dbrda':
            site_path = dir_path.rstrip('/') + '/db_rda_sites.xls'
            self.insert_table_detail(site_path, 'specimen', update_id=main_id)
            filelist = os.listdir(dir_path.rstrip('/') + '/')
            if 'db_rda_centroids.xls' in filelist:
                env_fac_path = dir_path.rstrip('/') + '/db_rda_centroids.xls'
                self.insert_table_detail(env_fac_path, 'factor', update_id=main_id)
            if 'db_rda_biplot.xls' in filelist:
                env_vec_path = dir_path.rstrip('/') + '/db_rda_biplot.xls'
                self.insert_table_detail(env_vec_path, 'vector', update_id=main_id)
            importance_path = dir_path.rstrip('/') + '/db_rda_importance.xls'
            self.insert_table_detail(importance_path, 'importance', update_id=main_id)
            if 'db_rda_plot_species_data.xls' in filelist:
                plot_species_path = dir_path.rstrip('/') + '/db_rda_plot_species_data.xls'
                self.insert_table_detail(plot_species_path, 'plot_species', update_id=main_id)
            if 'db_rda_envfit.xls' in filelist:
                envfit_path = dir_path.rstrip('/') + '/db_rda_envfit.xls'
                if len(open(envfit_path).readlines()) < 2:
                    os.remove(envfit_path)
                else:
                    self.insert_table_detail(envfit_path, 'envfit', update_id=main_id)
            self.bind_object.logger.info('beta_diversity:db_RDA分析结果导入数据库完成.')
        elif analysis == 'rda_cca':
            # if 'rda' in os.listdir(dir_path.rstrip('/') + '/Rda/')[1]:
            # if 'rda' in os.listdir(dir_path)[0]:  # modified @ 20181030
            #	rda_cca = 'rda'
            # else:
            #	rda_cca = 'cca'
            if "cca_sites.xls" in rda_files:
                rda_cca = "cca"
            else:
                rda_cca = "rda"
            site_path = dir_path.rstrip('/') + '/' + rda_cca + '_sites.xls'
            importance_path = dir_path.rstrip('/') + '/' + rda_cca + '_importance.xls'
            dca_path = dir_path.rstrip('/') + '/' + 'dca.xls'
            plot_species_path = dir_path.rstrip(
                '/') + '/' + rda_cca + '_plot_species_data.xls'  # add 4 lines by zhouxuan 20170123 20170401
            envfit_path = dir_path.rstrip('/') + '/' + rda_cca + '_envfit.xls'
            if os.path.exists(envfit_path):
                if len(open(envfit_path).readlines()) < 2:
                    os.remove(envfit_path)
                else:
                    self.insert_envfit_table(envfit_path, 'envfit', update_id=main_id)
            self.insert_table_detail(plot_species_path, 'plot_species', update_id=main_id)
            self.insert_table_detail(site_path, 'specimen', update_id=main_id)
            self.insert_table_detail(importance_path, 'importance', update_id=main_id, colls=['proportion_variance'])
            self.insert_table_detail(dca_path, 'dca', update_id=main_id)
            filelist = os.listdir(dir_path.rstrip('/') + '/')
            if (rda_cca + '_centroids.xls') in filelist:
                env_fac_path = dir_path.rstrip('/') + '/' + rda_cca + '_centroids.xls'
                self.insert_table_detail(env_fac_path, 'factor', update_id=main_id)
            if (rda_cca + '_biplot.xls') in filelist:
                env_vec_path = dir_path.rstrip('/') + '/' + rda_cca + '_biplot.xls'
                self.insert_table_detail(env_vec_path, 'vector', update_id=main_id, fileter_biplot=remove)
            self.bind_object.logger.info('beta_diversity:RDA/CCA分析结果导入数据库完成.')
        self.insert_main_tables(self._tables, update_id=main_id)
        return main_id

    def insert_table_detail(self, file_path, table_type, update_id,
                            coll_name='beta_diversity_detail',
                            main_coll='beta_diversity',
                            update_column=True, db=None, fileter_biplot=None, remove_key_blank=False,
                            split_fullname=False, colls=None):
        main_coll = self.main_col or main_coll
        coll_name = self.detail_col or coll_name
        self._tables.append(table_type)
        if not db:
            db = self.db
        collection = db[coll_name]
        main_collection = db[main_coll]
        beta_info = main_collection.find_one({'_id': update_id})
        task_name2id = beta_info["task_id"]
        self.sample_2_id = name2id(task_name2id, type="task")
        self.bind_object.logger.info(self.sample_2_id)
        with open(file_path, 'rb') as f:
            all_lines = f.readlines()
            columns = all_lines[0].rstrip().split('\t')[1:]
            if remove_key_blank:
                columns = [i.replace(' ', '') for i in columns]
            if colls:
                columns = colls
            data_temp = []
            for line in all_lines[1:]:
                values = line.rstrip().split('\t')
                if fileter_biplot:
                    if not isinstance(fileter_biplot, list):
                        self.bind_object.set_error('需要删除的匹配列必须为列表', code="52800306")
                    flag = 0
                    for i in fileter_biplot:
                        if re.match(r'^{}'.format(i), values[0]):
                            flag = 1
                    if flag:
                        continue
                else:
                    pass
                insert_data = {
                    self.main_id_name: update_id,
                    'type': table_type
                }
                if split_fullname:
                    insert_data['fullname'] = values[0]
                    if table_type == "specimen":
                        insert_data['name'] = self.sample_2_id[values[0].split(';')[-1].strip()]
                    else:
                        insert_data['name'] = values[0].split(';')[-1].strip()
                else:
                    if table_type == "specimen":
                        insert_data['name'] = self.sample_2_id[values[0]]
                    else:
                        insert_data['name'] = values[0]
                values_dict = dict(zip(columns, map(lambda x: round(float(x), 4), values[1:])))
                data_temp.append(dict(insert_data, **values_dict))
            if data_temp:
                try:
                    collection.insert_many(data_temp)
                except:
                    self.bind_object.set_error("导入detail表失败", code="52800307")
            else:
                return None
            if update_column:
                main_collection = db[main_coll]
                # default_column = {'specimen': 'detail_column', 'factor': 'factor_column',
                #                   'vector': 'vector_column',
                #                   'species': 'species_column', 'rotation': 'rotation_column'}
                default_column = {'specimen': 'specimen', 'factor': 'factor', 'vector': 'vector',
                                  'species': 'species', 'factor_stat': 'factor_stat',
                                  'vector_stat': 'vector_stat',
                                  'importance': 'importance', 'dca': 'dca', 'eigenvalues': 'eigenvalues',
                                  'plot_species': 'plot_species','envfit':'envfit'}
                if table_type in default_column:
                    main_collection.update_one({'_id': update_id},
                                               {'$set': {default_column[table_type]: ','.join(columns)}},
                                               upsert=False)
                else:
                    self.bind_object.set_error('错误的表格类型：%s不能在主表中插入相应表头', variables=(table_type), code="52800308")

    def insert_main_tables(self, tables, update_id, main='beta_diversity'):
        """
		"""
        main = self.main_col or main
        main_collection = self.db[main]
        main_collection.update_one({'_id': update_id},
                                   {'$set': {'tables': ','.join(tables)}},
                                   upsert=False)

    def insert_envfit_table(self, filepath, tabletype, update_id):
        """
		"""
        insert_data = []
        print filepath
        with open(filepath, 'rb') as r:
            # head = r.next().strip('\r\n')
            # head = re.split(' ', head)
            # new_head_old = head
            # new_head = []
            # for f in range(0, len(new_head_old)):
            #     pattern = re.compile('"(.*)"')
            #     new_head_ = pattern.findall(new_head_old[f])
            #     new_head.append(new_head_[0])
            # # m = re.match('/\"(.*)\"/', new_head[f])
            # # if m:
            # #	 new_head[f] = m.group(0)
            new_head = r.next().strip().split('\t')
            print(new_head)
            new_head[-1] = "p_values"
            new_head[-2] = "r2"
            self.new_head = new_head
            self.bind_object.logger.info(self.new_head)
            for line in r:
                line = line.rstrip("\r\n")
                line = re.split('\t', line)
                content = line[1:]
                # pattern = re.compile('"(.*)"')
                # env_name = pattern.findall(line[0])
                env_detail = dict()
                for i in range(1, len(content)+1):
                    env_detail[new_head[-i]] = content[-i]
                env_detail[self.main_id_name] = update_id
                env_detail['type'] = tabletype
                env_detail['name'] = line[0]
                insert_data.append(env_detail)
        try:
            collection = self.db[self.detail_col]
            collection.insert_many(insert_data)
        except Exception as e:
            self.bind_object.logger.error("导入beta_diversity_detail表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("导入beta_diversity_detail表格成功")
        self._tables.append(tabletype)
        # self.insert_main_tables(self._tables, update_id)
        main_collection = self.db['beta_diversity']
        main_collection.update_one({'_id': update_id},
                                   {'$set': {tabletype: ','.join(self.new_head)}},
                                   upsert=False)

    def insert_ellipse_table(self, infile, main_id):
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        insert_data = []
        with open(infile) as f:
            f.readline()
            name = ''
            tmp = {}
            for line in f:
                line = line.strip()
                spline = line.split('\t')
                if spline[0] != name:
                    insert_data.append(tmp)
                    name = spline[0]
                    # name_1 = spline[0].replace('PC','').replace('_','')
                    tmp = {'name': name}
                    tmp['type'] = 'circle'
                    tmp[self.main_id_name] = main_id
                k = spline[1]
                tmp[k] = ','.join(spline[2:])
            insert_data.append(tmp)

        try:
            detail_col = self.detail_col or 'beta_diversity_detail'
            collection = self.db[detail_col]
            collection.insert_many(insert_data[1:])
        except Exception as e:
            self.bind_object.logger.error("导入beta_diversity_detail表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("导入beta_diversity_detail表格成功")
        main_col = self.main_col or "beta_diversity"
        main_collection = self.db[main_col]
        main_collection.update_one({'_id': main_id}, {'$set': {'has_ci': 1}})

    def insert_group_ellipse(self, infile, main_id,type=None):
        if not isinstance(main_id, ObjectId):
            main_id = ObjectId(main_id)
        insert_data = []
        with open(infile) as f:
            a = f.readlines()
            name = ''
            tmp = {}
            for j in range(len(a)-1):
                if type:
                    name = "NMDS1_NMDS2"
                else:
                    if a[j + 1].split("\t")[0] == "12":
                        name = "PC1_PC2"
                    elif a[j + 1].split("\t")[0] == "13":
                        name = "PC1_PC3"
                    else:
                        name = "PC2_PC3"
                tmp = {'name': name}
                tmp['type'] = 'group_circle'
                tmp[self.main_id_name] = main_id
                for i in range(len(a[0].strip("\n").split("\t"))-1):
                    tmp[a[0].strip("\n").split("\t")[i+1]] = a[j+1].strip("\n").split("\t")[i+1]
                insert_data.append(tmp)
                tmp = {}
        try:
            detail_col = self.detail_col or "beta_diversity_detail"
            collection = self.db[detail_col]
            collection.insert_many(insert_data)
        except Exception as e:
            self.bind_object.logger.error("导入beta_diversity_detail表格信息出错:{}".format(e))
        else:
            self.bind_object.logger.info("导入beta_diversity_detail表格成功")
        main_col = self.main_col or "beta_diversity"
        main_collection = self.db[main_col]
        main_collection.update_one({'_id': main_id}, {'$set': {'has_group_circle': 1}})

    def update_check(self, dir_path, check_meth, main_id):
        if check_meth == "anosim":
            infos = {"plot": ["ANOSIM"],
                     "table": [['Method', 'Statistic R', 'P value', 'Permutation_num'],
                               ["ANOSIM"]]}
            file_path = os.path.join(dir_path, "format_results.xls")
            with open(file_path, 'r') as r:
                li = r.readlines()[1].strip().split('\t')
            infos["table"][1].extend(li[1:])
            infos["plot"].extend(li[1:3])
            mongo_data = {"ANOSIM": infos}
        elif check_meth == "adonis":
            infos = {"plot": ["Adonis"],
                     "table": [["Adonis", "Df", "Sums_of_sqs", "F.Model", "R2", "Pr(>F)"],
                               ["Group factors"],
                               ["Residuals"],
                               ["Total"]]}
            file_path = os.path.join(dir_path, "adonis_results.txt")
            lis = []
            with open(file_path, 'r') as r:
                r.readline()
                for l in r:
                    li = l.strip().split('\t')
                    lis.append([li[1], li[2], li[4], li[5], li[6]])
            for i in range(3):
                infos["table"][i + 1].extend(lis[i])
            infos["plot"].extend(lis[0][-2:])
            mongo_data = {"Adonis": infos}
        else:
            mongo_data = {"none": '-'}
        main_col = self.main_col or "beta_diversity"
        main_collection = self.db[main_col]
        main_collection.update_one({'_id': ObjectId(main_id)}, {'$set': {'diff_check': mongo_data}})
