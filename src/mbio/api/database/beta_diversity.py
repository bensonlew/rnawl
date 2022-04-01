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
            self.bind_object.set_error('错误的分析类型', code="51000201")

    @report_check
    def add_beta_diversity(self, dir_path, analysis, web_path=None,main_id=None, main=False, env_id=None, group_id=None,
                           task_id=None, remove=None, params=None, anno_type=None, name=None):
        self._tables = []  # 记录存入了哪些表格
        if task_id is None:
            task_id = self.bind_object.sheet.id
        if not isinstance(env_id, ObjectId) and env_id is not None:
            env_id = ObjectId(env_id)
        else:
            if 'env_id' in self.bind_object.sheet.data['options']:
                env_id = ObjectId(self.bind_object.option('env_id'))  # 仅仅即时计算直接绑定workflow对象
        _main_collection = self.db['beta_diversity']
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
                    self.bind_object.set_error('RDA/CCA分析没有生成正确的结果数据', code="51000202")
            main_id = _main_collection.insert_one(insert_mongo_json).inserted_id
        else:
            if not main_id:
                self.bind_object.set_error('不写入主表时，需要提供主表ID', code="51000203")
            if not isinstance(main_id, ObjectId):
                main_id = ObjectId(main_id)
            if analysis == 'rda_cca':  # 在主表中添加必要的rda或者是cca分类信息
                rda_files = os.listdir(dir_path.rstrip('/') + '/')
                if 'cca_sites.xls' in rda_files:
                    _main_collection.update_one({'_id': main_id}, {'$set': {'rda_cca': 'cca'}})
                elif 'rda_sites.xls' in rda_files:
                    _main_collection.update_one({'_id': main_id}, {'$set': {'rda_cca': 'rda'}})
                else:
                    self.bind_object.set_error('RDA/CCA分析没有生成正确的结果数据', code="51000202")
        result = _main_collection.find_one({'_id': main_id})
        if not result:
            self.bind_object.set_error('找不到主表id对应的表', code="51000204")
        if analysis == 'pca':
            site_path = dir_path.rstrip('/') + '/pca_sites.xls'
            rotation_path = dir_path.rstrip('/') + '/pca_rotation.xls'
            all_rotation_path = web_path.rstrip('/') + '/pca_rotation_all.xls'
            importance_path = dir_path.rstrip('/') + '/pca_importance.xls'
            _main_collection.update_one({'_id': main_id}, {'$set': {'download_file': all_rotation_path}})
            self.insert_table_detail(site_path, 'specimen', update_id=main_id)
            self.insert_table_detail(rotation_path, 'species', update_id=main_id, split_fullname=False)
            self.insert_table_detail(importance_path, 'importance', update_id=main_id, colls=['proportion_variance'])
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
            self.bind_object.logger.info('beta_diversity:PCoA分析结果导入数据库完成.')
        elif analysis == 'nmds':
            site_path = dir_path.rstrip('/') + '/nmds_sites.xls'
            self.insert_table_detail(site_path, 'specimen', update_id=main_id,
                                     colls=['NMDS1', 'NMDS2'])  # 这里暂且只有两个轴，故只写两个，后面如果修改进行多维NMDS，结果应该有多列，此种方式应该进行修改
            nmds_stress = round(float(open(dir_path.rstrip('/') + '/nmds_stress.xls').readlines()[1]),3)
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
            self.bind_object.logger.info('beta_diversity:db_RDA分析结果导入数据库完成.')
        elif analysis == 'rda_cca':
            # if 'rda' in os.listdir(dir_path.rstrip('/') + '/Rda/')[1]:
            if 'rda' in os.listdir(dir_path)[1]:
                rda_cca = 'rda'
            else:
                rda_cca = 'cca'
            site_path = dir_path.rstrip('/') + '/' + rda_cca + '_sites.xls'
            importance_path = dir_path.rstrip('/') + '/' + rda_cca + '_importance.xls'
            dca_path = dir_path.rstrip('/') + '/' + 'dca.xls'
            plot_species_path = dir_path.rstrip(
                '/') + '/' + rda_cca + '_plot_species_data.xls'  # add 4 lines by zhouxuan 20170123 20170401
            # envfit_path = dir_path.rstrip('/') + '/Rda/' + rda_cca + '_envfit.xls'
            # if os.path.exists(envfit_path):
            # self.insert_envfit_table(envfit_path, 'envfit', update_id=main_id)
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
                        self.bind_object.set_error('需要删除的匹配列必须为列表', code="51000205")
                    flag = 0
                    for i in fileter_biplot:
                        if re.match(r'^{}'.format(i), values[0]):
                            flag = 1
                    if flag:
                        continue
                else:
                    pass
                insert_data = {
                    'beta_diversity_id': update_id,
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
                collection.insert_many(data_temp)
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
                                  'plot_species': 'plot_species'}
                if table_type in default_column:
                    main_collection.update_one({'_id': update_id},
                                               {'$set': {default_column[table_type]: ','.join(columns)}},
                                               upsert=False)
                else:
                    self.bind_object.logger.error('错误的表格类型：%s不能在主表中插入相应表头' % table_type)
                    self.bind_object.set_error("主表更新失败", code="51000206")

    def insert_main_tables(self, tables, update_id, main='beta_diversity'):
        """
        """
        main_collection = self.db[main]
        main_collection.update_one({'_id': update_id},
                                   {'$set': {'tables': ','.join(tables)}},
                                   upsert=False)

    def insert_envfit_table(self, filepath, tabletype, update_id):
        """
        """
        insert_data = []
        with open(filepath, 'rb') as r:
            head = r.next().strip('\r\n')
            head = re.split(' ', head)
            new_head_old = head
            new_head = []
            for f in range(0, len(new_head_old)):
                pattern = re.compile('"(.*)"')
                new_head_ = pattern.findall(new_head_old[f])
                new_head.append(new_head_[0])
                # m = re.match('/\"(.*)\"/', new_head[f])
                # if m:
                #     new_head[f] = m.group(0)
            print(new_head)
            new_head[-1] = "p_values"
            new_head[-2] = "r2"
            self.new_head = new_head
            self.bind_object.logger.info(self.new_head)
            for line in r:
                line = line.rstrip("\r\n")
                line = re.split(' ', line)
                content = line[1:]
                pattern = re.compile('"(.*)"')
                env_name = pattern.findall(line[0])
                env_detail = dict()
                for i in range(0, len(content)):
                    env_detail[new_head[i]] = content[i]
                env_detail['beta_diversity_id'] = update_id
                env_detail['type'] = tabletype
                env_detail['env'] = env_name[0]
                insert_data.append(env_detail)
        try:
            collection = self.db['beta_diversity_detail']
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
