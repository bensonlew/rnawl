# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import os
import re
from biocluster.config import Config
from bson.objectid import ObjectId

client = Config().get_mongo_client(mtype="metaasv")
db = client[Config().get_mongo_dbname("metaasv")]


def export_env_table(data, option_name, dir_path, bind_obj=None):
    # db = Config().mongo_client[Config().MONGODB]
    file_path = os.path.join(dir_path, '%s_input_env.xls' % option_name)
    bind_obj.logger.debug('正在导出参数%s的环境因子表为文件，路径:%s' % (option_name, file_path))
    collection_main = db['env']
    result_main = collection_main.find_one({'_id': ObjectId(data)})
    if not result_main:
        raise Exception('环境因子id没有找到对应的表信息')
    all_envs = result_main['env_names'].strip().split(',')
    collection = db['env_detail']
    specimen_collection = db['specimen_detail']
    results = collection.find({'env_id': ObjectId(data)})
    if results.count() == 0:
        raise Exception('环境因子id没有找到对应的detail数据')
    bind_obj.logger.info('ALL ENVS:' + ' '.join(all_envs))
    with open(file_path, 'wb') as f:
        f.write('#SampleID\t' + '\t'.join(all_envs) + '\n')
        for one in results:
            specimen_name = specimen_collection.find_one({'_id': one['specimen_id']})
            if specimen_name:
                line_list = [specimen_name['specimen']]
                for env in all_envs:
                    if re.search(r' ', str(one[env])):
                        new_env = str(one[env]).replace(" ", "")
                        line_list.append(new_env)
                    else:
                        line_list.append(str(one[env]))
                f.write('\t'.join(line_list) + '\n')
            else:
                raise Exception('样本id:%s在specimen_detail中没有找到' % str(one['specimen_id']))
    return file_path


def export_env_table_by_detail(data, option_name, dir_path, bind_obj=None):
    """
    env_labs group_detail必须存在options中，依据这两个值筛选环境因子和样本
    """
    pass


def export_float_env(data, option_name, dir_path, bind_obj=None):
    # db = Config().mongo_client[Config().MONGODB]
    file_path = os.path.join(dir_path, '%s_input_env.xls' % option_name)
    bind_obj.logger.debug('正在导出参数%s的环境因子表为文件，路径:%s' % (option_name, file_path))
    collection_main = db['env']

    bind_obj.logger.debug(bind_obj.sheet.option("env_labs"))
    env_labs = bind_obj.sheet.option("env_labs").split(',')

    result_main = collection_main.find_one({'_id': ObjectId(data)})
    if not result_main:
        raise Exception('环境因子id没有找到对应的表信息')
    all_envs = result_main['env_names'].strip().split(',')
    collection = db['env_detail']
    specimen_collection = db['specimen_detail']
    results = collection.find({'env_id': ObjectId(data)})
    if results.count() == 0:
        raise Exception('环境因子id没有找到对应的detail数据')
    bind_obj.logger.info('ALL ENVS:' + ' '.join(all_envs))
    bind_obj.logger.info('SELECT ENVS:' + ' '.join(env_labs))
    # write_lines = ['#SampleID\t' + '\t'.join(env_labs) + '\n']
    # flit_envs = set()
    with open(file_path, 'wb') as f:
        f.write('#SampleID\t' + '\t'.join(env_labs) + '\n')
        for one in results:
            specimen_name = specimen_collection.find_one({'_id': one['specimen_id']})
            if specimen_name:
                line_list = [specimen_name['specimen']]
                for env in env_labs:
                    line_list.append(str(one[env]))
                f.write('\t'.join(line_list) + '\n')
            else:
                raise Exception('样本id:%s在specimen_detail中没有找到' % str(one['specimen_id']))
    return file_path
    # with open(file_path, 'wb') as f:
    #     for one in results:
    #         specimen_name = specimen_collection.find_one({'_id': one['specimen_id']})
    #         if specimen_name:
    #             line_list = [specimen_name['specimen_name']]
    #             for env in env_labs:
    #                 # line_list.append(str(one[env]))
    #                 if not re.match(r"\d+?\.?\d+$", str(one[env])):
    #                     print("tttttt")
    #                     print(env)
    #                     flit_envs.add(env)
    #                     continue
    #                 else:
    #                     line_list.append(str(one[env]))
    #             write_lines.append('\t'.join(line_list) + '\n')
    #         else:
    #             raise Exception('样本id:%s在sg_specimen中没有找到' % str(one['specimen_id']))
    #     first_write_line = write_lines[0].strip().split("\t")
    #     bind_obj.logger.info(first_write_line)
    #     bind_obj.logger.info(flit_envs)
    #     for fe in flit_envs:
    #         if fe in first_write_line[1:]:
    #             print(fe)
    #             first_write_line.remove(fe)
    #     bind_obj.logger.info("\t".join(first_write_line)+"\n")
    #     f.write("\t".join(first_write_line)+"\n")
    #     for line in write_lines[1:]:
    #         f.write(line)
    # return file_path


###满足排序回归的环境因子表的格式要求
def export_float_env_regression(data, option_name, dir_path, bind_obj=None):
    file_path = os.path.join(dir_path, '%s_input_env.xls' % option_name)
    bind_obj.logger.debug('正在导出参数%s的环境因子表为文件，路径:%s' % (option_name, file_path))
    collection_main = db['env']

    bind_obj.logger.debug(bind_obj.sheet.option("env_labs"))
    env_labs = bind_obj.sheet.option("env_labs").split(',')

    result_main = collection_main.find_one({'_id': ObjectId(data)})
    if not result_main:
        raise Exception('环境因子id没有找到对应的表信息')
    all_envs = result_main['env_names'].strip().split(',')
    collection = db['env_detail']
    specimen_collection = db['specimen_detail']
    results = collection.find({'env_id': ObjectId(data)})
    if results.count() == 0:
        raise Exception('环境因子id没有找到对应的detail数据')
    bind_obj.logger.info('ALL ENVS:' + ' '.join(all_envs))
    bind_obj.logger.info('SELECT ENVS:' + ' '.join(env_labs))
    with open(file_path, 'wb') as f:
        f.write('\t' + '\t'.join(env_labs) + '\n')   #和export_float_env 不同之处
        for one in results:
            specimen_name = specimen_collection.find_one({'_id': one['specimen_id']})
            if specimen_name:
                line_list = [specimen_name['specimen']]
                for env in env_labs:
                    line_list.append(str(one[env]))
                f.write('\t'.join(line_list) + '\n')   #和export_float_env 不同之处
            else:
                raise Exception('样本id:%s在specimen_detail中没有找到' % str(one['specimen_id']))
    return file_path