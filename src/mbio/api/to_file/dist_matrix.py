# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
import os
import copy
from biocluster.config import Config
from bson.objectid import ObjectId

client = Config().get_mongo_client(mtype="meta")
db = client[Config().get_mongo_dbname("meta")]


def export_distance_matrix(data, option_name, dir_path, bind_obj=None):
    # db = Config().mongo_client[Config().MONGODB]
    file_path = os.path.join(dir_path, "%s_input.matrix.xls" % option_name)
    bind_obj.logger.debug("正在导出参数%s的距离矩阵为文件，路径:%s" % (option_name, file_path))
    collection = db['sg_beta_specimen_distance_detail']
    results = collection.find({"specimen_distance_id": ObjectId(data)})
    if results.count() == 0:
        raise Exception('距离矩阵id没有找到对应的detail数据')
    samples = []
    copy_results = copy.deepcopy(results)
    bind_obj.logger.info(str(results))
    for result in results:
        samples.append(result["specimen_name"])
    bind_obj.logger.info('ALL SAMPLE:' + ' '.join(samples))
    copysamples = copy.deepcopy(samples)
    with open(file_path, "wb") as f:
        f.write("\t%s\n" % "\t".join(samples))
        for sample in samples:
            doc = {}
            value = []
            for result in copy_results:
                if result['specimen_name'] == sample:
                    doc = result
                    break
            for detail in copysamples:
                value.append(doc[detail])
            value = [str(i) for i in value]
            f.write(sample + '\t' + '\t'.join(value) + '\n')
    return file_path
