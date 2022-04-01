#/usr/bin/python
#__author__ = 'guanqing.zou'
#20181121


from biocluster.config import Config
from bson.objectid import ObjectId
import argparse
import copy

class FunctionDetail:
    def __init__(self,database,level,member_list, table_id, mg_dbversion=None):
        Config().DBVersion = mg_dbversion
        self.client = Config().get_mongo_client(mtype="metagenomic", ref=True)
        self.mongodb = self.client[Config().get_mongo_dbname("metagenomic", ref=True)]
        self.client_1 = Config().get_mongo_client(mtype="metagenomic")
        self.mongodb_func_set =  self.client_1[Config().get_mongo_dbname("metagenomic")]
        self.database = database
        self.level = level
        self.member = member_list
        #self.out_file = out_file
        if not isinstance(table_id, ObjectId):
            table_id = ObjectId(table_id)
        self.table_id = table_id
        self.collection = 'func_set'
        self.update_dic = None

    def get_collection_info(self, collection, search_key, values):
        ret_list = []
        for one in self.mongodb[collection].find({search_key: {"$in": values}}):
            ret_list.append(one)
        return ret_list

    def get_collection_mixture(self, collection, search_key, values):
        ret_list = []
        for one in self.mongodb[collection].find({search_key: {"$in": values}}):
            ret_list.append(one)
        return ret_list

    def get_set_values(self, ret_list,set_k):
        ret_dict = {}
        for k in set_k:
            ret_dict[k] = []
        for i in ret_list:
            for k in set_k:
                if k not in i:
                    continue
                if isinstance(i[k],list):
                    i[k] = [v for v in i[k] if v not in ['', '-']]
                    ret_dict[k].extend(i[k])
                else:
                    if i[k] not in ['', '-']:
                        ret_dict[k].append(i[k])
        ret_dict = {k: list(set(v)) for k, v in ret_dict.items()}
        return ret_dict


    def update_collection_info(self, update_dic):

        #print self.mongodb_func_set[self.collection]
        update_dic['anno_type'] = self.database
        update_dic['level'] = self.level + '_list'
        update_dic['member'] = self.member
        t = self.mongodb_func_set[self.collection].find_one({"_id":self.table_id})
        self.mongodb_func_set[self.collection].update({"_id":self.table_id},{"$set": update_dic},upsert=True)
        self.mongodb_func_set[self.collection].update({"_id":self.table_id},{"$set": {'status':'end'}})

    def deal_kegg(self):
        update_dic = {}
        if self.level == 'pathway':
            update_dic['pathway_list'] = self.member
            results = self.get_collection_info('kegg_pathway_info', "pathway_id", self.member)
            update_dic['module_list'] = self.get_set_values(results, ['module_list'])['module_list']
            update_dic['enzyme_list'] = self.get_set_values(results,['enzyme_id'])['enzyme_id']
            #k_result = self.get_collection_info('kegg_module_info','module_id', update_dic['module_list'])
            k_result = self.get_collection_mixture('kegg_ko_v1','pathway_id', self.member)
            update_dic['k_list'] = self.get_set_values(k_result,['ko_id'])['ko_id']
            #results = self.get_collection_info('kegg_ko_v1', "ko_id", update_dic['k_list'])


        elif self.level == 'module':
            update_dic['module_list'] = self.member
            results = self.get_collection_info('kegg_module_info', "module_id", self.member)
            set_result = self.get_set_values(results,['pathway_list','k_list', 'enzyme_id'])
            update_dic['pathway_list'] = set_result['pathway_list']
            update_dic['k_list'] = set_result['k_list']
            update_dic['enzyme_list'] = set_result['enzyme_id']
            #results = self.get_collection_info('kegg_ko_v1', "ko_id", update_dic['k_list'])
            #update_dic['enzyme_list'] = self.get_set_values(results,['enzyme_id'])['enzyme_id']


        elif self.level == 'k' or self.level == 'enzyme' or self.level == 'gene':
            if self.level == 'gene':
                update_dic['k_list'] = self.member
                results = self.get_collection_info('kegg_ko_v1', "ko_id", self.member)
                set_result = self.get_set_values(results,['enzyme_id','module_id','pathway_id'])
                update_dic['enzyme_list'] = set_result['enzyme_id']
                update_dic['pathway_list'] = set_result['pathway_id']
                update_dic['module_list'] = set_result['module_id']
            elif self.level == 'k':
                update_dic['k_list'] = self.member
                results = self.get_collection_info('kegg_ko_v1', "ko_id", self.member)
                set_result = self.get_set_values(results,['enzyme_id','module_id','pathway_id'])
                update_dic['enzyme_list'] = set_result['enzyme_id']
                update_dic['pathway_list'] = set_result['pathway_id']
                update_dic['module_list'] = set_result['module_id']
            else:
                update_dic['enzyme_list'] = self.member
                results = self.get_collection_info('kegg_enzyme_info', "enzyme_id", self.member)
                set_result = self.get_set_values(results,['module_list','pathway_list'])
                update_dic['pathway_list'] = set_result['pathway_list']
                update_dic['module_list'] = set_result['module_list']
                results = self.get_collection_mixture('kegg_ko_v1', "enzyme_id", self.member)
                update_dic['k_list'] = self.get_set_values(results,['ko_id'])['ko_id']

        self.update_dic = update_dic

    def deal_eggnog(self):
        update_dic = {}
        if self.level == 'category':
            update_dic['category_list'] = self.member
            results = self.get_collection_info('eggNOG4_category_info', "category", self.member)
            set_results = self.get_set_values(results, ['function','nog'])
            update_dic['function_list'] = set_results['function']
            update_dic['nog_list'] = set_results['nog']
        elif self.level == 'function':
            update_dic['function_list'] = self.member
            results = self.get_collection_info('eggNOG4_function_info', "function", self.member)
            set_results = self.get_set_values(results, ['category','nog'])
            update_dic['category_list'] = set_results['category']
            update_dic['nog_list'] = set_results['nog']
        elif self.level == 'nog':
            update_dic['nog_list'] = self.member
            results = self.get_collection_info('eggNOG4', "nog", self.member)
            set_results = self.get_set_values(results, ['category','function'])
            update_dic['category_list'] = set_results['category']
            update_dic['function_list'] = set_results['function']

        self.update_dic = update_dic


    def deal_cazy(self):
        update_dic = {}
        if self.level == 'class':
            update_dic['class_list'] = self.member
            results = self.get_collection_info('cazy_class_info', "class", self.member)
            set_results = self.get_set_values(results, ['family'])
            update_dic['family_list'] = set_results['family']

        elif self.level == 'family':
            update_dic['family_list'] = self.member
            results = self.get_collection_info('cazy', "family", self.member)
            set_results = self.get_set_values(results, ['class'])
            update_dic['class_list'] = set_results['class']

        self.update_dic = update_dic

    def deal_card(self):
        update_dic = {}
        if self.level == 'class':
            update_dic['class_list'] = self.member
            results = self.get_collection_info('card_class_info', 'class', self.member)
            set_results = self.get_set_values(results, ['aro_accession'])
            update_dic['aro_list'] = set_results['aro_accession']
        elif self.level == 'aro':
            update_dic['aro_list'] = self.member
            results = self.get_collection_info('card','aro_accession', self.member)
            set_results = self.get_set_values(results, ['class'])
            update_dic['class_list'] = set_results['class']

        self.update_dic =update_dic

    def deal_vfdb(self):
        update_dic = {}
        if self.level == 'level1':
            update_dic['level1_list'] = self.member
            results = self.get_collection_info('vfdb_level1_info','level1',self.member)
            set_results = self.get_set_values(results, ['level2','vfs'])
            update_dic['level2_list'] = set_results['level2']
            update_dic['vfs_list'] = set_results['vfs']
        elif self.level == 'level2':
            update_dic['level2_list'] = self.member
            results = self.get_collection_info('vfdb_level2_info','level2',self.member)
            set_results = self.get_set_values(results, ['level1','vfs'])
            update_dic['level1_list'] = set_results['level1']
            update_dic['vfs_list'] = set_results['vfs']
        elif self.level == 'vfs':
            update_dic['vfs_list'] = self.member
            results = self.get_collection_info('vfdb_vfs_info','vfs',self.member)
            set_results = self.get_set_values(results, ['level1','level2','vf'])
            update_dic['level1_list'] = set_results['level1']
            update_dic['level2_list'] = set_results['level2']
        self.update_dic = update_dic

    def deal_ardb(self):
        update_dic = {}
        if self.level == 'class':
            update_dic['class_list'] = self.member
            results = self.get_collection_mixture('ardb','class',self.member)
            set_results = self.get_set_values(results, ['type','arg', 'resistance'])
            update_dic['type_list'] = set_results['type']
            update_dic['arg_list'] = set_results['arg']
            update_dic['resistance_list'] = set_results['resistance']
        elif self.level == 'type':
            update_dic['type_list'] = self.member
            results = self.get_collection_mixture('ardb','type',self.member)
            set_results = self.get_set_values(results, ['class','arg', 'resistance'])
            update_dic['class_list'] = set_results['class']
            update_dic['arg_list'] = set_results['arg']
            update_dic['resistance_list'] = set_results['resistance']
        elif self.level == 'arg':
            update_dic['arg_list'] = self.member
            results = self.get_collection_mixture('ardb','arg',self.member)
            set_results = self.get_set_values(results, ['class','type', 'resistance'])
            update_dic['class_list'] = set_results['class']
            update_dic['type_list'] = set_results['type']
            update_dic['resistance_list'] = set_results['resistance']
        elif self.level == 'resistance':
            update_dic['resistance_list'] = self.member
            results = self.get_collection_info('ardb_resistance_info','resistance',self.member)

            set_results = self.get_set_values(results, ['type','class','arg'])
            update_dic['type_list'] = set_results['type']
            #type_results  = self.get_collection_mixture('ardb','type',set_results['type'])
            #type_set_results = self.get_set_values(type_results, ['class','arg'])
            update_dic['class_list'] = set_results['class']
            update_dic['arg_list'] = set_results['arg']

        self.update_dic = update_dic


    def main_fun(self):
        if self.database == 'kegg':
            self.deal_kegg()
        elif  self.database == 'cog':
            self.deal_eggnog()
        elif  self.database == 'cazy':
            self.deal_cazy()
        elif  self.database == 'card':
            self.deal_card()
        elif  self.database == 'vfdb':
            self.deal_vfdb()
        elif  self.database == 'ardb':
            self.deal_ardb()
        else:
            raise Exception('database name is wrong')
        self.update_collection_info(self.update_dic)





if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--database", required=True)
    parser.add_argument("-l", "--level", required=True)
    parser.add_argument("-m", "--member", required=True)
    #parser.add_argument("-o", "--outfile", required=True)
    parser.add_argument("-t", "--table", required=True)
    parser.add_argument("-v", "--db_version", default=None)  # mongo db_version
    args = parser.parse_args()

    database = args.database
    level = args.level.split('_list')[0]
    member_1 = args.member.strip('"')
    member = member_1.split(',')
    #outfile = args.outfile
    table = args.table
    cls = FunctionDetail(database, level, member, table, args.db_version)
    cls.main_fun()















