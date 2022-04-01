import requests
import re
from datetime import datetime
from urllib.request import urlopen
import json
import requests
import sys
import pandas as pd
import numpy as np
from  multiprocessing import Process,Pool
import os, time, random
import bs4

# ~/app/program/miniconda3/bin/python3
# ~/app/program/miniconda3/bin/pip install bs4
# ~/app/program/miniconda3/bin/pip install lxml
# ~/app/program/miniconda3/bin/pip install html5lib

# def ensembl_genesymbol_to_id(genename,speice):
#     # Homo sapiens, Mus musculus, Rattus norvegicus
#     speice_dic = {"Homo sapiens":"homo_sapiens", "Mus musculus":"mus_musculus", "Rattus norvegicus":"rattus_norvegicus"}
#     server = "https://rest.ensembl.org"
#     ext = "/lookup/symbol/{}/{}?".format(speice_dic[speice], genename)
#     r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
#     if not r.ok:
#         gene_id = "Not Found"
#         return gene_id
#         # r.raise_for_status()
#         # sys.exit()
#     else:
#         decoded = r.json()
#         # print(decoded['id'])
#         if "id" in decoded.keys():
#             gene_id = decoded['id']
#         else:
#             gene_id = "Not Found"
#         return gene_id

class mirnasearch(object):
    def __init__(self, speice, mirna, gene, outfile, miRTarBase_file, database):
        self.speice = speice
        self.mirna = mirna
        self.gene = gene
        self.outfile = outfile
        # self.miRTarBase_file = "/mnt/lustre/users/sanger-dev/sg-users/xuxi/tmp_2/miRTarBase_MTI.txt"
        self.miRTarBase_file = miRTarBase_file
        self.database = database

    def tar_base_search(self, speice="",mirna="",gene=""):
        if speice == "Human":
            speice_num = "1"
        if speice == "Mouse":
            speice_num = "5"
        if speice == "Rat":
            speice_num = "6"
        main_url = "http://carolina.imis.athena-innovation.gr/diana_tools/web/index.php?r=tarbasev8%2Findex&miRNAs%5B%5D={}&genes%5B%5D={}&species%5B%5D={}&sources%5B%5D=1&sources%5B%5D=7&sources%5B%5D=9&publication_year=&prediction_score=&sort_field=&sort_type=&query=".format(mirna,gene,speice_num)
        print(main_url)
        #################
        content = requests.get(main_url).content
        soup = bs4.BeautifulSoup(content, "html5lib")
        page_list = []
        page_items = soup.findAll("ul", {"class":"pagination"})
        for page_item in page_items:
            page_list = [int(child.text.strip()) for child in page_item.findAll('li') if child.text.strip().isdigit()]
        ####################
        all_raws = [["gene_name","mirna_name","experiments_throughput","publications","cell_lines","tissues","pred_score"]]
        export_items = soup.findAll("tr", {"class":"first-level"})
        for export_item in export_items:
            one_row = [child.text.strip() for child in export_item.findAll('td')][:7]
            all_raws.append(one_row)
        ####################
        if page_list:
            for page in page_list[1:]:
                detail_main_url = main_url + "&page=" +str(page)
                detail_content = requests.get(detail_main_url).content
                # print(content)
                detail_soup = bs4.BeautifulSoup(detail_content, "html5lib")
                export_items = detail_soup.findAll("tr", {"class":"first-level"})
                for export_item in export_items:
                    one_row = [child.text.strip() for child in export_item.findAll('td')][:7]
                    all_raws.append(one_row)
        #all_raws  结果形似 [['gene_name','mirna_name','experiments_throughput','publications','cell_lines','tissues','pred_score'],
        #                   ['RCOR1', 'hsa-miR-432-5p', 'low: 1 high: 0', '1', '1', '1', '-'],
        #                   ['IAPP', 'hsa-miR-432-5p', 'low: 1 high: 0', '1', '1', '1', '-']......]
        return all_raws
    # a = tar_base_search(speice="Human",mirna="hsa-miR-432-5p",gene="")

    def mirdb_search(self, speice="",mirna="",gene=""):
        # speice = "Human"
        # mirna = "" #"hsa-let-7e-5p"
        # gene = "HMGA2" #"HMGA2"
        if gene:
            searchType = "gene"
            searchBox = gene
        if mirna:
            searchType = "miRNA"
            searchBox = mirna
        with requests.Session() as session:
            headers_ = {
                "Host":"mirdb.org",
                "User-Agent":"Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:88.0) Gecko/20100101 Firefox/88.0",
                "Accept":"text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9",
                "Accept-Language":"zh-CN,zh;q=0.9,en-GB;q=0.8,en;q=0.7",
                "Accept-Encoding":"gzip, deflate",
                "Content-Type":"application/x-www-form-urlencoded",
                "Content-Length":"70",
                "Origin":"http://mirdb.org",
                "Connection":"keep-alive",
                "Referer":"http://mirdb.org/index.html",
                "Upgrade-Insecure-Requests":"1" 
            }
            paramsPost = {'species': speice, # Human, Mouse ,Rat
                    'searchBox': searchBox,
                    'submitButton': "Go",
                    "searchType": searchType
                    }
            if searchType == "gene":
                paramsPost["geneChoice"] = "symbol"
            response = session.post("http://mirdb.org/cgi-bin/search.cgi", data=paramsPost, headers=headers_)
            page_lists = []
            if response.status_code == 200:
                soup = bs4.BeautifulSoup(response.content, "html5lib")
                page_items = soup.findAll("tr")
                for page_item in page_items:
                    page_list = [child.text.strip() for child in page_item.findAll('td')]
                    if len(page_list) == 6:
                        page_lists.append(page_list)
            # page_lists 类似如下结构
            # [['Target Detail','Target Rank','Target Score','miRNA Name','Gene Symbol','Gene Description'],
            #  ['Details','1','100','hsa-miR-4458','HMGA2','high mobility group AT-hook 2']......]
            return page_lists

    def mirtarbase_search(self, speice="",mirna="",gene="",miRTarBase_file=""):
        speice_dic = {"Human":"Homo sapiens", "Mouse":"Mus musculus", "Rat":"Rattus norvegicus"}
        # speice = "Rat"
        # mirna = "hsa-let-7e-5p" #"hsa-let-7e-5p"
        # gene = "" #"HMGA2"
        mirtarbase_df = pd.read_table(miRTarBase_file,sep="\t",dtype={"Target Gene (Entrez ID)":"object"})
        if speice:
            mirtarbase_df_2 = mirtarbase_df.loc[mirtarbase_df['Species (miRNA)'] == speice_dic[speice]]
        else:
            mirtarbase_df_2 = mirtarbase_df
        if mirna:
            mirtarbase_df_3 = mirtarbase_df_2.loc[mirtarbase_df_2['miRNA'] == mirna]
        else:
            mirtarbase_df_3 = mirtarbase_df_2
        if gene:
            mirtarbase_df_4 = mirtarbase_df_3.loc[mirtarbase_df_3['Target Gene'] == gene]
        else:
            mirtarbase_df_4 = mirtarbase_df_3
        mirtarbase_df_5 = mirtarbase_df_4.iloc[:,[0,1,3,4]]
        mirtarbase_df_5.columns = ["mature_mirna_acc","mature_mirna_id","target_symbol","target_entrez"]
        mirtarbase_df_5["database"] = "mirtarbase"
        mirtarbase_df_6 = mirtarbase_df_5.drop_duplicates()
        return mirtarbase_df_6

    # def ncbi_genesymbol_to_id(self, genename,speice):
    #     # Homo sapiens, Mus musculus, Rattus norvegicus
    #     # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=HOXA9+Mus%20musculus&sort=relevance&retmax=1&retmode=json
    #     url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={}+{}&sort=relevance&retmax=1&retmode=json".format(genename,speice.replace(" ","%20"))
    #     html_dict = json.loads(urlopen(url).read())
    #     if int(html_dict["esearchresult"]["count"]) > 0:
    #         gene_id = html_dict["esearchresult"]["idlist"][0]
    #     else:
    #         gene_id = "Not Found"
    #     return gene_id

    # def ncbi_genesymbol_to_id(self,genename,speice):
    #     url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={}+{}&sort=relevance&retmax=1&retmode=json".format(genename,speice.replace(" ","%20"))
    #     r = requests.get(url, headers={ "Content-Type" : "application/json"})
    #     if not r.ok:
    #         gene_id = "Not Found1"
    #         r.raise_for_status()
    #         # sys.exit()
    #         return gene_id
    #     else:
    #         decoded = r.json()
    #         if int(decoded["esearchresult"]["count"]) > 0:
    #             gene_id = decoded["esearchresult"]["idlist"][0]
    #         else:
    #             gene_id = "Not Found2"
    #         return gene_id

    def ncbi_genesymbol_to_id(self,genename,speice,df=""):
        # gene = "hsd3b7"
        # speice = "Rat"
        # df = pd.read_table(human_mouse_rat_id_db_file,sep="\t")
        a = df[(df['external_gene_name'].str.lower() == genename.lower()) &(df['speice'].str.lower() == speice.lower())]
        if a.shape[0] > 0:
            ncbi_id = str(a.iloc[0,:].to_dict()["entrezgene_id"])
        else:
            ncbi_id = "None"
        return ncbi_id

    def ensembl_genesymbol_to_id(self,genename,speice,df=""):
        # gene = "hsd3b7"
        # speice = "Rat"
        # df = pd.read_table(human_mouse_rat_id_db_file,sep="\t")
        a = df[(df['external_gene_name'].str.lower() == genename.lower()) &(df['speice'].str.lower() == speice.lower())]
        if a.shape[0] > 0:
            ensembl_id = str(a.iloc[0,:].to_dict()["ensembl_gene_id"])
        else:
            ensembl_id = "None"
        return ensembl_id

    # def ensembl_genesymbol_to_id(self, genename,speice):
    #     # Homo sapiens, Mus musculus, Rattus norvegicus
    #     speice_dic = {"Homo sapiens":"homo_sapiens", "Mus musculus":"mus_musculus", "Rattus norvegicus":"rattus_norvegicus"}
    #     server = "https://rest.ensembl.org"
    #     ext = "/lookup/symbol/{}/{}?".format(speice_dic[speice], genename)
    #     r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    #     if not r.ok:
    #         gene_id = "Not Found"
    #         return gene_id
    #         # r.raise_for_status()
    #         # sys.exit()
    #     else:
    #         decoded = r.json()
    #         # print(decoded['id'])
    #         if "id" in decoded.keys():
    #             gene_id = decoded['id']
    #         else:
    #             gene_id = "Not Found"
    #         return gene_id

    def run(self):
        if self.database == "all":
            tar_base_result = self.tar_base_search(speice=self.speice,mirna=self.mirna,gene=self.gene)
            headers = tar_base_result.pop(0)
            df = pd.DataFrame(tar_base_result, columns=headers)
            tarbase_df = df.iloc[:,[0,1]]
            tarbase_df.columns = ["target_symbol","mature_mirna_id"]
            tarbase_df["database"] = "tarbase"
            tarbase_df = tarbase_df.drop_duplicates()
            #
            mirdb_result = self.mirdb_search(speice=self.speice,mirna=self.mirna,gene=self.gene)
            headers = mirdb_result.pop(0)
            df = pd.DataFrame(mirdb_result, columns=headers)
            mirdb_df = df.iloc[:,[3,4]]
            mirdb_df.columns = ["mature_mirna_id","target_symbol"]
            mirdb_df["database"] = "mirdb"
            mirdb_df = mirdb_df.drop_duplicates()
            #
            mirtarbase_df = self.mirtarbase_search(speice=self.speice,mirna=self.mirna,gene=self.gene,miRTarBase_file=self.miRTarBase_file)
            #
            web_result_df = pd.concat([tarbase_df,mirdb_df,mirtarbase_df])
            web_result_df["target_ensembl"] = np.nan
            return web_result_df
            # #
            # ensembl_id_result = {}
            # speice_dic = {"Human":"Homo sapiens", "Mouse":"Mus musculus", "Rat":"Rattus norvegicus"}
            # pool = Pool(5) #创建一个5个进程的进程池
            # for i in set(web_result_df["target_symbol"]):
            #     ensembl_id_result[i] = pool.apply_async(func=ensembl_genesymbol_to_id, args=(i,speice_dic[self.speice],))
            # pool.close()
            # pool.join()
            # print('并行结束')
            # ensembl_id_result_all = {i:ensembl_id_result[i].get() for i in ensembl_id_result}
            # #
            # ncbi_id_result = {}
            # for i in set(web_result_df["target_symbol"]):
            #     time.sleep(0.5)
            #     ncbi_id_result[i] = ncbi_genesymbol_to_id(i,"Homo sapiens")
            # #
            # web_result_df['target_ensembl'] = web_result_df['target_symbol'].map(ensembl_id_result_all)
            # web_result_df['target_entrez'] = web_result_df['target_symbol'].map(ncbi_id_result)
            # web_result_df = web_result_df.replace(np.nan, 'None', regex=True)
            # web_result_df = web_result_df.replace("Not Found", 'None', regex=True)
            # web_result_df = web_result_df[['database','mature_mirna_acc','mature_mirna_id','target_symbol','target_entrez','target_ensembl']]
            # web_result_df.to_csv(self.outfile, sep='\t',index=False)
        if self.database == "tarbase":
            tar_base_result = self.tar_base_search(speice=self.speice,mirna=self.mirna,gene=self.gene)
            headers = tar_base_result.pop(0)
            df = pd.DataFrame(tar_base_result, columns=headers)
            tarbase_df = df.iloc[:,[0,1]]
            tarbase_df.columns = ["target_symbol","mature_mirna_id"]
            tarbase_df["database"] = "tarbase"
            tarbase_df = tarbase_df.drop_duplicates()
            #
            web_result_df = tarbase_df
            web_result_df["target_ensembl"] = np.nan
            web_result_df["mature_mirna_acc"] = np.nan
            web_result_df["target_entrez"] = np.nan
            return web_result_df
        if self.database == "mirdb":
            mirdb_result = self.mirdb_search(speice=self.speice,mirna=self.mirna,gene=self.gene)
            headers = mirdb_result.pop(0)
            df = pd.DataFrame(mirdb_result, columns=headers)
            mirdb_df = df.iloc[:,[3,4]]
            mirdb_df.columns = ["mature_mirna_id","target_symbol"]
            mirdb_df["database"] = "mirdb"
            mirdb_df = mirdb_df.drop_duplicates()
            #
            web_result_df = mirdb_df
            web_result_df["target_ensembl"] = np.nan
            web_result_df["mature_mirna_acc"] = np.nan
            web_result_df["target_entrez"] = np.nan
            return web_result_df
        if self.database == "mirtarbase":
            mirtarbase_df = self.mirtarbase_search(speice=self.speice,mirna=self.mirna,gene=self.gene,miRTarBase_file=self.miRTarBase_file)
            #
            web_result_df = mirtarbase_df
            web_result_df["target_ensembl"] = np.nan
            return web_result_df


if __name__ == "__main__":
    import argparse
    print("参数解析开始")
    parser = argparse.ArgumentParser(description='This script is used to search miRNA-Gene from 5 database.')
    parser.add_argument('-speice', type=str, default="Human", help='speice name: Human/Mouse/Rat')
    parser.add_argument('-mirna', type=str, default="", help='miRNA name')
    parser.add_argument('-gene', type=str, default="",help='gene name')
    parser.add_argument('-outfile', type=str, help='output file path')
    parser.add_argument('-mirtarbase', type=str, help='miRTarBase file path')
    parser.add_argument('-idmappdb', type=str, help='idmapping database file path，parsed from linbinxu database sanger_biodb sgdb_*_gene')
    parser.add_argument('-database', type=str, default="all", help='which database: tarbase/mirdb/mirtarbase/all')
    args = parser.parse_args()
    mirnasearch = mirnasearch(args.speice, args.mirna, args.gene, args.outfile, args.mirtarbase, args.database)
    print("参数解析结束："+str(args))
    print("run函数开始")
    web_result_df = mirnasearch.run()
    print("run函数结束")

    #
    iddf = pd.read_table(args.idmappdb,sep="\t")
    #
    print('ensembl开始搜索')
    ensembl_id_result = {}
    for i in set(web_result_df["target_symbol"]):
        ensembl_id_result[i] = mirnasearch.ensembl_genesymbol_to_id(i,args.speice,df=iddf)
    print('ensembl搜索结束')
    #
    #
    print('ncbi开始搜索')
    ncbi_id_result = {}
    for i in set(web_result_df["target_symbol"]):
        ncbi_id_result[i] = mirnasearch.ncbi_genesymbol_to_id(i,args.speice,df=iddf)
    print('ncbi搜索结束')
    #
    # #
    # print('ensembl并行开始')
    # ensembl_id_result = {}
    # speice_dic = {"Human":"Homo sapiens", "Mouse":"Mus musculus", "Rat":"Rattus norvegicus"}
    # pool = Pool(5) #创建一个5个进程的进程池
    # for i in set(web_result_df["target_symbol"]):
    #     ensembl_id_result[i] = pool.apply_async(func=ensembl_genesymbol_to_id, args=(i,speice_dic[args.speice],))
    # pool.close()
    # pool.join()
    # print('ensembl并行结束')
    # ensembl_id_result_all = {i:ensembl_id_result[i].get() for i in ensembl_id_result}
    # #
    # print('ncbi串行开始')
    # ncbi_id_result = {}
    # for i in set(web_result_df["target_symbol"]):
    #     time.sleep(1)
    #     ncbi_id_result[i] = mirnasearch.ncbi_genesymbol_to_id(i,"Homo sapiens")
    # print('ncbi串行结束')
    # #
    web_result_df['target_ensembl'] = web_result_df['target_symbol'].map(ensembl_id_result)
    web_result_df['target_entrez'] = web_result_df['target_symbol'].map(ncbi_id_result)
    web_result_df = web_result_df.replace(np.nan, 'None', regex=True)
    web_result_df = web_result_df.replace("Not Found", 'None', regex=True)
    web_result_df = web_result_df[['database','mature_mirna_acc','mature_mirna_id','target_symbol','target_entrez','target_ensembl']]
    web_result_df.to_csv(args.outfile, sep='\t',index=False)
    print('成功导出结果文件')