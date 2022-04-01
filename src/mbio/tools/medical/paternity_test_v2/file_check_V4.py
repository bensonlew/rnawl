#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
作者：kefei.huang
时间：20171226

这个程序是为了合并三个数据表格而写的

1. 先把两个表格合并。线上表格是必须会有的，线下可能没有
2。随后对数据进行检查
3. 再调用批次表，得到加急非加急
4. 再调用数据库进行查重
5. 最后导表

注意下 @非哥 ：由于该文件你那边好像用是tab代替空格 所以 我都统一改成了空格。

"""
from __future__ import print_function
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
from bson.objectid import ObjectId
import os
import re
import codecs
import chardet
import MySQLdb as mdb
import datetime


class FileCheckV4Agent(Agent):
    def __init__(self,parent):
        super(FileCheckV4Agent,self).__init__(parent)
        self.logger.info("run")
        options = [
            {"name":"up","type":"string"},
            {"name":"down","type":"string","default":"empty"},
            {"name":"date","type":"string"},
            {"name":"batch","type":"string"},
            {"name":"flowcell","type":"string"},
            {"name":"main_id","type":"string"}
        ]
        self.add_option(options)
        self.logger.info("options finished")

    def check_options(self):
        if not self.option("up"):
            raise OptionError("请输入线上表格")
        if not self.option("date"):
            raise OptionError("请输入结果文件名称，以时间为单位，格式请参照20171023")
        if not self.option("batch"):
            raise OptionError("请输入实验批次表")
        if not self.option("flowcell"):
            raise OptionError("需要输入板号")
        if not self.option("main_id"):
            raise OptionError("需要输入主表ID")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = "1G"

    def end(self):
        super(FileCheckV4Agent,self).end()


class FileCheckV4Tool(Tool):
    def __init__(self, config):
        super(FileCheckV4Tool,self).__init__(config)
        self.up = self.option("up")
        self.down = self.option("down")
        self.date = self.option("date")
        self.batch = self.option("batch")
        self.flowcell = self.option("flowcell")
        self.update_info = self.option("main_id")
        self.interval = self.date + ".csv"
        self.interval_add_emergency = self.date + ".interval.csv"
        self.split = self.date + ".split.csv"
        self.errorlog = self.date + ".error.log"
        self.guest = self.date + ".guest.csv"
        self.customer = self.date + ".customer.csv"

        self.dbconfigure={'host':'10.100.200.135','name':'mjlimsbak','passwd':'Q6k9pU9ZN5FTow','database':'mjlims'}
        self.libraryType = {"NIPT文库构建实验流程":"WS", 
                            "WQ gDNA建库实验流程_多重":"WQ", 
                            "WQ cfDNA建库实验流程_杂捕":"WQ",
                            "FA gDNA建库实验流程_杂捕":"FA",
                            "FA ctDNA建库实验流程_杂捕":"FA",
                            "FA gDNA建库实验流程":"FA",
                            "FA cfDNA建库实验流程_杂捕":"FA",
                            "YCZL 建库实验流程":"YCZL",
                            "QA 建库实验流程":"QA",
                            "CT 建库实验流程":"CT",
                            "HLY 建库实验流程":"HLY,ZL",
                            "RXA 建库实验流程":"RXA",
                            "snp 建库实验流程":"SNP",
                            "XA RNA建库实验流程":"XA",
                            "XA DNA建库实验流程":"XA",
                            "甲基化 建库实验流程":"METH",
                            "转录组文库 建库实验流程":"TRSC",
                            "单细胞DNA 建库实验流程":"scWGS,sc,MOS",
                            "small RNA建库实验流程":"small",
                            "YGY 建库实验流程": "YG",
                            "全外白细胞建库": "QBC"}
        self.sampletype = {"全血": "QX", 
                            "血浆": "XJ", 
                            "蜡块": "SL", 
                            "石蜡": "SL", 
                            "石蜡切片": "SL", 
                            "石蜡切片(白片)": "SL,FP", 
                            "胸腹水" : "XS", 
                            "手术标本": "XZ,SS", 
                            "穿刺标本": "XZ", 
                            "穿刺样本": "XZ,CC",  # add by HD 20180131 添加CC
                            "组织标本": "XZ", 
                            "新鲜组织": "XZ",
                            "蜡卷": "SL,LJ",
                            "精斑": "JB",
                            "亲子父本全血":"QQ",
                            "指甲": "ZJ"}
        self.analysistype = {"NIPT文库构建实验流程":"nipt", 
                            "WQ gDNA建库实验流程_多重":"dcpt", 
                            "WQ cfDNA建库实验流程_杂捕":"pt",
                            "FA gDNA建库实验流程_杂捕":"ctdna",
                            "FA ctDNA建库实验流程_杂捕":"ctdna",
                            "FA gDNA建库实验流程":"ctdna",
                            "FA cfDNA建库实验流程_杂捕":"ctdna",
                            "YCZL 建库实验流程":"genetic_tumor",
                            "QA 建库实验流程":"QA_str",
                            "CT 建库实验流程":"ct_str",
                            "HLY 建库实验流程":"dc_str",
                            "XA RNA建库实验流程":"blood_cancer",
                            "XA DNA建库实验流程":"blood_cancer",
                            "甲基化 建库实验流程":"",
                            "RXA 建库实验流程":"",
                            "snp 建库实验流程":"",
                            "转录组文库 建库实验流程":"",
                            "单细胞DNA 建库实验流程":"",
                            "small RNA建库实验流程":"",
                            "YGY 建库实验流程": "",
                            "全外白细胞建库": ""}
        print("initialization finished")

    def run(self):
        super(FileCheckV4Tool, self).run()
        self.logger.info("tools start run")
        self.DicConvertUTF8()
        self.logger.info("convert finished")
        self.mergeupdown()
        self.logger.info("mergeupdown finished")

        Checkbarcode_Error = self.Checkbarcode()
        NamePattern_Error = self.NamePattern()
        Checkemergency_Error = self.Checkemergency()
        self.GenerateSplitTable()
        CheckSqlLIMS_Error = self.CheckSqlLIMS()
        #CheckMongoDupName_Error = self.CheckMongoDupName()
        ###判断是否可以导表
        self.logger.info(Checkbarcode_Error)
        self.logger.info(NamePattern_Error)
        self.logger.info(Checkemergency_Error)
        self.logger.info(CheckSqlLIMS_Error)
        Error_sum = Checkbarcode_Error + NamePattern_Error + Checkemergency_Error + CheckSqlLIMS_Error
        tempout = "Checkbarcode_Error = {}\nNamePattern_Error = {}\nCheckemergency_Error = {}\nCheckSqlLIMS_Error = {}\n".format(Checkbarcode_Error,NamePattern_Error,Checkemergency_Error,CheckSqlLIMS_Error)
        print (tempout)
        MAINID = self.update_info
        if Error_sum == 0:
            self.dump_error("end")
            self.dump_guest()
            self.dump_experiment_batch()
            self.dump_sample_info()
            self.dump_split()
            self.dump_sg_customer(MAINID)
        else:
            self.dump_error("failed")

        self.end()

    def DicConvertUTF8(self):
        for k,v in self.sampletype.items():
            new = k.encode("utf-8").decode("utf-8")
            self.sampletype.pop(k)
            self.sampletype[new] = v

        for k,v in self.libraryType.items():
            new = k.encode("utf-8").decode("utf-8")
            self.libraryType.pop(k)
            self.libraryType[new] = v

        for k,v in self.analysistype.items():
            new = k.encode("utf-8").decode("utf-8")
            self.analysistype.pop(k)
            self.analysistype[new] = v

    def Checkencoding(self,filepath):
        with open(filepath,"r") as f:
            title = f.readline()
            result = chardet.detect(title)
            if result['encoding'] != "utf-8":
                code = "gbk"
            else:
                code = "utf-8"
            return(code)

    def linestrip(self,line):
        line = line.strip()
        line = line.strip("\r")
        return(line)

    def mergeupdown(self):
        ###简单的用字典合并
        SC = "生产"
        SC = SC.encode("utf-8").decode("utf-8")
        mergeDic = {}
        count = 1
        if(self.up != "empty"):
            code_type = self.Checkencoding(self.up)
            with codecs.open (self.up,"r",code_type) as IN_UP: 
                title = self.linestrip(IN_UP.readline())
                titletemp = title.split(",")[0:14]
                titletemp[-1] = "备注"
                titletemp = ",".join(titletemp)
                mergeDic["title"] = titletemp
                for line in IN_UP:
                    line = self.linestrip(line)
                    if len(line) > 1:
                        linetemp = line.split(",")[0:14]
                        # if not re.search(r'^TEST',linetemp[6]):  # modified by hd 问题由黄克非提出，并与20180202电话与宣红东确认， 同意修改
                        linetemp[-1] = SC
                        linetemp = ",".join(linetemp)
                        mergeDic[count] = linetemp
                        count += 1

        if(self.down != "empty"):
            code_type = self.Checkencoding(self.down)
            with codecs.open (self.down,"r",code_type) as IN_DOWN:
                IN_DOWN.next()
                IN_DOWN.next() ###前面三行没有作用
                IN_DOWN.next()
                        
                for line in IN_DOWN:
                    line = self.linestrip(line)
                    if len(line) > 1:
                        linetemp = line.split(",")[0:11]
                        if linetemp[0] == "" and linetemp[1] == "":
                            break

                        # if not re.search(r'^TEST',linetemp[2]):  # modified by hd 问题由黄克非提出，并与20180202电话与宣红东确认， 同意修改
                        pushtemp = ["","",""]
                        pushtemp.extend(linetemp)
                        mergeDic[count] = ",".join(pushtemp)
                        count += 1

        self.logger.info("t table finished")

        with codecs.open(self.interval,"w","utf-8") as out:
            print(mergeDic["title"],file=out)
            for k in sorted(mergeDic.keys()):
                if k == "title" or mergeDic[k] == None:
                    continue
                else:
                    print(mergeDic[k],file = out)

    def Checkbarcode(self):
        """
            这一部分要处理简单的查重
            1. 首先检查barcode列是不是对的。
            2. 其次检查里面有没有重复
            3. 如果有重复，就要检查样品名是不是一致，如果一致就跳过。否则就要报错
        """
        error_count = 0
        BarcodeDic = {}
        code_type = self.Checkencoding(self.interval)
        ERROR = codecs.open(self.errorlog,"w","utf-8")
        with codecs.open(self.interval,"r",code_type) as IN:
            title = IN.readline()
            for line in IN:
                line = self.linestrip(line)
                linetemp = line.split(",")
                if re.match(r'^[ATCG]+$',linetemp[-3]) == None:
                    errorlog = "barcode consist of ATCG only,so barcode column may be error,check your file"
                    print (errorlog,file=ERROR)
                    error_count += 1
                if linetemp[-3] not in BarcodeDic:
                    BarcodeDic[linetemp[-3]] = [linetemp[6]]
                else:
                    BarcodeDic[linetemp[-3]].append(linetemp[6])
                    samplename = list(set(BarcodeDic[linetemp[-3]]))
                    if len(samplename) >= 1:
                        samplenameout = ",".join(samplename)
                        errorlog = samplenameout + " has same barcode " + linetemp[-3] + " but differ sample name"
                        error_count += 1
                        print (errorlog,file=ERROR)
        ERROR.close()
        return(error_count)

    def NamePattern(self):
        """
            这一部分就是处理杂捕，多重之类的命名规则。
            字典放置在程序的最上层。考虑到以后很多东西都会重新命名
        """
        pt_db = Config().get_mongo_client(mtype="pt_v2")[Config().get_mongo_dbname("pt_v2")]
        nipt_db = Config().get_mongo_client(mtype="nipt_v2")[Config().get_mongo_dbname("nipt_v2")]
        error_count = 0
        ERROR = codecs.open(self.errorlog,"a","utf-8")
        code_type = self.Checkencoding(self.interval)
        with codecs.open(self.interval,"r",code_type) as IN:
            title = IN.readline()
            for line in IN:
                line = self.linestrip(line)
                linetemp = line.split(",")

                ###检查样品名称中是否有多余的"."符号、空格符号
                if '.' in linetemp[6] or '.' in linetemp[7]:
                    error_count += 1
                    errorlog = '{}-{} has . in name!'.format(linetemp[6], linetemp[7])
                    print(errorlog, file=ERROR)
                if ' ' in linetemp[6] or ' ' in linetemp[7]:
                    error_count += 1
                    errorlog = '{}-{} has space symbol in name!'.format(linetemp[6], linetemp[7])
                    print(errorlog, file=ERROR)
                ###检查样品是否已经分析过
                sample_id = '{}-{}'.format(linetemp[6], linetemp[7])
                if sample_id.startswith('WQ'):
                    tab_result = pt_db["sg_sample_tab"].find_one({"sample_id": sample_id})
                    tab_result_problem = pt_db["sg_sample_tab_problem"].find_one({"sample_id": sample_id})
                    if tab_result or tab_result_problem:
                        error_count += 1
                        errorlog = '{}：样品已分析过，请检查是否命名错误！'.format(sample_id)
                        print(errorlog, file=ERROR)
                elif sample_id.startswith('WS'):
                    bed_result = nipt_db["sg_bed"].find_one({"sample_id": sample_id})
                    if bed_result:
                        error_count += 1
                        errorlog = '{}：样品已分析过，请检查是否命名错误！'.format(sample_id)
                        print(errorlog, file=ERROR)
                ###检查样品类型是不是对的
                if self.sampletype.has_key(linetemp[8]):
                    pattern = "|".join(self.sampletype[linetemp[8]].split(","))
                    if re.match(r'^TEST', linetemp[6]):
                        continue
                    if re.search(r'^'+pattern,linetemp[5]) == None:
                        errorlog = linetemp[5] + " do not have right sampletype " + linetemp[8] +" , check your files"
                        print (errorlog,file=ERROR)
                        error_count += 1
                # else:    # modified by hd 解决线下上机表 没有样本类型 流程就不会进行后面的计算 bug是由黄克非提出，并与20180202电话与宣红东确认，同意修改
                #     continue
                ###检查建库类型是不是对的
                if self.libraryType.has_key(linetemp[3]):
                    pattern = "|".join(self.libraryType[linetemp[3]].split(","))
                    if re.match(r'^TEST', linetemp[6]):
                        continue
                    if re.search(r'^'+pattern,linetemp[6]) == None:
                        errorlog = linetemp[5] + " with differ librarytype " + linetemp[3]
                        print (errorlog,file=ERROR)
                        error_count += 1
                else:
                    errorlog = linetemp[3] + " do not have right librarytype, check your files"
                    print (errorlog,file=ERROR)
                    error_count += 1
        ERROR.close()
        return(error_count)



    def Checkemergency(self):
        """
            从样本批次表里面获取加急信息。
            1. 如果样本里没有批次信息，报错
            2. 批次信息数目只能 >= 样品数目。如果不是，就需要把这个表格生成
        """
        error_count = 0
        Dic_PICI={}   ###保存批次表所有结果，在循环sample表的时候不断删除，最后得到批次表中有，sample表中没有的结果
        Dic_PICI_pass = {} ###保存sample表中有，batch表中也有的最终结果
        Dic_PICI_NOTIN_sample = {} ###保存sample表中有，而批次表中没有的结果，这个会报错

        interval_code = self.Checkencoding(self.interval)
        pici_code = self.Checkencoding(self.batch)
        INTERVAL = codecs.open(self.interval,"r",interval_code) 
        PICI = codecs.open(self.batch,"r",pici_code)
        ERROR = codecs.open(self.errorlog,"a","utf-8")

        PICI_title = self.linestrip(PICI.readline())
        INTERVAL_title = self.linestrip(INTERVAL.readline())
        for line in PICI:
            line = self.linestrip(line)
            linetemp = line.split(",")
            if Dic_PICI.has_key(linetemp[0]):
                errorlog = linetemp[0] + " is replicated in batch table!!!check your files!!!"
                print (errorlog,file=ERROR)
                error_count += 1
            else:
                Dic_PICI[linetemp[0]] = line
                

        ###开始循环interval表
        for line in INTERVAL:
            line = self.linestrip(line)
            linetemp = line.split(",")
            if linetemp[7] == "":
                pattern = linetemp[6]
            else:
                pattern = "-".join(linetemp[6:8])

            if Dic_PICI.has_key(pattern):
                linetemp[-1] = "" ##有时候拆分表会备注一些奇奇怪怪的东西。
                linetemp[-1] = (Dic_PICI[pattern].split(","))[-1]
                Dic_PICI_pass[pattern] = ",".join(linetemp)
                Dic_PICI.pop(pattern)
            else:
                Dic_PICI_NOTIN_sample[pattern] = line

        ###开始判断结果
        if len(Dic_PICI_NOTIN_sample) > 0:
            errorlog = "error:there are samples not in the batch table!!!"
            error_count += 1
            print (errorlog,file=ERROR)
            for k,v in Dic_PICI_NOTIN_sample.items():
                errorlog = k + " not in the batch table"
                print (errorlog,file=ERROR)

        if len(Dic_PICI) > 0:
            errorlog = "warning:there are samples in the batch table not sequenced!!!"
            print (errorlog,file=ERROR)
            for k in Dic_PICI.keys():
                errorlog = k + " not sequenced"
                print (errorlog,file=ERROR)

        INTERVAL.close()
        PICI.close()
        ERROR.close()

        ###输出结果
        with codecs.open(self.interval_add_emergency,"w","utf-8") as OUT:
            print (INTERVAL_title,file=OUT)
            for k in sorted(Dic_PICI_pass.keys()):
                print (Dic_PICI_pass[k],file=OUT)

        return(error_count)

    def GenerateSplitTable(self):
        out_title = ["sample_name","analysis_type","emergency","department","M reads","index","length"]
        IN = codecs.open(self.interval_add_emergency,"r","utf-8")
        OUT = codecs.open(self.split,"w","utf-8")
        IN_title = IN.readline()
        print(",".join(out_title),file=OUT)
        for line in IN:
            line = self.linestrip(line)
            linetemp = line.split(",")
            if linetemp[7] == "":
                sample_name = linetemp[6]
            else:
                sample_name = "-".join(linetemp[6:8])

            reformat = [sample_name,self.analysistype[linetemp[3]],linetemp[-1],"MED",linetemp[-4],linetemp[-3],linetemp[-2]]
            print(",".join(reformat),file=OUT)

        IN.close()
        OUT.close()

    def CheckSqlLIMS(self):
        """
            这个模块主要执行3个功能
            1. 检查WQ的命名问题，必须包含F/M/S
            2. 检查亲子样品是否都有case
            3. 检查产筛样品是否都有case
        """
        DicWQ = {}
        DicWS = {}
        WQ_sample_name = []
        WQ_error_count = 0
        with codecs.open(self.split,"r","utf-8") as IN:
        #with codecs.open(self.interval_add_emergency,"r","utf-8") as IN:
            for line in IN:
                line = self.linestrip(line)
                linetemp = line.split(",")
                if linetemp[0].startswith("WQ"):
                    if re.search(r'-S|-M|-F',linetemp[0]) == None: ## -s -f -m
                        WQ_error_count += 1
                    case_name = linetemp[0].split("-")[0]
                    DicWQ[case_name] = 1
                    WQ_sample_name.append(linetemp[0])
                if linetemp[0].startswith("WS"):
                    DicWS[linetemp[0]] = 1

        ###链接数据库
        ERROR = codecs.open(self.errorlog,"a","utf-8")
        con = mdb.connect(self.dbconfigure["host"], self.dbconfigure["name"], self.dbconfigure["passwd"], self.dbconfigure["database"], charset = 'utf8')
        cur = con.cursor()

        if len(DicWQ) > 0:
            WQ_sample = "','".join(DicWQ.keys())
            sql = "SELECT i.note, o.create_date, o.name, i.patient_name, st.name, i.code, i.parented, o.sjrq, " \
                    "i.accept_date, o.isurgent, o.gestational_weeks, t.name, t.principal FROM sample_info i " \
                    "LEFT JOIN sample_order o ON i.sample_order = o.id LEFT JOIN primary_task t ON t.id = o.advance " \
                    "LEFT JOIN dic_sample_type st ON st.id = i.sample_kind WHERE o.id IS NOT NULL AND i.note LIKE 'WQ%'" \
                    "and i.note in ('{}') ORDER BY o.create_date DESC".format(WQ_sample)
            cur.execute(sql)
            case_infor = cur.fetchall()
            ###循环核对样品名
            Dic_case = {}
            circle_error = 0
            for each in case_infor:
                each = [str(i).strip() for i in each]
                eachtemp = "{}-{}".format(each[0],each[6])
                Dic_case[eachtemp] = 1
                ###做一些特例,这部分的原则是，尽量让LIMS录入改好。以后碰到什么问题再修改对应的规则
                if not each[6]:        # 解决lims系统sample_info表parented字段为Null，文件检查无法通过
                    continue
                if len(each[6]) == 1:  ##如果长度只有1，只会出现两种情况，M或者F，把F1和F-1或者M1和M-1加入考虑
                    eachtemp1 = "{}-{}1".format(each[0],each[6])
                    eachtemp2 = "{}-{}-1".format(each[0],each[6])
                    Dic_case[eachtemp1] = 1
                    Dic_case[eachtemp2] = 1
                if len(each[6]) == 2:  ##如果长度是2的话，那么一般情况下是F1和M1，那么把F-1和M-1加入考虑
                    eachtemp3 = "{}-{}-{}".format(each[0],each[6][0],each[6][1])
                    Dic_case[eachtemp3] = 1

            for each in WQ_sample_name:
                if re.search(r'-S',each) == None:
                    if re.search(r'-T',each) != None:
                        each = each.split("-T")[0]
                    if not Dic_case.has_key(each):
                        errorlog = each + " do not found in LIMS,check your files"
                        print(errorlog, file=ERROR)
                        circle_error += 1

            ERROR.close()
            return circle_error
        else:      # add by hd 20180118 line 493-494 修复没有亲子样本时的bug
            return 0

    def CheckMongoDupName(self):
        error_count = 0
        db_pt = Config().get_mongo_client(mtype="pt")[Config().get_mongo_dbname("pt")]
        db_nipt = Config().get_mongo_client(mtype="nipt_v2")[Config().get_mongo_dbname("nipt_v2")]

        ERROR = codecs.open(self.errorlog,"a","utf-8")
        with codecs.open(self.split,"r","utf-8") as IN:
            title = IN.readline()
            for line in IN:
                line = self.linestrip(line)
                linetemp = line.split(',')
                name = linetemp[0]
                if name.startswith('WQ'): ##修改原则，以后如果需要检查重命名，在这里修改。
                    name_return = db_pt['sg_pt_qc'].find_one({'sample_id':name})
                elif name.startswith('WS'):
                    name_return = db_nipt['sg_sample_qc'].find_one({'sample_id':name})
                else:
                    continue

                if name_return != None:
                    errorlog = name + " is already used,check your files and change name"
                    print (errorlog,file=ERROR)
                    error_count += 1
        return(error_count)
        ERROR.close()

    def dump_error(self, status):
        with codecs.open(self.errorlog,"r","utf-8") as ERROR:
            errorlog = ERROR.readlines()
            errorlog = "".join(errorlog)

        # infor1 = "实验批次表:" + self.batch
        # infor2 = "线上上机表(LIMS系统生成):" + self.up
        # infor3 = "线下上机表(返工/测试样品):" + self.down
        # infor_table = [infor1,infor2,infor3]
        # infor_table = "\n".join(infor_table)

        db = Config().get_mongo_client(mtype="pt_v2")[Config().get_mongo_dbname("pt_v2")]

        # insert_data = {
        #     'created_ts':str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
        #     'status':status,
        #     'desc':"",
        #     'log':errorlog,
        #     'member_id':'4545',
        #     'info_table': infor_table,
        # }
        insert_data = {
            "log": errorlog
        }
        db['sg_file_check'].update({"_id": ObjectId(self.option("main_id"))}, {"$set": insert_data}, upsert=True)
        if status == 'failed':
            self.set_error("error")  # 这里error不能变，是用于后面更新的
            raise Exception("文件检查中发现错误，流程终止！")
        # collection = db['sg_file_check']
        # MAINID = collection.insert_one(insert_data).inserted_id
        print ("add sg_file_check!")
        #print (MAINID)
        # return MAINID

    def dump_guest(self):
        DicWQ = {}
        DicWS = {}
        WQ_sample_name = []
        WQ_error_count = 0
        with codecs.open(self.split,"r","utf-8") as IN:
            for line in IN:
                line = self.linestrip(line)
                linetemp = line.split(",")
                if linetemp[0].startswith("WQ"):
                    if re.search(r'-S|-M|-F',linetemp[0]) == None: ## -s -f -m
                        WQ_error_count += 1
                    case_name = linetemp[0].split("-")[0]
                    DicWQ[case_name] = 1
                    WQ_sample_name.append(linetemp[0])
                if linetemp[0].startswith("WS"):
                    DicWS[linetemp[0].split("-")[0]] = 1    # 产筛样品在lims系统客户信息中内部订单编号没有-1 -2 等后缀

        ###链接数据库
        con = mdb.connect(self.dbconfigure["host"], self.dbconfigure["name"], self.dbconfigure["passwd"], self.dbconfigure["database"], charset = 'utf8')
        cur = con.cursor()

        WQ_sample = "','".join(DicWQ.keys())
        WS_sample = "','".join(DicWS.keys())
        sql = "SELECT i.note, o.create_date, o.name, i.patient_name, st.name, i.code, i.parented, o.sjrq, " \
                "i.accept_date, o.isurgent, o.gestational_weeks, t.name, t.principal FROM sample_info i " \
                "LEFT JOIN sample_order o ON i.sample_order = o.id LEFT JOIN primary_task t ON t.id = o.advance " \
                "LEFT JOIN dic_sample_type st ON st.id = i.sample_kind WHERE o.id IS NOT NULL AND i.note LIKE 'WQ%'" \
                "and i.note in ('{}') ORDER BY o.create_date DESC".format(WQ_sample)
        #WQ17093590 2017-09-20 14:25:37 张瑶瑶 黄小明 精斑  JB1709200003    F1  2017-09-16 14:30:23 2017-09-20 14:24:57 0   26  曾蓉（网络客服）    曾蓉  示例
        cur.execute(sql)
        case_infor = cur.fetchall()

        db = Config().get_mongo_client(mtype="pt_v2")[Config().get_mongo_dbname("pt_v2")]
        collection = db['sg_guest_info']

        for each in case_infor:
            each = [str(i).strip() for i in each]
            if not each[6]:           # 解决：lims系统sample_info表parented字段为Null，文件检查无法通过
                continue
            if each[6] == 'M':        # add by hongdong @20180102 修改客户信息表中M就是M-1，在api中依旧进行下穷举
                sample = 'M-1'
                name = each[0] + "-M-1"
                sample_id = each[0] + "-M-1"
            else:
                sample = each[6]
                name = each[0] + "-" + each[6]
                sample_id = each[0] + "-" + each[6]
            insert_data = {
              "case_name": each[0],
              "create_ts": str(each[1]).split(" ")[0],
              "ask_person": each[2],
              "sample_name": each[3],
              "sample_type": each[4],
              "sample_number": each[5],
              "sample": sample,    # modified by hongdong 20171229  规范化样本名字
              "sample_id": sample_id,
              "ask_time": str(each[7]).split(" ")[0],
              "accept_time": str(each[8]).split(" ")[0],
              "gestation_week": each[10],
              "company_name": each[11],
              "contacts_people": each[12],
              "name": name,
              'insert_time': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            }
            collection.update({"name":each[0] + "-" + each[6]},{'$set':insert_data},upsert=True)

        # #导入WS的信息表  modified by hd 20180108 修改limus获取信息中年龄与采样时间为空，以及孕产史格式问题 615-646
        sql = "SELECT i.order_num, i.note, i.patient_name, i.receipt_date, i.accept_date, o.age, o.pregnancy_num, " \
              "o.last_menstrual_period, o.gestational_weeks, o.ivf, o.phone, o.medical_institutions, " \
              "o.attending_doctor, i.sample_type2, o.syndrome21, t.name, t.principal, o.yield_num FROM sample_info i " \
              "LEFT JOIN sample_order o ON i.sample_order = o.id LEFT JOIN primary_task t ON t.id = o.advance LEFT " \
              "JOIN dic_sample_type st ON st.id = i.sample_kind WHERE o.id IS NOT NULL AND i.note LIKE 'WS%' " \
              "and i.note in('{}')".format(WS_sample)

        cur.execute(sql)
        case_infor = cur.fetchall()

        db = Config().get_mongo_client(mtype="nipt_v2")[Config().get_mongo_dbname("nipt_v2")]
        collection = db['sg_customer']

        for each in case_infor:
            each = [str(i).strip() for i in each]
            insert_data = {
                "status": "正常",  # modified by hd 20171229 英文改成中文
                "number": each[0],  # 订单编号
                "report_num": each[1],  # 订单内部编号
                "patient_name": each[2] if each[2] else '/',  # 检测人姓名
                "sample_date": str(each[3]).split(" ")[0] if each[3] else '/',  # 采样日期
                "accepted_date": str(each[4]).split(" ")[0] if each[4] else '/',  # 接受日期
                "age": each[5] if each[5] else '/',
                "gestation": "孕{}产{}".format(each[6], each[17]) if each[6] and each[17] else '/',   # 孕产史
                "final_period": str(each[7]).split(" ")[0] if each[7] else '/',  # 末次月经
                "gestation_week": each[8] if each[8] else '/',  # 孕周
                "IVFET": "否" if each[9] == '0' else "是",   # modified by hd 20180108  IVF-ET妊娠
                "tel": each[10] if each[10] else each[10],
                "hospital": each[11] if each[11] else "/",
                "doctor": each[12] if each[12] else "/",
                "sample_type": each[13] if each[13] else "/",   # 样本类型
                "sample_id": each[1],
                "pregnancy": "单胎" , # modified by hd 20171229 英文改成中文
                "company_name": each[15],
                "contacts_people": each[16]
            }
            collection.update({"name": each[1]}, {'$set': insert_data}, upsert=True)

    def dump_experiment_batch(self):
        db = Config().get_mongo_client(mtype="pt_v2")[Config().get_mongo_dbname("pt_v2")]
        collection = db['sg_sample_batch']

        pici_code = self.Checkencoding(self.batch)
        with codecs.open(self.batch,"r",pici_code) as PICI:
            title = PICI.readline()
            for line in PICI:
                line = self.linestrip(line)
                linetemp = line.split(",")
                if linetemp[0].startswith('WQ') or linetemp[0].startswith('WS'):
                    insert_data = {
                        'sample_id':linetemp[0],
                        'extract_batch':linetemp[1],
                        'library_batch':linetemp[2],
                        'urgency_degree':linetemp[3],
                        'board_batch':self.flowcell,
                        'insert_time':datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    }
                    status = collection.update({"sample_id":linetemp[0]},{'$set':insert_data},upsert=True)

    def dump_sample_info(self):

        db = Config().get_mongo_client(mtype="pt_v2")[Config().get_mongo_dbname("pt_v2")]
        collection = db['sg_sample_info']

        with codecs.open(self.interval,"r","utf-8") as INTER:
            title = INTER.readline()
            for line in INTER:
                line = self.linestrip(line)
                linetemp = line.split(",")
                if linetemp[7] == "":
                    pattern = linetemp[6]
                else:
                    pattern = linetemp[6] + "-" + linetemp[7]
                    
                insert_table = {
                    'number': linetemp[0],
                    'sample_accept_time': linetemp[1],
                    'sample_name': linetemp[2],
                    'library_type': linetemp[3],
                    'analysis_type': self.analysistype[linetemp[3]],
                    'product_type': linetemp[4],
                    'sample_number': linetemp[5],
                    'case_name': linetemp[6],
                    'sample_id': pattern,   # modified by hd 20180111 统一sample_id 与余果姐后面家系组建衔接起来
                    'sample': linetemp[7],
                    'sample_type': linetemp[8],
                    'library_name': linetemp[9],
                    's_sequence_num': linetemp[10],
                    'index': linetemp[11],
                    'sequence_size': linetemp[12],
                    'board_batch': self.flowcell,
                    'pattern': pattern
                }
                collection.update({'sample_id': pattern}, {'$set': insert_table}, upsert=True)

    def dump_split(self):

        db = Config().get_mongo_client(mtype="pt_v2")[Config().get_mongo_dbname("pt_v2")]
        collection = db['sg_datasplit']

        WQ_sample = []
        WS_sample = []
        # split_type_count = 0
        with codecs.open(self.split,"r","utf-8") as IN:
            title = IN.readline()
            for line in IN:
                line = self.linestrip(line)
                linetemp = line.split(",")
                if linetemp[0].startswith("WS"):
                    WS_sample.append(linetemp[0])
                elif linetemp[0].startswith("WQ"):
                    WQ_sample.append(linetemp[0])
                    # split_type_count += 1
                # else:
                #     split_type_count += 1

        WQ_num = len(WQ_sample)
        WS_num = len(WS_sample)
        # if split_type_count > 0:
        #     split_type = "PE"
        # else:
        #     split_type = "SE"

        split_abs = os.getcwd() + "/" + self.split
        insert_table = {
            'board_batch':self.flowcell,
            "created_ts":datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),  # modified by hongdong
            #  20171230 create_time -> created_ts
            "wq_sample":WQ_sample,
            "ws_sample":WS_sample,
            "wq_num":WQ_num,
            "ws_num":WS_num,
            "start_time":"",
            "status":"unproduced",
            "desc":"",
            "end_time":"",
            'split_type': '',
            "split_info_path":split_abs,
            'board_path': self.get_board_path(self.flowcell)
        }
        collection.update({'board_batch':self.flowcell},{'$set':insert_table},upsert = True)

    def get_board_path(self, board_batch):
        """
        用于添加具体的板子所在路径，目前是nextseq 与nextseq1  add by hongdong 20180108
        :param board_batch:
        :return:
        """
        if os.path.exists("/mnt/clustre/upload/nextseq1/" + board_batch):
            board_path = "/mnt/clustre/upload/nextseq1/" + board_batch
        elif os.path.exists("/mnt/clustre/upload/nextseq/" + board_batch):
            board_path = "/mnt/clustre/upload/nextseq/" + board_batch
        else:
            board_path = ''
            self.set_error("{}不在nextseq与nextseq1中".format(board_batch))
        self.logger.info("board_path{}".format(board_path))
        return board_path

    def dump_sg_customer(self,MAINID):
        PICI_dic = {}
        SG = codecs.open(self.customer,"w","utf-8")
        pt_db = Config().get_mongo_client(mtype="pt_v2")[Config().get_mongo_dbname("pt_v2")]
        nipt_db = Config().get_mongo_client(mtype="nipt_v2")[Config().get_mongo_dbname("nipt_v2")]
        pt_collection = pt_db['sg_file_info']
        nipt_collection = nipt_db['sg_customer']
        code_type = self.Checkencoding(self.batch)
        nipt_count = 1    # modified by hongdong 20171229 747-749 删除每行后的； 因为python 语法不对
        pt_count = 1
        count = 0
        
        with codecs.open(self.batch,"r",code_type) as PICI:
            title = PICI.readline()
            for line in PICI:
                line = self.linestrip(line)
                linetemp = line.split(",")
                PICI_dic[linetemp[0]] = line
        
        with codecs.open(self.interval,"r","utf-8") as INTER:
            title = INTER.readline()
            for line in INTER:
                line = self.linestrip(line)
                linetemp = line.split(",")
                if linetemp[7] == "":
                    sample_id = linetemp[6]
                else:
                    sample_id = "-".join(linetemp[6:8])
                ##重构输出数据
                PICI_result = PICI_dic[sample_id].split(",")
                #tempout = [[sample_id,PICI_result[1:3],linetemp[3],linetemp[9:14]]]
                tempout = [sample_id]
                tempout.extend(PICI_result[1:3])
                tempout.append(linetemp[3])
                tempout.extend(linetemp[9:14])
                tempout.append(PICI_result[3])

                if sample_id.startswith("WQ"):
                    count = pt_count
                    pt_count += 1
                elif sample_id.startswith("WS"):
                    count = nipt_count
                    nipt_count += 1

                insert_table = {
                    "number" : count,
                    "sample_id" : tempout[0],
                    "extract_batch" : tempout[1],
                    "library_batch" : tempout[2],
                    "library_type" : tempout[3],
                    "library_name" : tempout[4],
                    "index" : tempout[6],
                    "s_sequence_num" : tempout[5],
                    "sequence_size" : tempout[7],
                    "urgence" : tempout[9],
                    "use_type" : tempout[8],
                    "upload_time" : datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    "board_batch" : self.flowcell,
                    "check_id" : ObjectId(MAINID)
                }
                
                if sample_id.startswith("WQ"):
                    pt_collection.update({'sample_id': tempout[0]}, {'$set': insert_table}, upsert=True)
                elif sample_id.startswith("WS"):
                    pt_collection.update({'sample_id': tempout[0]}, {'$set': insert_table}, upsert=True)
                    # nipt_collection.update({'sample_id':tempout[0]},{'$set':insert_table},upsert = True)

                tempout = ",".join(list(tempout))
                print(tempout,file=SG)
        SG.close()
