# -*_ coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171115
from biocluster.core.exceptions import FileError
from biocluster.iofile import File
import os
import json
import MySQLdb
import xml.etree.ElementTree as ET


class SampleInfoFile(File):
    """
    定义高通量数据拆分的输入json文件的格式
    """
    def __init__(self):
        super(SampleInfoFile, self).__init__()

    def dump_json(self):
        """解析json文件"""
        f = open(self.path, 'rb')
        try:
            json_dict = json.loads(f.read())
        except:
            raise FileError('json格式不正确')
        return json_dict

    def get_info(self):
        """
        获取文件属性
        """
        super(SampleInfoFile, self).get_info()
        self.json_dict = self.dump_json()
        self.set_property('data_path', self.json_dict['data_path'])
        self.set_property('lane_info', self.json_dict['lane_info'])
        self.set_property('library_info', self.json_dict['library_info'])
        self.set_property('sample_info', self.json_dict['sample_info'])
        seq = self.create_base_mask()
        if self.json_dict['seq_mode'] != seq[1]:
            raise FileError('json里的测序模式和RunInfo.xml里的测序模式不同，请检查json里的index序列、测序模式')
        self.set_property('seq_type', seq[0])
        self.set_property('seq_mode', seq[1])

    def check(self):
        """
        检测文件是否满足要求，发生错误时应该触发FileError异常
        """
        super(SampleInfoFile, self).check()
        self.json_dict = self.dump_json()
        prop = ['board_name', 'lane_info', 'library_info', 'sample_info', 'data_path', 'seq_mode']
        for p in prop:
            if p not in self.json_dict:
                raise FileError('json中缺少属性：' + p)
        if not os.path.exists(self.json_dict['data_path']):
            raise FileError('json中data_path文件不存在')
        self.bcl_path = os.path.join(self.json_dict['data_path'], 'Data/Intensities/BaseCalls/')
        if not os.path.exists(self.bcl_path):
            raise FileError("bcl文件夹：{}缺失，请检查".format(self.bcl_path))
        self.run_info = os.path.join(self.json_dict['data_path'], 'RunInfo.xml')
        if not os.path.exists(self.run_info):
            raise FileError("RunInfo.xml文件:{}缺失，请检查".format(self.run_info))
        lane_info = ['lane_name', 'lane_data_size', 'pre_lane_data_size']
        for lane in self.json_dict['lane_info']:
            for p in lane_info:
                if p not in lane.keys():
                    raise FileError('json中缺少lane的属性：' + p)
        library_info = ['library_number', 'lane_name', 'library_type', 'index', 'index_seq', 'index2', 'index2_seq', 'convert_data_size', 'order_data_size', 'setup_seq_size', 'library_insert', 'library_denisity', 'has_child']
        for library in self.json_dict['library_info']:
            for p in library_info:
                if p not in library.keys():
                    raise FileError('json中缺少library的属性：' + p)
        sample_info = ['seq_type', 'library_type', 'ticket_number', 'project_sn', 'client_name', 'library_number', 'sample_name', 'p1_barcode', 'primer', 'barcode']
        for s in self.json_dict['sample_info']:
            for p in sample_info:
                if p not in s.keys():
                    raise FileError('json中缺少sample的属性：' + s)
        self.indexs, self.index2s = [], []
        for library in self.json_dict['library_info']:
            if library['index_seq'] == '':
                raise FileError('文库的index缺失，请检查')
            if library['index_seq'] not in self.indexs:
                self.indexs.append(library['index_seq'])
            else:
                raise FileError('同一个板中文库的index：{}重复，请检查'.format(library['index_seq']))
            if library['index2_seq'] not in self.index2s:
                self.index2s.append(library['index2_seq'])
            elif library['index2_seq'] != "":
                raise FileError('同一个板中文库的index2：{}重复，请检查'.format(library['index2_seq']))
            else:
                pass

    def create_sample_sheet(self, seq_type, csv):
        csv = open(csv, 'w')
        if seq_type == 'PE':
            head = '[Data],,,,,,,,,,\nLane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well\
                    ,I7_Index_ID,index,I5_index_ID,index2,Sample_Project,Description\n'
            csv.write(head)
            for library in self.json_dict['library_info']:
                line = library['lane_name'] + ',' + library['library_number'] +',' +\
                       library['library_number'] + ',,,' + library['index2'] + ',' +\
                       library['index2_seq'] + ',' + library['index'] + ',' + library['index_seq'] + ',Fastq,\n'
                csv.write(line)
        if seq_type == 'SE':
            head = '[Data],,,,,,,,\nLane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,\
                    I7_Index_ID,index,Sample_Project,Description\n'
            csv.write(head)
            for library in self.json_dict['library_info']:
                line = library['lane_name'] + ',' + library['library_number'] +',' +\
                       library['library_number'] + ',,,' + library['index'] + ',' + library['index_seq'] + ',Fastq,\n'
                csv.write(line)

    def create_base_mask(self):
        """
        解析"RunInfo.xml"文件，获取机器模式中的index长度,生成测序模式
        """
        e = ET.parse(self.run_info).getroot()
        machineIndexLen = 0
        machineIndexLen2 = 0
        for ele in e.iter("Read"):
            if ele.attrib["Number"] == "1" and ele.attrib['IsIndexedRead'] == "N":
                sqMethod1 = int(ele.attrib['NumCycles'])
            if ele.attrib["Number"] == "4":  # 双端
                self.seq_type = "PE"
                if ele.attrib["Number"] == "4" and ele.attrib['IsIndexedRead'] == "N":
                    sqMethod2 = int(ele.attrib['NumCycles'])
                if ele.attrib["Number"] == "2" and ele.attrib['IsIndexedRead'] == "Y":
                    machineIndexLen = int(ele.attrib['NumCycles'])
                if ele.attrib["Number"] == "3" and ele.attrib['IsIndexedRead'] == "Y":
                    machineIndexLen2 = int(ele.attrib['NumCycles'])
            else:  #  单端
                self.seq_type = "SE"
                if ele.attrib["Number"] == "2" and ele.attrib['IsIndexedRead'] == "Y":
                    machineIndexLen = int(ele.attrib['NumCycles'])
                if ele.attrib["Number"] == "3" and ele.attrib['IsIndexedRead'] == "N":
                    sqMethod2 = int(ele.attrib['NumCycles'])
        if machineIndexLen == 0:
            raise FileError("RunInfo.xml解析出错, 无法解析上机时的index长度，请检查{}是否正确".format(self.run_info))
        if self.seq_type == "PE" and machineIndexLen2 == 0:
            raise FileError("RunInfo.xml解析出错, 无法解析上机时的index长度，请检查{}是否正确".format(self.run_info))
        minIndexLen = 9999  # 单端
        # myCode2Index = dict()  # 为了减少数据库的查询次数
        # for index_code in self.indexs:
        #     index = self.code2index(index_code)[0]
        #     myCode2Index[index_code] = index
        #     if len(index) > machineIndexLen:
        #         raise FileError("文库的index长度大于机器模式中的index长度")
        #     if len(index) < minIndexLen:
        #         minIndexLen = len(index)
        for index in self.indexs:
            if len(index) > machineIndexLen:
                raise FileError("文库的index长度大于机器模式中的index长度")
            if len(index) < minIndexLen:
                minIndexLen = len(index)
        iContent = str(minIndexLen) + "n" * int(machineIndexLen - minIndexLen)
        self.base_mask = "y{},i{},y{}".format(sqMethod1, iContent, sqMethod2)  # 单端
        if self.seq_type == "PE":  # 双端
            minIndexLen2 = 9999
            # for index_code2 in self.index2s:
            #     index2 = self.code2index(index_code2)[0]
            #     myCode2Index[index_code2] = index2
            #     if len(index2) > machineIndexLen2:
            #         raise FileError("文库的index2长度大于机器模式中的index2长度")
            #     if len(index2) < minIndexLen2:
            #         minIndexLen2 = len(index2)
            for index2 in self.index2s:
                if len(index2) > machineIndexLen2:
                    raise FileError("文库的index2长度大于机器模式中的index2长度")
                if len(index2) < minIndexLen2:
                    minIndexLen2 = len(index2)
            iContent1 = str(minIndexLen) + "n" * int(machineIndexLen - minIndexLen)
            iContent2 = str(minIndexLen2) + "n" * int(machineIndexLen2 - minIndexLen2)
            self.base_mask = "y{},i{},i{},y{}".format(sqMethod1, iContent1, iContent2, sqMethod2)
        return self.seq_type, self.base_mask

    # def code2index(self, code):
    #     """
    #     根据一个index的代码，查询mysql数据库，获取具体的index序列
    #     """
    #     db = MySQLdb.connect("192.168.10.51", "mydb", "mydb", "mjanalysis")
    #     cursor = db.cursor()
    #     sql = "SELECT * FROM sg_index_info WHERE label=\"{}\"".format(code)
    #     try:
    #         cursor.execute(sql)
    #         results = cursor.fetchall()
    #     except:
    #         raise Exception("无法从index数据库中获得信息")
    #     if len(results) == 0:
    #         raise ValueError("未找到该index代码: {}".format(code))
    #     left_index = results[0][3]
    #     right_index = results[0][4]
    #     varbase = results[0][2]
    #     if len(varbase) == 1:
    #         f_varbase = varbase
    #         r_varbase = varbase
    #     elif len(varbase) == 2:
    #         f_varbase = varbase[0]
    #         r_varbase = varbase[1]
    #     try:
    #         f_varbase = int(f_varbase)
    #         r_varbase = int(r_varbase)
    #     except:
    #         pass
    #     return (left_index, right_index, f_varbase, r_varbase)



if __name__ == '__main__':
    a = SampleInfoFile()
    a.set_path('/mnt/ilustre/users/sanger-dev/sg-users/zengjing/datasplit/new/meta_sample_info.json')
    a.check()
    a.get_info()
    # csv = '/mnt/ilustre/users/sanger-dev/sg-users/zengjing/datasplit/new/test1_sample.csv'
    # seq_type = 'SE'
    # a.create_sample_sheet(seq_type, csv)
    # s1 = a.create_base_mask()
    # print s1[0]
