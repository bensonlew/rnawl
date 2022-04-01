# -*- coding: utf-8 -*-
# __author__ = 'xuting'
# __version__ = 'v1.0'
# __last_modified__ = '20160303'
import re
import MySQLdb


def reverse_complement(string):
    """
    将一个序列进行互补
    :param string: 输入的序列
    """
    length = len(string)
    newstring = ""
    for i in range(length):
        if string[i] == "A":
            newstring = newstring + "T"
        elif string[i] == "a":
            newstring = newstring + "t"
        elif string[i] == "T":
            newstring = newstring + "A"
        elif string[i] == "t":
            newstring = newstring + "a"
        elif string[i] == "C":
            newstring = newstring + "G"
        elif string[i] == "G":
            newstring = newstring + "C"
        elif string[i] == "c":
            newstring = newstring + "G"
        elif string[i] == "g":
            newstring = newstring + "c"
        else:
            newstring = newstring + string[i]
    return newstring


def find_seq_len(sequencing_id):
    """
    根据一个测序版的id，获取测序长度， 用来生成base_mask
    """
    sequencing_id = int(sequencing_id)
    db = MySQLdb.connect("192.168.10.51", "mydb", "mydb", "mjanalysis")
    cursor = db.cursor()
    sql = "SELECT * FROM sg_sequencing WHERE sequencing_id =\"{}\"".format(sequencing_id)
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
    except:
        raise Exception("无法从测序板表sg_sequencing中获得信息")
    if len(results) == 0:
        raise ValueError("测序版id错误, 找不到id为{}的测序板".format(sequencing_id))
    sqMethod = results[0][11]
    return sqMethod


def code2index(code):
    """
    根据一个index的代码，获取具体的index序列
    """
    db = MySQLdb.connect("192.168.10.51", "mydb", "mydb", "mjanalysis")
    cursor = db.cursor()
    sql = "SELECT * FROM sg_index_info WHERE label=\"{}\"".format(code)
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
    except:
        raise Exception("无法从index数据库中获得信息")
    if len(results) == 0:
        raise ValueError("未找到该index代码: {}".format(code))
    left_index = results[0][3]
    right_index = results[0][4]
    varbase = results[0][2]
    if len(varbase) == 1:
        f_varbase = varbase
        r_varbase = varbase
    elif len(varbase) == 2:
        f_varbase = varbase[0]
        r_varbase = varbase[1]
    try:
        f_varbase = int(f_varbase)
        r_varbase = int(r_varbase)
    except:
        pass
    return (left_index, right_index, f_varbase, r_varbase)


def code2primer(code):
    """
    根据一个primer的代码，获取具体的primer
    """
    db = MySQLdb.connect("192.168.10.51", "mydb", "mydb", "mjanalysis")
    cursor = db.cursor()
    code = re.split('_', code)
    for my_code in code:
        sql = "SELECT * FROM sg_primer_info WHERE label=\"{}\"".format(my_code)
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
        except:
            raise Exception("无法从primer数据库中获得信息")
        if len(results) == 0:
            raise ValueError("未找到该primer代码")
        if re.search("F$", my_code):
            f_primer = results[0][2]
        if re.search("R$", my_code):
            r_primer = results[0][2]
    return (f_primer, r_primer)


def str_check(real_str, list_str):
    """
    比较两个index，返回两个字符串不同字符的个数
    """
    length = len(real_str)
    count = 0
    ABBR = {"A": "AA", 'C': "CC", 'T': "TT", 'G': "GG",
            'M': "AC", 'R': "AG", 'W': "AT", 'S': "CG",
            'Y': "CT", 'K': "GT", 'V': "ACG", 'H': "ACT",
            'D': "AGT", 'B': "CGT", 'X': "ACGT", 'N': "ACGT"}
    for i in range(length):
        realbase = real_str[i]
        if list_str[i] in ABBR.keys():
            indexbase = ABBR[list_str[i]]
        else:
            indexbase = list_str[i]
        if not re.search(realbase, indexbase):
            count = count + 1
    return count

if __name__ == "__main__":
    real_str = "GTGCCAGCCGCCGCGG"
    list_str = "GTGCCAGCMGCCGCGG"
    test = str_check(real_str, list_str)
    test2 = code2index("NB67", "miseq")
    test3 = code2primer("1737F_2043R")
    print test
    print test2
    print test3
