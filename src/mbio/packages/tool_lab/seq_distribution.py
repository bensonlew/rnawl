

"""
对fa文件中的序列进行处理
# 获取序列的id和序列信息
# 统计每个id对应的序列的长度
# 对序列长度进行统计
# 绘制直方图
"""
import os
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import sys
def read_seq(file_path):
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                seq_id.append(line.rstrip('\n').replace('>', ''))
            else:
                seq.append(line.rstrip('\n'))
    seq_dic = dict(zip(seq_id,seq))  # 将序列编号和序列信息两个列表合并成字典
    for v in seq_dic.values():
        seq_len.append(len(v))
    return seq_id , seq, seq_len


def write_count(seq_id, seq, seq_len):
    # pandas 写入序列编号、序列信息、序列长度
    data1 = pd.DataFrame({"Seq_ID": seq_id})
    data2 = pd.DataFrame({"Seq_Info": seq})
    data3 = pd.DataFrame({"Seq_Len": seq_len})

    writer = pd.ExcelWriter(abs_path + '\\' + "test1.xlsx") # windows 下使用
    data1.to_excel(writer,sheet_name="data",startcol=0,index=False)
    data2.to_excel(writer,sheet_name="data",startcol=1,index=False)
    data3.to_excel(writer,sheet_name="data",startcol=2,index=False)
    #writer.save()  # 数据保存为excel文件

def count_bar(seq_len):
    """
    # 根据上一步获得的序列长度信息，对其进行sort/uniq，matplotlib 处理并绘制直方图
    # 首先对数据进行排序统计对相同长度进行计数
    # 数据清洗后进行画bar图
    """
    len_count = Counter(seq_len)  # 提取上一步获得第三列长度数据进行清洗并统计每个长度的个数
    # matplotlib绘图
    x = []
    y = []
    for k ,v in len_count.items():
        x.append(k)
        y.append(v)
    print(x ,y)
    plt.bar(x,y)
    plt.xlabel("Sequence Length")
    plt.ylabel("Sequence Numbers")
    plt.show()

#def main():


if __name__ == "__main__":
    abs_path = os.getcwd()  # 获取当前目录路径
    print(abs_path)
    file_name = sys.argv[0]
    file_path = abs_path + '\\' + '*.fa'  # 获取当前目录下的文件信息
    seq_id = []  # 新建列表存储fasta文件中序列编号信息
    seq = []  # 新建列表存储对应fasta文件中序列信息
    seq_len = []  # 新建列表存储对应序列的长度信息
    read_seq(file_path)
    write_count(seq_id, seq, seq_len)
    count_bar(seq_len)
