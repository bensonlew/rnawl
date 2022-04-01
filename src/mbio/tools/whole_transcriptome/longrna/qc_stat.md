
> sg_specimen
### fadsf
# 我的尝试
> ##### **代码**
'''

import numpy as np

import pandas as pd

'''
### '左对齐 &&&&'
我的[](https://www.jianshu.com/p/3115d2260ec2)


'import collections'

====================
>Collection: **fdasf**  区块引用
*  注意： about_qc 用来是否区分质控， before 代表rawdata， after代表cleandata <br>
   该表表示四个页面展示<br>
   (fdasf)

>__1.1节__  or **给我变粗**

~~发生发生地方~~

###### ![](D:\wiki\lnc_rna.wiki\qc\biaoda.jpg)



| 列名 | 值 | 类型 | 描述 |
| :--- | :--- | :--- | :---|
|_id|ObjectId("5c2404e2v6g2fg16272367b3h")|ObjectId|None|
|created_ts|2019-02-19 16:43:26|string|创建时间|
|task_id|lnc_rna|String|任务id|
|old_name|A1|String|原始名称|
|new_name|A1|String|更改后名称|
|alias|S1|String|别名，用于传递原始名称和更改后名称的唯一标志|
|error_rate|0.0201 |Double |错误率 |
|gc_rate|0.62 |Double | GC含量|
|q20_rate|0.90 |Double |Q20 | 
|q30_rate|0.86 |Double |Q30 |
|a_rate|0.25 |Double |碱基A含量 | 
|t_rate|0.25 |Double | 碱基T含量|
|g_rate|0.25 |Double | 碱基G含量| 
|c_rate| 0.25| Double| 碱基C含量|
|n_rate| 0.001| Double| 碱基N含量|
|reads_with_ns|16354 |Int32 | 含有N的reads数|
|n_reads_rate|0.06 |Double |含有N的reads百分比|
|total_reads|1325666|Int32|reads总数|
|total_bases|165372223 |Int64 | 碱基总数|
|type|pe|String|测序类型，单端SE，双端PE|
|about_qc|before|String|原始数据before，质控数据after|
|group|A|String|组名|
|desc|肝脏|String|样品描述|
|rRNA_Ratio|1.631|Double|核糖体rna占比|

*  注：about_qc用来区分是否质控，before代表rawdata，after代表cleandata <br>
    该表对应四个页面展示:<br>

| 页面 | 字段 | 说明 |
| :--- | :--- | :---|
| 项目背景下的样品信息表 | alias, old_name, new_name, group, desc |  |
| 测序数据质控下的测序数据统计表| new_name, total_reads, total_bases, total_reads, total_bases, error_rate, q20_rate, q30_rate, gc_rate,rRNA_Ratio | 第一个total_reads和total_bases取about_qc=before条件下的值，其余的取about_qc=after条件下的值|
| 测序数据质控下的原始数据统计表| new_name, total_reads, total_bases, error_rate, q20_rate, q30_rate, gc_rate | 取about_qc=before条件下的值|
| 测序数据质控下的质控数据统计表| new_name, total_reads, total_bases, error_rate, q20_rate, q30_rate, gc_rate,rRNA_Ratio | 取about_qc=after条件下的值|


> Collection: **sg_specimen_graphic**

| 列名          |                                   值 | 类型     | 描述                                        |
| :---          |                                 :--- | :---     | :---                                        |
| _id           | ObjectId("5b695095a4e1af1846489340") | ObjectId | ?                                           |
| q1            |                                   32 | integer  |                               |
| A             |                     8.99184814119884 | float    | 碱基组成分布图                              |
| q3            |                                   32 | integer  |                           |
| G             |                     46.5237521648634 | float    | 碱基组成分布图                              |
| about_qc      |                               before | string   | before代表质控前的统计结果，after代表质控后 |
| min           |                                   32 | integer  |                             |
| column        |                                    1 | integer  | ?                                           |
| specimen_id   | ObjectId("5b695094a4e1af184648933a") | ObjectId | 与sg_specimen关联的字段                     |
| error         |                   0.0698232404077172 | float    | 碱基错误率分布图                            |
| median        |                                   32 | integer  |                             |
| C             |                     38.6280311321645 | float    | 碱基组成分布图                              |
| N             |                    0.358162431731681 | float    | 碱基组成分布图                              |
| T             |                     5.49820613004164 | float    | 碱基组成分布图                              |
| max           |                                   32 | integer  | |
| type          |                                 left | string   | 判断是哪端                                  |
| specimen_name |                                 WT_1 | string   | ?                                           |