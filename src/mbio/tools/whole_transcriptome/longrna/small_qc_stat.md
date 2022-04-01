质控数据统计
===========


-----------------------------------

PATH
-----------------------------------
/质控/质控-开发文档.md

-----------------------------------
## 参数设置
+ 展示样式
<br>`参数设置`
<br>不需要设置

## 质控数据统计


### 1、原始数据统计表
###    质控数据统计表

+ 展示样式
<br>![原型图](images/qc_tab1.png)
<br>不需要筛选


+ 结果表数据：sg_specimen

| 原型表头字段        | 对应mongo表字段  | 条件       | 数据库表名  |
| :---               | :---            | :---       |   :---      |
| Quality Assessment |                 | 画图的链接  | sg_specimen |
| Sample Name        | new_name        |            | sg_specimen |
| Raw reads          | total_reads     |            | sg_specimen |
| Raw Bases (bp)     | total_bases     |            | sg_specimen |
| Raw Error Rate (%) | error_rate      |            | sg_specimen |
| Raw Q20 (%)        | q20_rate        |            | sg_specimen |
| Raw Q30 (%)        | q30_rate        |            | sg_specimen |
| GC content(%)      | qc_rate         |            | sg_specimen |

注， 仅筛选 about_qc==‘before’ 的记录


### 2、测序质量分布图

+ 展示样式

![碱基质量质量分布图](images/qc_raw_pic1.png)

+ 数据表及筛选 sg_specimen

	- 数据表：
	- 筛选字段：根据about_qc==‘before’ or about_qc==‘after’ 获取数据);type(left/right,成对出现并绘制在一张图上，用淡灰色虚线区分)
	- 页面筛选项：无 
	- 页面动态数据筛选项: 无

图的参数文件：

|    来源    | 页面名称 |                                                           值                                                           | Mongo字段 |   值范围   |                          默认值                           |
|------------|----------|------------------------------------------------------------------------------------------------------------------------|-----------|------------|-----------------------------------------------------------|
| 页面筛选项 | 颜色设置 | \<input type=”text” name=”color” value=”#A020F0”/>                                                                     |           |            | #A020F0                                                   |
| 页面筛选项 | 主标题   | \<input type=”text” name=”title” value=”Base quality distribution of raw data summarized by reads positon：样品名称”/> |           |            | "Base quality distribution of raw data summarized by reads positon：" + 样品名称 |
| 页面筛选项 | 主标题   | \<input type=”checkbox” name=”title”/>                                                                                 |           | TRUE/FALSE | TRUE                                                      |
| 页面筛选项 | X轴标题  | \<input type=”text” name=”xlimtitle” value=”reads position”/>                                                           |           |            | "reads position"                                           |
| 页面筛选项 | Y轴标题  | \<input type=”text” name=”ylimtitle” value=”Quality Score”/>                                                           |           |            | "Quality Score"                                           |

```
注：1. 颜色设置根据固定的调色板选择，默认使用原型中的颜色；
    2. 横坐标标签步长按照图片大小自动适应。
```


+ 该图所需参数

|        图的key        |         Value          | 是否可调 | 是否从数据库中取出 |                                                           说明                                                          |
|-----------------------|------------------------|----------|--------------------|-------------------------------------------------------------------------------------------------------------------------|
| data                  | [max,min,q1,q3,median] | 否       | 是                 | 从sg_specimen_graphic表取出的碱基质量信息                                                                                 |
| categories            | column                 | 是       | 是                 | 碱基位点                                                                                                                |
| size[width]           | 400                    | 是       | 否                 | 图宽                                                                                                                    |
| size[height]          | 300                    | 是       | 否                 | 图高                                                                                                                    |
| params[rotation]      | 20                     | 否       | 否                 |                                                                                                                         |
| params[x_label]       | "read position"        | 是       | 否                 | 页面文本框传入                                                                                                          |
| params[y_label]       | "Quality Score"        | 是       | 否                 | 页面文本框传入                                                                                                          |
| params[title]         | "Base quality distribution of raw data summarized by reads positon：样品名称"                     | 是       | 否                 | 页面文本框传入，文本框为空值时此参数值为"Base quality distribution of raw data summarized by reads positon:" + 样品名称 |
| params[stack]         | 1                      | 否       | 否                 |                                                                                                                         |
| params[colors]        | "#A020F0"              | 是       | 否                 | 页面调色工具传入，默认#"A020F0"                                                                                         |
| params[tooltip_names] | "all"                  |          | 否                 |                                                                                                                         |
| params[show_legend]   | true                   | 否       | 否                 | 是否显示legend                                                                                                          |



### 3、测序错误率分布图

+ 展示样式

![碱基错误率分布图](images/qc_raw_pic2.png)

+ 数据表及筛选 

	- 数据表：datastat_graphic
	- 筛选字段：根据about_qc==‘before’ or about_qc==‘after’ 获取数据);type(left/right,成对出现并绘制在一张图上，用淡灰色虚线区分)
	- 页面筛选项：无 
	- 页面动态数据筛选项: 无

图的参数文件：

|    来源    | 页面名称 |                                                           值                                                           | Mongo字段 |   值范围   |                          默认值                           |
|------------|----------|------------------------------------------------------------------------------------------------------------------------|-----------|------------|-----------------------------------------------------------|
| 页面筛选项 | 颜色设置 | \<input type=”text” name=”color” value=”#A020F0”/>                                                                     |           |            | #A020F0                                                   |
| 页面筛选项 | 主标题   | \<input type=”text” name=”title” value=”Base quality distribution of raw data summarized by reads positon：样品名称”/> |           |            | "Base quality distribution of raw data summarized by reads positon：" + 样品名称 |
| 页面筛选项 | 主标题   | \<input type=”checkbox” name=”title”/>                                                                                 |           | TRUE/FALSE | TRUE                                                      |
| 页面筛选项 | X轴标题  | \<input type=”text” name=”xlimtitle” value=”reads position”/>                                                           |           |            | "reads position"                                           |
| 页面筛选项 | Y轴标题  | \<input type=”text” name=”ylimtitle” value=”Quality Score”/>                                                           |           |            | "Quality Score"                                           |

```
注：1. 颜色设置根据固定的调色板选择，默认使用原型中的颜色；
    2. 横坐标标签步长按照图片大小自动适应。
```


+ 该图所需参数

| 图的key               | Value                                                                         | 是否可调 | 是否从数据库中取出 | 说明                                                                                                                    |
|-----------------------+-------------------------------------------------------------------------------+----------+--------------------+-------------------------------------------------------------------------------------------------------------------------|
| data                  | error                                                                         | 否       | 是                 | 从sg_specimen_graphic表取出的碱基质量信息                                                                               |
| categories            | column                                                                        | 是       | 是                 | 碱基位点                                                                                                                |
| size[width]           | 400                                                                           | 是       | 否                 | 图宽                                                                                                                    |
| size[height]          | 300                                                                           | 是       | 否                 | 图高                                                                                                                    |
| params[rotation]      | 20                                                                            | 否       | 否                 |                                                                                                                         |
| params[x_label]       | "read position(bp)"                                                           | 是       | 否                 | 页面文本框传入                                                                                                          |
| params[y_label]       | "Error rate %"                                                                | 是       | 否                 | 页面文本框传入                                                                                                          |
| params[title]         | "Base quality distribution of raw data summarized by reads positon：样品名称" | 是       | 否                 | 页面文本框传入，文本框为空值时此参数值为"Base quality distribution of raw data summarized by reads positon:" + 样品名称 |
| params[stack]         | 1                                                                             | 否       | 否                 |                                                                                                                         |
| params[colors]        | "#A020F0"                                                                     | 是       | 否                 | 页面调色工具传入，默认#"A020F0"                                                                                         |
| params[tooltip_names] | "all"                                                                         |          | 否                 |                                                                                                                         |
| params[show_legend]   | true                                                                          | 否       | 否                 | 是否显示legend                                                                                                          |


### 4、质控数据统计表

+ 展示样式
<br>
![原型图](images/qc_clean_tab1.png)
<br>不需要筛选


+ 结果表数据：

| 原型表头字段            | 对应mongo表字段 | 条件       | 数据库表名  |
|-------------------------+-----------------+------------+-------------|
| Quality Assessment      |                 | 画图的链接 | sg_specimen |
| Sample Name             | new_name        |            | sg_specimen |
| Raw reads               | total_reads     |            | sg_specimen |
| Raw Bases (bp)          | total_bases     |            | sg_specimen |
| Raw Error Rate (%)      | error_rate      |            | sg_specimen |
| Raw Q20 (%)             | q20_rate        |            | sg_specimen |
| Raw Q30 (%)             | q30_rate        |            | sg_specimen |
| GC content(%)           | gc_rate         |            | sg_specimen |
| Useful reads(18nt-32nt) | useful_reads    |            | sg_specimen |

注， 仅筛选 about_qc==‘after’ 的记录

### 5、分布图

+ 展示样式

![质控分布图](images/qc_clean_pic1.png)

+ 数据表及筛选 

	- 数据表：sg_qclen sg_qclen_detail
	- 筛选字段：根据about_qc==‘before’ or about_qc==‘after’ 获取数据);type(left/right,成对出现并绘制在一张图上，用淡灰色虚线区分)
	- 页面筛选项：无 
	- 页面动态数据筛选项: 无

| 来源       | 页面名称 | 值  | Mongo表和字段   |
| ---        | ---      | --- | ---             |
| 页面筛选项 | Sample   |     | sg_qclen_detail |



| 图的key             | Value                         | 是否可调 | 是否数据库中取出 | 说明     |
| ---                 | ---                           | ---      | ---              | ---      |
| data                |                               | 是       | 是               | len_data |
| color               |                               | 是       |                  |          |
| params[text]        | Sequence lengths Distribution | 是       |                  |          |
| params[x_labe]      | Length                        | 是       |                  |          |
| params[y_label]     | Num                           | 是       |                  |          |
| params[legend]      |                               | 是       |                  |          |
| params[show_legend] | TRUE                          |          |                  |          |
| size[width]         | 600                           |          |                  |          |
| size[height]        | 800                           |          |                  |          |
| Show_legend         | TRUE                          |          |                  |          |

展示样式：
![原型图](images/qc_clean_tab2.png)
数据表 sg_qclen_detail


| 原型表头字段 | 对应mongo表字段  | 数据库表名      | 条件 |
|--------------+------------------+-----------------+------|
| length       | len_data 的key   | sg_qclen_detail |      |
| number       | len_data 的value | sg_qclen_detail |      |

