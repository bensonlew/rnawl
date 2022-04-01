# -*- coding: utf-8 -*-
from lxml import etree
import base64
import re
from biocluster.config import Config

def get_mongo():
    client = Config().get_mongo_client(mtype="meta")
    db = client[Config().get_mongo_dbname("meta")]
    return client, db

def save_params(out_file, run_id, params=None):
    mongo, db = get_mongo()
    sg_params = db["sg_params"]
    task_id = '_'.join(run_id.split('_')[:2])
    if not params:
        sg_status = db["sg_status"]
        run_info = sg_status.find_one({"task_id": task_id, "run_id": run_id})
        params = run_info["params"]
    params_info = sg_params.find_one({"task_id": task_id, "params": params})
    if not params_info:
        pass
        #raise Exception("缺失参数信息")
    else:
        parse_string(base64.decodestring(params_info["html"]), out_file+"/运行参数.txt")
    # sg_params.delete_one({"task_id": task_id, "params": run_info["params"]})


def parse_string(html_str, out_file):
    html = etree.HTML(html_str.decode('utf-8'))
    all_opts = html.xpath("//div[@class='wrap']/form/ul") or html.xpath("//div[@class='wrap']/ul")
    all_opts2 = html.xpath("//ul[@class='edit-data display-inline']/form/ul") or html.xpath("//ul[@class='edit-data display-inline']/ul")
    with open(out_file, "wb") as w:
        index_list =[]
        for one_line in all_opts:
            opts = one_line.xpath("li")
            for one_opt in opts:
                name = ''
                name1 = ''
                if one_opt.get("class") and "select-group-parent" in one_opt.get("class"):
                    if list(one_opt.itertext())[4].split("\n")[0]:
                        name = list(one_opt.itertext())[4].split("\n")[0]
                    elif list(one_opt.itertext())[4].strip():
                        name = list(one_opt.itertext())[4].strip()
                    elif list(one_opt.itertext())[1].strip():
                        name = list(one_opt.itertext())[1].strip()
                    else:
                        name = "分组方案："
                    content = __get_group_info(one_opt)
                elif one_opt.get("class") == "index_id":
                    if one_opt.xpath("a[@class='checked']"):
                        index_list.append(list(one_opt.itertext())[2])
                elif one_opt.get("class") == "show-pair-group":
                    pass
                elif one_opt.text and one_opt.text.strip():
                    ## 指数类型：
                    if one_opt.text.strip() == u'\u6307\u6570\u7c7b\u578b\uff1a':
                        for xx in one_opt.xpath("//li[@class='index_id']/a[@class='checked']"):
                            index_list.append(xx.xpath("*")[0].text)
                    ## 挑选特征物种：
                    elif one_opt.text.strip() == u'\u6311\u9009\u7279\u5f81\u7269\u79cd\uff1a' or one_opt.text.strip() == u'\u9009\u62e9\u7269\u79cd\u4fe1\u606f\uff1a':
                        opt_infos = one_opt.xpath("*")
                        tmp_opt = ""
                        for each in opt_infos:
                            tag_name = each.tag
                            if tag_name == "div":
                                tmp_opt = each
                        if tmp_opt and len(tmp_opt.xpath("ul[@class='chosen-choices']")[0].xpath("li")) > 1:
                            name = one_opt.text.strip()
                            tmp_list = []
                            for y in tmp_opt.xpath("ul[@class='chosen-choices']")[0].xpath("li")[:-1]:
                                tmp_list.append(y.xpath("*")[0].text)
                            content = ";".join(tmp_list)
                        else:
                            name = ''
                            content = ""
                    ## 选择分组：
                    elif one_opt.text.strip() == u'\u9009\u62e9\u5206\u7ec4\uff1a':
                        name = one_opt.text.strip()
                        opt_infos = one_opt.xpath("*")
                        for each in opt_infos:
                            tag_name = each.tag
                            if tag_name == "ul":
                                cat_sp = ''
                                group_infos = each.xpath("li[@class='group_name current']")[0].xpath("*")[1].xpath(
                                    "*")
                                for each_all in group_infos:
                                    this_group = each_all.get("category_name")
                                    sps = [i.xpath("*")[0].text for i in each_all.xpath("*")[1].xpath("*")]
                                    cat_sp += "\n\t\t{}: {}".format(this_group, ", ".join(sps))
                                group_cat = each.xpath("li[@class='group_name current']")[0].xpath("*")[0].text
                                content = group_cat + '\t' + cat_sp
                    ## 选择分类水平：从
                    elif one_opt.text.strip() == u'\u9009\u62e9\u5206\u7c7b\u6c34\u5e73\uff1a\u4ece':
                        name = one_opt.text.strip()
                        opt_infos = one_opt.xpath("*")
                        sel1 = opt_infos[0].xpath("option[@selected]")[0].text
                        sel2 = opt_infos[1].xpath("option[@selected]")[0].text
                        content = sel1 + " 到 " + sel2
                    ## 软件'：
                    elif one_opt.text.strip() == u'\u8f6f\u4ef6\uff1a':
                        pass
                    ## 多样性类型：'
                    elif one_opt.text.strip() == u'\u591a\u6837\u6027\u7c7b\u578b\uff1a':
                        name = one_opt.text.strip()
                        opt_infos = one_opt.xpath("*")
                        sel1 = opt_infos[0].xpath("option[@selected]")[0].text
                        if "Alpha" in opt_infos[0].xpath("option[@selected]")[0].text:
                            sel2 = opt_infos[1].xpath("option[@selected]")[0].text
                        else:
                            sel2 = opt_infos[2].xpath("option[@selected]")[0].text
                        content = sel1 + "\t" + sel2
                    ## 检验方法
                    elif one_opt.text.strip() == u'\u68c0\u9a8c\u65b9\u6cd5\uff1a':
                        name = one_opt.text.strip()
                        opt_infos = one_opt.xpath("*")
                        content = ""
                        if len(opt_infos) > 1:
                            for info in opt_infos[1:]:
                                if info.xpath("option[@selected]"):
                                    if info.get("style") not in ["display:none;","display: none;"]:
                                        content = opt_infos[0].xpath("option[@selected]")[0].text + "\t" + info.xpath("option[@selected]")[0].text
                            if not content:
                                content = __get_content(one_opt)
                        else:
                            content = __get_content(one_opt)
                    ## 控制单位：
                    elif one_opt.text.strip() == u'\u63a7\u5236\u5355\u4f4d\uff1a':
                        name = one_opt.text.strip()
                        opt_infos = one_opt.xpath("*")
                        content = ""
                        for x in opt_infos[0].xpath("span[@class='current']"):
                            if x.get("title"):
                                content += x.get("title") + "\t"
                    ## 物种选择：
                    elif one_opt.text.strip() == u'\u7269\u79cd\u9009\u62e9\uff1a':
                        name = one_opt.text.strip()
                        content = __get_content(one_opt)
                        name1 = one_opt.xpath("span")[-1].text.strip()
                        content1 = __get_content2(one_opt)
                        content1 += u" \u7684\u7269\u79cd"
                    ## 置换次数：
                    elif one_opt.text.strip() == u'\u7f6e\u6362\u6b21\u6570\uff1a':
                        if one_opt.get("style") in ["display:none;","display: none;"]:
                            pass
                        else:
                            name = one_opt.text.split("\n")[1].strip()
                            if name:
                                content = __get_content(one_opt)
                            else:
                                name = one_opt.text.strip()
                                if name:
                                    content = __get_content(one_opt)
                    ## CI计算方法：
                    elif one_opt.text.strip() == u'CI\u8ba1\u7b97\u65b9\u6cd5\uff1a':
                        if one_opt.get("style") in ["display:none;","display: none;"]:
                            pass
                        else:
                            name = one_opt.text.strip()
                            content = one_opt.xpath("*")[0].xpath("option[@selected]")[0].text + "\t" + one_opt.xpath("*")[1].xpath("option[@selected]")[0].text
                    ## Post-hoc检验：
                    elif one_opt.text.strip() == u'Post-hoc\u68c0\u9a8c\uff1a':
                        if one_opt.get("style") in ["display:none;","display: none;"]:
                            pass
                        else:
                            name = one_opt.text.strip()
                            content = one_opt.xpath("*")[0].xpath("option[@selected]")[0].text + "\t" + one_opt.xpath("*")[1].xpath("option[@selected]")[0].text
                    elif one_opt.text.split("\n")[0].strip():
                        if one_opt.get("style") in ["display:none;","display: none;"]:
                            pass
                        else:
                            name = one_opt.text.split("\n")[0].strip()
                            if name:
                                content = __get_content(one_opt)
                    elif one_opt.text.split("\n")[1].strip():
                        if one_opt.get("style") in ["display:none;","display: none;"]:
                            pass
                        else:
                            name = one_opt.text.split("\n")[1].strip()
                            if name:
                                if name == u'\u7f6e\u4fe1\u533a\u95f4\u9608\u503c\uff1a':
                                    content = __get_content(one_opt.xpath("span")[0])
                                else:
                                    content = __get_content(one_opt)
                else:
                    if one_opt.xpath("span"):
                        if one_opt.xpath("span")[0].text:
                            name1 = one_opt.xpath("span")[0].text.strip()
                            content1 = __get_content2(one_opt)
                        elif one_opt.xpath("span")[-1].get("class") in ["top_n_id"]:
                            name1 = one_opt.xpath("span")[-1].text.strip()
                            content1 = __get_content2(one_opt)

                    #else:
                    #    if len(one_opt.text.split("\n")) == 2:
                    #        name = one_opt.text.split("\n")[1]
                    #        content = __get_content(one_opt)
                if name:
                    w.write("{} {}\n".format(name.strip(), content))
                if name1:
                    w.write("{} {}\n".format(name1.strip(), content1))
            if index_list:
                w.write("{} {}\n".format("指数类型：", ",".join(index_list)))
        if all_opts2:
            for one_line in all_opts2:
                opts = one_line.xpath("li")
                for one_opt in opts:
                    name = ''
                    if one_opt.text:
                        name = one_opt.text.split("\n")[0].strip()
                        if name:
                            content = __get_content(one_opt)
                    if name:
                        w.write("{} {}\n".format(name.strip(), content))


def __get_group_info(opt):
    cat_sp = ''
    group_infos = opt.xpath("ul[@class='choose-list']")[0]
    for each in group_infos:
        this_group = each.xpath("*/span")[0].text
        sps = [i.text for i in each.xpath("*/div[@class='popup']/a[@class='checked']")]
        if sps:
            cat_sp += "\n\t\t{}: {}".format(this_group, ", ".join(sps))
        #else:
        #    cat_sp = ""
    group_cat = opt.xpath("*/div[@class='select-input']/input[@type='text']")[0].get("value")
    if cat_sp:
        content = group_cat + '\t' + cat_sp
    elif opt.xpath("ul") and opt.xpath("ul")[0].get("class") == "choose-list":
        group_infos = opt.xpath("ul")[0].xpath("li")
        for each in group_infos:
            this_group = each.xpath("*/span")[0].text
            sps = [i.text for i in each.xpath("*/div[@class='popup']/a[@class='checked']")]
            if sps:
                cat_sp += "\n\t\t{}: {}".format(this_group, ", ".join(sps))
        content = group_cat + '\t' + cat_sp
    else:
        content = "无"
    return content

def __get_env_list(div):
    env_spans = div.xpath("span[@class='show_span current']")
    contents = []
    if env_spans:
        pass
    else:
        env_spans = div.xpath("span[@class='show_span current unclick']")
    for env in env_spans:
        contents.append(env.text)
    if not contents:
        contents = u"\u65e0\u000d\u000a"
    return "\t".join(contents)

def __get_content(opt):
    content = ""
    content1 = ""
    content2 = ""
    env = False
    opt_infos = opt.xpath("*")
    for each in opt_infos:
        tag_name = each.tag
        if tag_name == "select":
            sel = each.xpath("option[@selected]") or each.xpath("option")
            content = sel[0].text
        elif tag_name == 'input':
            content = each.get("value")
        elif tag_name == 'label':
            checked = each.xpath("input[@checked]")
            if checked:
                content = list(each.itertext())[1].split("\n")[0]
            else:
                checked = each.xpath("input")
                if len(checked) >1:
                    if checked[1].get("value"):
                        content = checked[1].get("value")
        elif tag_name == "div":
            class1 = each.get("class")
            if class1 == "env-list":
                content = __get_env_list(each)
                content2 = content
            elif class1 == "chosen-container chosen-container-multi":
                each.xpath("a[@class='chosen-search-input']")
            else:
                content = __multi_filter(opt) or (each.xpath("*//input[@type='text']")[0].get("value") if each.xpath("*//input[@type='text']") else "")
                try:
                    if each.xpath("*//input[@type='text']")[0].get("name") == "env_name":
                        content1 =content
                        env = True
                    elif each.xpath("*//input[@type='text']")[0].get("id") == "env_name":
                        content1 =content
                        env = True
                except:
                    pass
                #print each.xpath("*//input[@type='text']")[0].get("name")
            print content
            print env
        if content in ["Scheffe", "Welch's (uncorrected)", "Tukey-kramer", "Games-Howell","bootstrap"
                       ,"DP：AsymptoticCC","DP：Asymptotic","DP：Newcombe-Wilson"]:
            content1 = content
        elif content1 and content2:
            content = content1 + ":(" + content2 + ")"
            break
        elif content and (content2 == "" and env ==False):
            if content1 and content1 in ["Scheffe", "Welch's (uncorrected)", "Tukey-kramer", "Games-Howell","bootstrap"
                                         ,"DP：AsymptoticCC","DP：Asymptotic","DP：Newcombe-Wilson"]:
                content = content1 + ";" + content
            break
    return content

def __get_content2(opt):
    content = ''
    opt_infos = opt.xpath("*")
    for aa in opt_infos:
        for each in aa:
            tag_name = each.tag
            if tag_name == "select":
                sel = each.xpath("option[@selected]") or each.xpath("option")
                content = sel[0].text
            elif tag_name == 'input':
                content = each.get("value")
            elif tag_name == 'label':
                checked = each.xpath("input[@checked]")
                if checked:
                    content = list(each.itertext())[1].split("\n")[0]
            elif tag_name == "div":
                class1 = each.get("class")
                if class1 == "env-list":
                    content = __get_env_list(each)
                else:
                    content = __multi_filter(opt) or (
                        each.xpath("*//input[@type='text']")[0].get("value") if each.xpath(
                            "*//input[@type='text']") else "")
            if content and content != "origin_env_table":
                break
    return content

def __multi_filter(opt):
    print(opt)
    mul_filter = opt.xpath("*//div[@class='siftings-item clear current']")
    content = ''
    if mul_filter:
        for one_filter in mul_filter:
            spans = one_filter.xpath("*/*//span")
            content += "\t\t"
            for span in spans:
                span_text = __all_text(span.itertext())
                span_infos = span.xpath('*')
                conts = []
                for one in span_infos:
                    c = ''
                    one_text = __all_text(one.itertext())
                    span_text = span_text.replace(one_text, '', 1).replace("\t\t", "\t")
                    tag_name = one.tag
                    if tag_name == "select":
                        sel = one.xpath("option[@selected]") or one.xpath("option")
                        c = "{} {}".format(__get_text(one.text), __get_text(sel[0].text))
                    elif tag_name == "input":
                        c = "{} {}".format(one.get("value"), __get_text(one.text))
                    conts.append(c)
                span_text = span_text.strip().split('\t')
                content += " {}".join(span_text)
                if len(span_text) == 1:
                    content += " {} "
                content = content.format(*conts)
            content += '\n'
    return content


def __get_text(text):
    return (text or '') and text.split("\n")[0]


def __all_text(textiter):
    span_text = [re.sub(r'[^\S]', '', i) for i in textiter]
    span_text = filter(lambda x: x, span_text)
    return "\t".join(span_text)


