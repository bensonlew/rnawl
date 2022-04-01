# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
import os,pandas as pd
from biocluster.core.exceptions import FileError
import pickle
import time
from bson.son import SON
from mbio.packages.metagenomic.common import link_file, link_dir


def set_step(event):
    if event['data'].has_key("start"):
        event['data']['start'].start()
    if event['data'].has_key("end"):
        event['data']['end'].finish()
    event['data']['bind_step'].update()


def set_run(bind_object, opts, module, event, step, start=True):
    module.set_options(opts)
    module.on('start', set_step, {'start': step, 'bind_step': bind_object.step})
    module.on('end', set_step, {'end': step, 'bind_step': bind_object.step})
    module.on('end', bind_object.set_output, event)
    if start:
        module.run()

def mv_file(old_file, new_dir):
    file_name = os.path.basename(old_file)
    link_file(old_file, os.path.join(new_dir, file_name))

def check_raw_dir(bind_object):
    raw_dir_path = os.path.dirname(bind_object.option("raw_dir_json").prop["path"]) + "/raw_dir"
    if not os.path.isdir(raw_dir_path):
        os.mkdir(raw_dir_path)
    list_path = os.path.join(raw_dir_path, "list.txt")
    # 将raw_dir_json的文件全部转移到raw_dir路径下
    with open(list_path, "w") as f:
        f.write("Sample Name\tFile Name\tInsert Size(bp)\tReads Length(bp)\tGenome Size(Mb)\tLibrary\tLibrary Name\n")
        for sample in bind_object.option("raw_dir_json").samples:
            f.write(sample.name + "\t")
            bind_object.logger.info(sample.insert_size)
            bind_object.logger.info(sample.read_length)
            bind_object.logger.info(sample.genome_size)
            bind_object.logger.info(sample.library)
            bind_object.logger.info(sample.library_name)
            if sample.r1_file:
                path = sample.r1_file.prop["path"]
                f.write(os.path.basename(path))
                mv_file(path, raw_dir_path)
                if sample.r2_file:
                    path = sample.r2_file.prop["path"]
                    f.write("," + os.path.basename(path))
                    mv_file(path, raw_dir_path)
            else:
                file_list = []
                for file in sample.files:
                    path = file.prop["path"]
                    file_list.append(os.path.basename(path))
                    mv_file(path, raw_dir_path)
                f.write(",".join(file_list))
            for i in [sample.insert_size, sample.read_length, sample.genome_size, sample.library, sample.library_name]:
                i = i if i else "-"
                f.write("\t" + i)
            f.write("\n")
    bind_object.option("raw_dir").set_path(os.path.dirname(list_path))
    bind_object.option("raw_dir").check()


def usable_file(path):
    if os.path.isfile(path) and os.path.getsize(path) != 0:
        return True
    else:
        return False


def mark_major(data, genome_info, output):
    """
    获取pe reads，找到不同文库间碱基数最高的那一组
    :param data:
    :param output:
    :return:
    """
    data = data[data["Library"] == "PE"]
    if len(data) == 0:
        return False
    new_data = pd.DataFrame(columns=data.columns)
    samples = set(data["Sample Name"])
    for sample in samples:
        sample_data = data[data["Sample Name"] == sample]
        sample_data.loc[:,"genome"] = genome_info[sample][1]
        if len(sample_data) == 1:
            new_data = new_data.append(sample_data)
        elif len(sample_data) > 1:
            new_data = new_data.append(sample_data.sort_values("total bases(bp)", ascending=False)[:1])
    new_data.reindex(columns=["Sample_lib", "total bases(bp)", "genome"]).to_csv(output, sep="\t", index=False)
    return True

def mark_qc(data, genome_info, output):
    # 二代数据增加预计的基因组大小
    def get_genome(df):
        df["genome"] = genome_info[df["Sample Name"]][1]
        return df
    data = data.apply(get_genome, axis=1)
    data.to_csv(output, sep="\t", index=False)

def change_qc_list(data, genome_info, output):
    def get_genome(df):
        df["size"] = genome_info[df["Sample Name"]][1]
        return df
    if os.path.isfile(output):
        os.remove(output)
    data["file"] = data.Sample_lib + ".clean.1.fq," + data.Sample_lib + ".clean.2.fq"
    data = data.apply(get_genome, axis=1)
    data = data.reindex(columns=["Sample Name", "file", "Insert Size(bp)", "Read Len(bp)", "size", "Library"])
    data.columns = ["sample", "file", "insert", "length", "size", "lib"]
    data.to_csv(output, sep="\t", index=False)

def add_genome(table, genome_info, output):
    # 三代数据表增加预计的基因组大小
    def get_genome(df):
        df["genome"] = genome_info[df["Sample Name"]][1]
        return df
    data = pd.read_table(table, header=None, names=["Sample Name", "Library", "File"])
    data = data.apply(get_genome, axis=1)
    data.to_csv(output, sep="\t", index=False)

def add_all_genome(genome_info, output):
    genome_data = {}
    for sample in genome_info:
        genome_data[sample] = genome_info[sample][1]
    data = pd.DataFrame({"genome": genome_data})
    data.index.name = "Sample Name"
    data.to_csv(output, sep="\t", index=True)

def save_sample_info(reads_qc, path):
    sample_info = {
        "all_samp": dict(zip(reads_qc.all_samp, reads_qc.all_samp)),
        "pe_samp": reads_qc.pe_samp,
        "chr_samp": reads_qc.chr_samp
    }
    with open(path, "w") as file:
        pickle.dump(sample_info, file)

def load_pickle(path):
    with open(path, "r") as file:
        file_json = pickle.load(file)
        for i,j in file_json.items():
            if i == "library":
                continue
            elif i == "all_samp":
                sorted_j = sorted(j.items(), key=lambda d:d[0])
                file_json[i] = SON(sorted_j)
            else:
                file_json[i] = sorted(j)
    return file_json

def wait_end(bind_object, samples, max_wait_time=10):
    if max_wait_time <= 0:
        if bind_object.draft_assess.is_end:
            bind_object.logger.info("Active set draft_assess output")
            link_dir(bind_object.draft_assess.output_dir, bind_object.output_dir + "/genomic_assessment")
        bind_object.logger.info("Active exit wait_end after 10 times...")
        return
    dir_path = bind_object.output_dir + "/genomic_assessment"
    if not os.path.isdir(dir_path) or set(samples) != set(os.listdir(dir_path)):
        bind_object.logger.info("waiting end 5s...")
        bind_object.logger.info(samples)
        if os.path.isdir(dir_path):
            bind_object.logger.info(set(os.listdir(bind_object.output_dir + "/genomic_assessment")))
        time.sleep(5)
        wait_end(bind_object, samples, max_wait_time-1)
