# -*- coding: utf-8 -*-
# !/usr/bin/env python
# __author__ = 'guhaidong'

import argparse
import numpy as np
import pandas as pd
from sklearn import model_selection
# from ml import MlGroup


def option():
    par = argparse.ArgumentParser()
    par.add_argument("-m", metavar='[method]',
                     choices=['lr', 'lda', 'svm', 'bayes', 'tree', 'knn', 'randomforest', 'gradient', 'sgd', 'ada'],
                     required=True, help='classifier method: svm,bayes,tree,knn,randomforest,gradient等')
    par.add_argument("-i", metavar='[train_file]', required=True, help='input train file')
    par.add_argument("-g", metavar='[group_table]', help='input group_table for train, txt format')
    par.add_argument("-model", metavar='[model]', required=True, help='train result model for classify')
    par.add_argument("-disease", metavar='[disease]', help='disease')
    args = par.parse_args()
    return args


class Train(object):
    def __init__(self, **argv):
        super(Train, self).__init__()
        self.clf = ""  # 模型
        self.x = ""
        self.y = ""
        self.x_train = ""
        self.y_train = ""
        self.x_test = ""
        self.y_test = ""
        self.method = argv['method']
        self.file = argv['file']
        self.model = argv['model']
        self.disease = argv['disease']
        print "=======CHECK GROUP FILE========"
        if 'group' in argv.keys():
            self.group = argv['group']
            print "yes,%s" % self.group
        else:
            self.group = False
            print "no group file,%s" % self.group
        marker_crc = [
            "d__Bacteria;k__unclassified_d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Roseburia;s__Roseburia_hominis",
            "d__Bacteria;k__unclassified_d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium_sp._CAG:465",
            "d__Bacteria;k__unclassified_d__Bacteria;p__Fusobacteria;c__Fusobacteriia;o__Fusobacteriales;f__Fusobacteriaceae;g__Fusobacterium;s__Fusobacterium_nucleatum",
            "d__Bacteria;k__unclassified_d__Bacteria;p__Verrucomicrobia;c__Verrucomicrobiae;o__Verrucomicrobiales;f__Akkermansiaceae;g__Akkermansia;s__Akkermansia_muciniphila",
            "d__Bacteria;k__unclassified_d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides_fragilis",
            "d__Bacteria;k__unclassified_d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium_hathewayi_CAG:224",
            "d__Bacteria;k__unclassified_d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Lachnoclostridium;s__[Clostridium]_symbiosum"
        ]  # 需配置,CRC的marker
        marker_t2d = [
            "d__Bacteria;k__unclassified_d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Ruminococcus;s__Ruminococcus_obeum_CAG:39",
            "d__Bacteria;k__unclassified_d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Roseburia;s__Roseburia_inulinivorans_CAG:15",
            "d__Bacteria;k__unclassified_d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Roseburia;s__Roseburia_inulinivorans",
            "d__Bacteria;k__unclassified_d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Roseburia;s__Roseburia_intestinalis_CAG:13",
            "d__Bacteria;k__unclassified_d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides_pectinophilus_CAG:437",
            "d__Bacteria;k__unclassified_d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Roseburia;s__Roseburia_intestinalis",
            "d__Bacteria;k__unclassified_d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Faecalibacterium;s__Faecalibacterium_prausnitzii",
            "d__Bacteria;k__unclassified_d__Bacteria;p__Spirochaetes;c__Spirochaetia;o__Spirochaetales;f__Spirochaetaceae;g__Spirochaeta;s__Spirochaeta_smaragdinae"
        ]
        marker_obesity = [
            "d__Bacteria;k__unclassified_d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Faecalibacterium;s__Faecalibacterium_prausnitzii",
            "d__Bacteria;k__unclassified_d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia_coli",
            "d__Bacteria;k__unclassified_d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Ruminiclostridium;s__[Clostridium]_leptum",
            "d__Bacteria;k__unclassified_d__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Bifidobacteriales;f__Bifidobacteriaceae;g__Bifidobacterium;s__Bifidobacterium_pseudocatenulatum"
        ]
        self.marker_list = {
            'crc': marker_crc,
            't2d': marker_t2d,
            'obesity': marker_obesity
        }
        self.index_list = []  # 储存用来建模的物种

    def get_train_data(self):
        data = self.get_txt()
        # self.x,self.y = np.split(data,[len(data[0])-1,],axis=1)  # 最后一列为分类标签，此时以数字表示
        self.y, self.x = np.split(data, [1, ], axis=1)  # 第一列为分类标签，此时以数字显示
        self.x_train, self.x_test, self.y_train, self.y_test = model_selection.train_test_split(self.x, self.y, test_size=0.2)

    def get_txt(self):
        from ml import MlGroup
        raw = pd.read_table(self.file, header=None,low_memory=False, index_col=0)
        raw_mark = pd.read_table(self.file, low_memory=False, index_col=0)
        ml_group = MlGroup(self.group)
        # raw_top = raw.replace(ml_group.desease_hash).drop([0], axis=1).astype(float)
        raw_top = raw.replace(ml_group.desease_hash).astype(float)
        if self.marker_list[self.disease] != []:
            # tttttt = ml_group.get_sample_list(['Healthy'])
            # yyyyyy = self.marker_list[self.disease]
            # uuuuuu = raw_mark[tttttt]
            # raw_mark = uuuuuu.loc[yyyyyy]
            raw_mark = raw_mark[ml_group.get_sample_list(['Healthy'])].reindex(self.marker_list[self.disease], fill_value=0)
            mark_file = self.model.split("_")[0] + "_mark.xls"
            raw_mark.to_csv(mark_file, sep="\t")
        # raw_top_marker = raw_top.loc[self.marker_list]
        self.index_list = raw_top.index.tolist()[1:]
        data = raw_top.values.T
        # print data
        return data

    def train(self):
        if self.method == 'lr':
            from sklearn.linear_model import LogisticRegression
            self.clf = LogisticRegression()
        elif self.method == 'lda':
            from sklearn.decomposition import LatentDirichletAllocation
            self.clf = LatentDirichletAllocation(n_topics=30, max_iter=50, learning_method='batch')
        elif self.method == 'svm':
            from sklearn import svm
            self.clf = svm.SVC(C=0.8, kernel='rbf', gamma=20, decision_function_shape='ovr')
        elif self.method == 'bayes':
            from sklearn import naive_bayes
            self.clf = naive_bayes.GaussianNB()
        elif self.method == 'tree':
            from sklearn import tree
            self.clf = tree.DecisionTreeClassifier(
                criterion='gini')  # criterion = ['gini', 'information gain', 'chi-square']
        elif self.method == 'knn':
            from sklearn.neighbors import KNeighborsClassifier
            self.clf = KNeighborsClassifier(n_neighbors=2)
        elif self.method == 'randomforest':
            from sklearn.ensemble import RandomForestClassifier
            self.clf = RandomForestClassifier(n_estimators=600, oob_score=True, criterion="gini", min_samples_split=10, max_depth=10)
        elif self.method == 'gradient':
            from sklearn.ensemble import GradientBoostingClassifier
            self.clf = GradientBoostingClassifier(n_estimators=100, learning_rate=1.0, max_depth=1, random_state=0)
        elif self.method == 'sgd':
            from sklearn.linear_model import SGDClassifier
            self.clf = SGDClassifier(loss="log", alpha=0.01, max_iter=200, fit_intercept=True)
        elif self.method == 'ada':
            from sklearn.ensemble import AdaBoostClassifier
            self.clf = AdaBoostClassifier(n_estimators=100)
        fit = self.clf.fit(self.x_train, self.y_train.ravel())
        if self.method == 'randomforest':
            self.export_index(fit)

    def export_index(self, return_fit):
        print "Export Feature Importance:"
        importance_table = pd.DataFrame(data={"gene": self.index_list,"importances": return_fit.feature_importances_.tolist()})
        importance_table.sort_values(by='importances', ascending=False, inplace=True)
        index_file = self.model.split("_")[0] + "_index"
        importance_table.to_csv(index_file, sep="\t", index=False)

    def test_score(self):
        from sklearn import metrics
        y1_predprob = self.clf.predict_proba(self.x_train)[:,1]
        y2_predprob = self.clf.predict_proba(self.x_test)[:, 1]
        y1_pred = self.clf.predict(self.x_train)
        y2_pred = self.clf.predict(self.x_test)
        print "=========Model Test Report==========="
        print "AUC Score (Train): %f" % metrics.roc_auc_score(self.y_train.ravel(), y1_predprob)
        print "AUC Score (Test): %f" % metrics.roc_auc_score(self.y_test.ravel(), y2_predprob)
        fpr,tpr,thresholds = metrics.roc_curve(self.y_test.ravel(), y2_predprob, pos_label=1)
        # fpr,tpr,thresholds = metrics.roc_curve(self.y_train.ravel(), y1_predprob, pos_label=1)
        roc_coord = pd.DataFrame({"x": fpr, "y": tpr})
        roc_coord.to_csv(self.model + '.auc', sep="\t", index=0)
        report = metrics.classification_report(self.y_test.ravel(),y2_pred, target_names=["class 0", "class 1"])
        with open(self.model + ".report", "w") as f:
            f.write(report)
            f.write("AUC Score (Train): %f" % metrics.roc_auc_score(self.y_train.ravel(), y1_predprob))
            f.write("\n")
            f.write("AUC Score (Test): %f" % metrics.roc_auc_score(self.y_test.ravel(), y2_predprob))

    def save_model(self):
        from sklearn.externals import joblib
        joblib.dump(self.clf, self.model, compress=3)

    def save_pic(self):
        # 做图片，这里不要传数据了
        import ml
        save_file = self.model + '.draw'
        # ml.combine(self.clf.predict_proba(self.x_train)[:,1], self.y_train.ravel(), save_file)
        ml.combine(self.clf.predict_proba(self.x_train), self.y_train.ravel(), save_file)

    def open_model(self):
        from sklearn.externals import joblib
        self.clf = joblib.load(self.model)

    def predict(self):
        result = [self.clf.predict(self.x_test), self.clf.predict_proba(self.x_test)[:,1]]
        return result

    def run(self):
        self.get_train_data()
        self.train()
        self.test_score()
        self.save_model()
        self.save_pic()

if __name__ == "__main__":
    opts = option()
    train_object = Train(method=opts.m, file=opts.i, model=opts.model, group=opts.g, disease=opts.disease)
    print train_object.run()
