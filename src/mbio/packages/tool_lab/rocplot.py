
# ~/app/bioinfo/tool_lab/miniconda3/bin/python
# ~/app/bioinfo/tool_lab/miniconda3/bin/pip install matplotlib -i https://pypi.tuna.tsinghua.edu.cn/simple
# ~/app/bioinfo/tool_lab/miniconda3/bin/pip install pandas -i https://pypi.tuna.tsinghua.edu.cn/simple
# ~/app/bioinfo/tool_lab/miniconda3/bin/pip install sklearn -i https://pypi.tuna.tsinghua.edu.cn/simple
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from itertools import cycle
import pandas as pd
import os

def plotroc(infile,outdir):
	# colors = cycle(['blue', 'red', 'green'])
	cycol = cycle('bgrcmk')
	lw = 1.75
	df = pd.read_table(infile,sep="\t")
	n_classes = df.shape[1]-1
	df2_array = pd.DataFrame(zip(*[df.iloc[:,0]]*int(df.shape[1]-1))).to_numpy()
	df3_array = df.iloc[:,1:].to_numpy()
	df3_columns = df.columns.tolist()[1:]
	fpr = dict()
	tpr = dict()
	roc_auc = dict()
	for i in range(n_classes):
	    fpr[i], tpr[i], _ = roc_curve(df2_array[:, i], df3_array[:, i])
	    roc_auc[i] = auc(fpr[i], tpr[i])
	for k,v in roc_auc.items():
	    roc_auc[k] = 1-v
	# for i, color in zip(range(n_classes), colors):
	plt.figure(figsize=(9, 7))
	for i in range(n_classes):
	    plt.plot(tpr[i],fpr[i], color=next(cycol), lw=lw,
	             label='ROC of class {0} (AUC = {1:0.2f})'
	             ''.format(df3_columns[i], roc_auc[i]))
	plt.plot([0, 1], [0, 1], 'k--', lw=lw)
	plt.xlim([-0.05, 1.0])
	plt.ylim([0.0, 1.05])
	plt.xlabel('1-Specificity')
	plt.ylabel('Sensitivity')
	plt.title('Receiver operating characteristic')
	plt.legend(loc="lower right")
	# plt.show()
	plt.savefig(os.path.join(outdir,"roc.pdf"))  # or you can pass a Figure object to pdf.savefig
	plt.close()


if __name__ == "__main__":
    import argparse
    print("参数解析开始")
    parser = argparse.ArgumentParser(description='This script is used to plot roc curve')
    parser.add_argument('-infile', type=str, default="", help='infile')
    parser.add_argument('-outdir', type=str, default="", help='outdir')
    args = parser.parse_args()
    plotroc(args.infile, args.outdir)
    print('成功绘制了ROC曲线')