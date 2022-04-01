# -*- coding: utf-8 -*-
# __author__ = 'xuxi'

import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
import os
import sys

def peperr_image(dmass_infile, outdir):
    dmass_df = pd.read_table(dmass_infile,sep="\t")
    x = dmass_df.iloc[:,0].array
    y = dmass_df.iloc[:,1].array
    mpl.pyplot.figure(figsize=(16,10))
    mpl.rcParams['axes.spines.right'] = False
    mpl.rcParams['axes.spines.top'] = False
    mpl.pyplot.scatter(x, y, s=20, c="#5ca35f", alpha=0.5)
    mpl.pyplot.xlabel("m/z(Da)\n", fontsize=18)
    mpl.pyplot.ylabel("DeltaM(ppm)\n", fontsize=18)
    mpl.pyplot.title("The distribution of peptide matching error\n", fontsize=17,fontweight='bold')
    mpl.pyplot.savefig(os.path.join(outdir, "dMass.pdf"))
    mpl.pyplot.savefig(os.path.join(outdir, "dMass.png"), dpi=100)


if __name__ == "__main__":
    peperr_image(sys.argv[1], sys.argv[2])