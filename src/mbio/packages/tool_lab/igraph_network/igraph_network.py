# coding=utf-8
import subprocess
import argparse
import pandas as pd

"""

"""
# get arguments
parser = argparse.ArgumentParser(description="igraph:")
parser.add_argument('-nodes_file', type=str, metavar="nodes_file", required=True,
                    help="node infos detail file")
parser.add_argument('-edges_file', type=str, metavar="edges_df", required=True,
                    help="edge infos detail file")
parser.add_argument('-layout', type=str, required=True, metavar="layout",default="layout_in_circle",
                    help="plot_layout")
parser.add_argument('-node_style', type=str, metavar="node_style", default="circle",
                    help="node_style, one in ['none','circle','square','csquare','rectangle','crectangle', 'vrectangle','pie','raster','sphere']")
parser.add_argument('-line_style',  type=int, metavar="line_style", default=1,
                    help="line_style,0:'blank',1:'solid',2:'dashed',3:'dotted',4:'dotdash',5:'longdash',6:'twodash'")
args = parser.parse_args()



r_cmds = \
"""
# igraph network
library('ggplot2')
library('igraph')
nodes_df = read.table('{nodes_file}', header=T, sep="\\t")
edges_df = read.table('{edges_file}', header=T, sep="\\t")
g <- graph_from_data_frame(d = edges_df, vertices = nodes_df,directed = F)
# nodes_types_num <-  length(levels(as.factor(nodes_df$type.label)))
# vcolor<- categorical_pal(nodes_types_num)
print(g)
V(g)$color = V(g)$plot_color
# V(g)$color<-vcolor[factor(V(g)type)]
V(g)$size =  V(g)$size
E(g)$width <- E(g)$weight 
E(g)$color <- E(g)$plot_color
V(g)$label.dist = V(g)$size/7
pdf("igraph.pdf")
plot(g, layout = {layout},
     vertex.shape ='{node_style}',
     edge.lty = {line_style},
     )
legend(x=-1.5,y=-1.1,rle(V(g)$group)$values,pch=21,col="#777777",pt.bg=rle(V(g)$color)$values)
dev.off()
png("igraph.png")
plot(g, layout = {layout},
     vertex.shape ='{node_style}',
     edge.lty = {line_style},
     )
legend(x=-1.5,y=-1.1,rle(V(g)$group)$values,pch=21,col="#777777",pt.bg=rle(V(g)$color)$values)
dev.off()
svg("igraph.svg")
plot(g, layout = {layout},
     vertex.shape ='{node_style}',
     edge.lty = {line_style},
     )
legend(x=-1.5,y=-1.1,rle(V(g)$group)$values,pch=21,col="#777777",pt.bg=rle(V(g)$color)$values)
dev.off()
 """.format(
        nodes_file = args.nodes_file,
        edges_file = args.edges_file,
        layout=args.layout,
        node_style=args.node_style,
        line_style=args.line_style,
    )

with open('wgcna_relate_analysis.r', 'w') as f:
    f.write(r_cmds)
subprocess.check_call("Rscript wgcna_relate_analysis.r", shell=True)


