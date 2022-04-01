
importance = c(${importance})
pvalue = ${pvalue}
df = read.table('merge',sep='\t',row.names=1,header=T)
group_header = "${group_header}"


procruste <- function(df, importance, pc_a, pc_b, monte_p, group=NULL) {{
# .libPaths('/mnt/ilustre/users/hongming.wang/R_libs')
library( ggplot2 )
library( grid )
library(RColorBrewer)

my_theme <- theme( panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   axis.line = element_line( colour = '#3B3B3B',
                                             size=0.5)) +
                   theme( text = element_text(size=12)) +
                   theme(plot.margin = unit(c(1,1,1,1), "cm")) +
                   theme(axis.text=element_text(colour="#3B3B3B"))

pdf('Procrustes.pdf', width=7.2, height=7)
names(importance) = c('PC1', 'PC2', 'PC3')
df['tmp1'] = (df[,pc_a] + df[, paste( pc_a, '_t', sep='' )])/2
df['tmp2'] = (df[,pc_b] + df[, paste( pc_b, '_t', sep = '' )])/2
p <- ggplot( data = df )
#p <- p + geom_segment(aes_string(x=pc_a, y=pc_b, xend='tmp1',yend='tmp2'), color='#FFC1C1', size=0.6)
#p <- p + geom_segment(aes_string(x=paste( pc_a, '_t', sep='' ), y=paste( pc_b, '_t', sep = '' ), xend='tmp1',yend='tmp2'), color='#CDC9C9', size=0.6)

   if (is.null(group)) {{
   p <- p + geom_segment(aes_string(x=pc_a, y=pc_b, xend=paste( pc_a, '_t', sep='' ),yend=paste( pc_b, '_t', sep = '' ) ), color='#CDC9C9', size=0.6)
   p <- p + geom_point(aes_string(x=pc_a, y=pc_b), size=3, colour = '#00BAAE' )
#   p <- p + geom_point(aes_string(x= paste( pc_a, '_t', sep='' ), y= paste( pc_b, '_t', sep = '' )), size=3, colour = '#00BAAE' )
   }}
   else {{
   p <- p + geom_point(aes_string(x=pc_a, y=pc_b, colour = group), size=3)
   p <- p + geom_segment(aes_string(x=pc_a, y=pc_b, xend=paste( pc_a, '_t', sep='' ),yend=paste( pc_b, '_t', sep = '' ),color = group ), size=0.6)
#   p <- p + geom_point(aes_string(x= paste( pc_a, '_t', sep='' ), y= paste( pc_b, '_t', sep = '' ), colour = group), size=3)
   }}

p <- p + ylab( paste( pc_b,'  ', importance[pc_b], '%\n', sep='' ) )

pValue <- as.numeric(monte_p)
   if (pValue < 0.001) {{
     p <- p + xlab(paste( '\n', pc_a,'  ', importance[pc_a], '%\n\n', 'Monte Carlo Label Permutations P<0.001', sep='' ) )
 	}}
   else if (pValue < 0.01) {{
     p <- p + xlab(paste( '\n', pc_a,'  ', importance[pc_a], '%\n\n', 'Monte Carlo Label Permutations P<0.01', sep='' ) )
	}}
   else if (pValue < 0.05) {{
     p <- p + xlab(paste( '\n', pc_a,'  ', importance[pc_a], '%\n\n', 'Monte Carlo Label Permutations P<0.05', sep='' ) )
	}}
   else {{
     p <- p + xlab(paste( '\n', pc_a,'  ', importance[pc_a], '%\n\n', 'Monte Carlo Label Permutations\nP=', monte_p, sep='' ) )
	}}
{theme}
p <- p + my_theme
    if ( (! is.null(group) ) && (length( unique( df[,group] ) ) <= 9) )
        {{
	p <- p + scale_colour_brewer(palette = "Set1")
        }}

p<- p + scale_x_continuous( limits = c( min( c(df[,pc_a], df[, paste( pc_a, '_t', sep = '' )]) ), 1.4*max( c(df[,pc_a], df[, paste( pc_a, '_t', sep = '' )]) ) ) )
p<- p + scale_y_continuous( limits = c( min( c(df[,pc_b], df[, paste( pc_b, '_t', sep = '' )]) ), 1.4*max( c(df[,pc_b], df[, paste( pc_b, '_t', sep = '' )]) ) ) )
print(p)
dev.off()
}}

procruste(df, importance, 'PC1', 'PC2', pvalue, group_header)