# # # # ~/app/bioinfo/rna/miniconda2/bin/R

# options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# install.packages("ggplot2")
# install.packages("reprex")
# install.packages("tidyverse")


library('getopt');

spec = matrix(c(
  'infile','i',1,'character',
  'outdir','o',1,'character',
  'help','h',0,'logical'
  ), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat(" 
Usage:
  --infile(-i)   the input file
  --outdir(-o)   the output directory
  --help(-h)    usage
\n")
  q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$infile)) { print_usage(spec)}
if ( is.null(opt$outdir)) { print_usage(spec)}
# if ( is.null(opt$method)){ opt$method="BCCC()"}
# if ( is.null(opt$preprocess)){ opt$preprocess=FALSE}else{opt$preprocess=TRUE}
# if ( is.null(opt$heatmap)){ opt$heatmap=FALSE }else{opt$heatmap=TRUE}
# if ( is.null(opt$ngenes)){ opt$ngenes=1000}
# if ( is.null(opt$missing)){ opt$missing="geneMedian"}
# if ( is.null(opt$format)){ opt$format=1}
# if ( is.null(opt$transform)){ opt$transform=1}

inFile = opt$infile
outdir = opt$outdir

print(paste("inFile： ",inFile))
print(paste("outdir ",outdir))

#参考https://blog.csdn.net/qq_42374697/article/details/111110555

library(ggplot2)
library(tidyverse)
# install.packages('ragg')
# library(ragg)

data = as.data.frame((read.table(inFile,header=T,sep="\t")))

if (length(colnames(data)) == 3){
  origin_colnames = colnames(data)
  colnames(data) <- c("individual", "group","value")
  origin_max_value <- max(data$value)
  origin_min_value <- min(data$value)
  data$value <- data$value/(max(data$value)/100)
  # 对数据做转换
  datagroup <- data$group %>% unique() #"A" "B" "C" "D"

  allplotdata <- tibble('group' = datagroup,
                        'individual' = paste0('empty_individual_', seq_along(datagroup)),
                        'value' = 0) %>% 
    bind_rows(data) %>% arrange(group) %>% mutate(xid = 1:n()) %>% 
    mutate(angle = 90 - 360 * (xid - 0.5) / n()) %>% 
    mutate(hjust = ifelse(angle < -90, 1, 0)) %>% 
    mutate(angle = ifelse(angle < -90, angle+180, angle)) 

  # 提取出空的数据，做一些调整
  firstxid <- which(str_detect(allplotdata$individual, pattern = "empty_individual")) # 1 12 43 58

  segment_data <- data.frame('from' = firstxid + 1,
                             'to' = c(c(firstxid - 1)[-1], nrow(allplotdata)),
                             'label' = datagroup) %>% 
    mutate(labelx = as.integer((from + to)/2))


  # 自定坐标轴 
  coordy <- tibble('coordylocation' = seq(from = min(allplotdata$value), to = max(allplotdata$value), (max(allplotdata$value)-min(allplotdata$value))/10),
         'coordylocation_origin' = seq(from = 0, to = origin_max_value, origin_max_value/10),
         'coordytext' = as.character(round(coordylocation_origin, 2)),
         'x' = 1)

  # 自定义坐标轴的网格
  griddata <- expand.grid('locationx' = firstxid[-1], 'locationy' = coordy$coordylocation)

  # 绘图
  p <- ggplot() + 
    geom_bar(data = allplotdata, aes(x = xid, y = value, fill = group), stat = 'identity') + 
    geom_text(data = allplotdata %>% filter(!str_detect(individual, pattern = "empty_individual")), 
              aes(x = xid, label = individual, y = value+10, angle = angle, hjust = hjust),
              color="black", fontface="bold",alpha=0.6, size=2.5) + 
    geom_segment(data = segment_data, aes(x = from, xend = to), y = -5, yend=-5) + 
    geom_text(data = segment_data, aes(x = labelx, label = label), y = -15) + 
    geom_text(data = coordy, aes(x = x, y = coordylocation, label = coordytext),
              color="grey", size=3 , angle=0, fontface="bold") + 
    geom_segment(data = griddata, 
                 aes(x = locationx-0.5, xend = locationx + 0.5, y = locationy, yend = locationy),
                 colour = "grey", alpha=0.8, size=0.6) + 
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(limits = c(-50,130)) + 
    coord_polar() +
    theme_void() +
    guides(fill=guide_legend(title=origin_colnames[[2]])) +
    theme(legend.position = 'right',plot.margin = unit(c(0,1,0,0), "cm"))

  # p
  # 保存
  ggsave(filename = file.path(outdir, 'circular_bar.pdf'), plot = p, width = 13, height = 13)
}
if (length(colnames(data)) == 2){
  colnames(data) <- c("individual","value")
  data$value <- data$value/(max(data$value)/100)
  # ----- This section prepare a dataframe for labels ---- #
  # Get the name and the y position of each label
  data$id <- seq.int(nrow(data))
  label_data <- data
   
  # calculate the ANGLE of the labels
  number_of_bar <- nrow(label_data)
  angle <-  90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
   
  # calculate the alignment of labels: right or left
  # If I am on the left part of the plot, my labels have currently an angle < -90
  label_data$hjust<-ifelse( angle < -90, 1, 0)
   
  # flip angle BY to make them readable
  label_data$angle<-ifelse(angle < -90, angle+180, angle)
  # ----- ------------------------------------------- ---- #
   
   
  # Start the plot
  p <- ggplot(data, aes(x=as.factor(id), y=value)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
    
    # This add the bars with a blue color
    geom_bar(stat="identity", fill=alpha("skyblue", 0.7)) +
    
    # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
    ylim(-100,120) +
    
    # Custom the theme: no axis title and no cartesian grid
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm")      # Adjust the margin to make in sort labels are not truncated!
    ) +
    
    # This makes the coordinate polar instead of cartesian.
    coord_polar(start = 0) +
    
    # Add the labels, using the label_data dataframe that we have created before
    geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) 
   
    # p
    # 保存
    ggsave(filename = file.path(outdir, 'circular_bar.pdf'), plot = p, width = 13, height = 13)
}
