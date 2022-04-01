#plot KEGG pic
cols=c("#313695", "#4575B4", "#74ADD1", "#ABD9E9" ,"#E0F3F8" ,"#FEE090", "#FDAE61", "#F46D43" ,"#D73027","#A50026")

parse_XMLfile<-function(pathway_id,org="ko",database_dir=getwd()) {
  if (pathway_id=="01100" | pathway_id=="01110" | pathway_id=="01120" | pathway_id=="01130"| pathway_id=="01200"| pathway_id=="01210"| pathway_id=="01212"| pathway_id=="01230"| pathway_id=="01220") {
    print(paste("Skip global maps:",pathway_id,sep=""))
    return(NULL)
  }
  require(XML)
  # inter1<-xmlTreeParse(paste(database_dir,"/",org,pathway_id,".xml",sep=""),useInternalNodes=TRUE)
  inter1<-xmlTreeParse(paste(xml_path, sep=""),useInternalNodes=TRUE)
  inter2<-getNodeSet(inter1,"//entry")
  inter3<-lapply(inter2,  function(xxx) xmlGetAttr(xxx,  "name"))
  inter4<-lapply(inter2,  function(xxx) xmlGetAttr(xxx,  "type"))
  inter5<-sapply(inter2,  function(xxx) getNodeSet(xxx,".//graphics"))
  inter_graphic_type<-sapply(inter5,  function(xxx) xmlGetAttr(xxx,  "type"))
  inter6<-sapply(inter5,  function(xxx) xmlGetAttr(xxx,  "x"))
  inter7<-sapply(inter5,  function(xxx) xmlGetAttr(xxx,  "y"))
  inter8<-sapply(inter5,  function(xxx) xmlGetAttr(xxx,  "width"))
  inter9<-sapply(inter5,  function(xxx) xmlGetAttr(xxx,  "height"))
  inter10<-sapply(inter5,  function(xxx) xmlGetAttr(xxx,  "name"))
  result<-NULL
  for (i in 1:length(inter4)) {
    if ((inter4[[i]]=="ortholog" | inter4[[i]]=="compound") & inter_graphic_type[[i]]!="line") {
      temp<-strsplit(inter3[[i]]," ")[[1]]
      name<-strsplit(inter10[[i]],",")[[1]][1]
      name<-gsub('\\.\\.\\.',"",name)
      for (n in 1:length(temp)) {
        result<-rbind(result,c(temp[n],inter6[[i]],inter7[[i]],inter8[[i]],inter9[[i]],name))
      }
    }
  }
  return(result)
}


plot_polygon<-function(choose_xml,i,i_max,col,height,err_x_location,err_y_location,magnify,border,...) {
	x<-choose_xml$x
	y<-choose_xml$y
	w<-choose_xml$width+2
	h<-choose_xml$height+1
	polygon(c(
		x-w/2*magnify+w/i_max*(i-1)*magnify+err_x_location,
		x-w/2*magnify+w/i_max*i*magnify+err_x_location,
		x-w/2*magnify+w/i_max*i*magnify+err_x_location,
		x-w/2*magnify+w/i_max*(i-1)*magnify+err_x_location),
			c(height-y-h/2*magnify-err_y_location,height-y-h/2*magnify-err_y_location,height-y+h/2*magnify-err_y_location,height-y+h/2*magnify-err_y_location)
			,col=col[i],border=border[i],lwd = 2,...)
}

plot_segments<-function(choose_xml,i,i_max,col,height,err_x_location,err_y_location,magnify,border,...) {
  x<-choose_xml$x
  y<-choose_xml$y
  w<-choose_xml$width+2
  h<-choose_xml$height+1
	if(i==1){
		x0=c(x-w/2*magnify+w/i_max*i*magnify+err_x_location,x-w/2*magnify+w/i_max*(i-1)*magnify+err_x_location,x-w/2*magnify+w/i_max*(i-1)*magnify+err_x_location)
		x1=c(x-w/2*magnify+w/i_max*(i-1)*magnify+err_x_location,x-w/2*magnify+w/i_max*(i-1)*magnify+err_x_location,x-w/2*magnify+w/i_max*i*magnify+err_x_location)
		y0=c(height-y+h/2*magnify-err_y_location,height-y+h/2*magnify-err_y_location,height-y-h/2*magnify-err_y_location)
		y1=c(height-y+h/2*magnify-err_y_location,height-y-h/2*magnify-err_y_location,height-y-h/2*magnify-err_y_location)
	}else if(i==i_max){
		x0=c(x-w/2*magnify+w/i_max*(i-1)*magnify+err_x_location,x-w/2*magnify+w/i_max*i*magnify+err_x_location,x-w/2*magnify+w/i_max*i*magnify+err_x_location)
		x1=c(x-w/2*magnify+w/i_max*i*magnify+err_x_location,x-w/2*magnify+w/i_max*i*magnify+err_x_location,x-w/2*magnify+w/i_max*(i-1)*magnify+err_x_location)
		y0=c(height-y+h/2*magnify-err_y_location,height-y+h/2*magnify-err_y_location,height-y-h/2*magnify-err_y_location)
		y1=c(height-y+h/2*magnify-err_y_location,height-y-h/2*magnify-err_y_location,height-y-h/2*magnify-err_y_location)
	}else{
		x0=c(x-w/2*magnify+w/i_max*i*magnify+err_x_location,x-w/2*magnify+w/i_max*i*magnify+err_x_location)
		x1=c(x-w/2*magnify+w/i_max*(i-1)*magnify+err_x_location,x-w/2*magnify+w/i_max*(i-1)*magnify+err_x_location)
		y0=c(height-y+h/2*magnify-err_y_location,height-y-h/2*magnify-err_y_location)
		y1=c(height-y+h/2*magnify-err_y_location,height-y-h/2*magnify-err_y_location)
	}
	segments(x0,y0,x1,y1,col=border[i],lwd = 2)
}

# plot_KO <-function(xml_KO,bg_col,fg_col,height,err_x_location,err_y_location){
plot_KO <-function(xml_KO,height,err_x_location,err_y_location){
  # if(!is.na(bg_col[1])){bg_col = paste(bg_col,"99",sep="")}
  # i_max<-length(bg_col)
  i_max = 1
  #if(i_max==1){
  #  for (i in 1:i_max) {
      # plot_polygon(xml_KO,i=i,i_max=i_max,col=bg_col,height=height,err_x_location=err_x_location,err_y_location=err_y_location,magnify=1,border = fg_col)
      plot_polygon(xml_KO,i=i,i_max=i_max,col=NA,height=height,err_x_location=err_x_location,err_y_location=err_y_location,magnify=1,border = NA)
  #  }
  #}else{
  #  for (i in 1:i_max) {
  #    plot_polygon(xml_KO,i=i,i_max=i_max,col=bg_col,height=height,err_x_location=err_x_location,err_y_location=err_y_location,magnify=1,border = NA)
  #  }
  #  i_max<-length(fg_col)
  #  for (i in 1:i_max) {
  #    plot_segments(xml_KO,i=i,i_max=i_max,col=NA,border=fg_col,height=height,err_x_location=err_x_location,err_y_location=err_y_location,magnify=1)
  #  }
  #}
}



args <- commandArgs(T)
pathway_name = args[1]
kos = args[2]
pdf_path = args[3]
# database_dir = args[4]
xml_path = args[4]
png_path = args[5]


path_KOs <- read.table(kos,head = F,comment.char = "")
rownames(path_KOs)=path_KOs[,1]


# xml <-parse_XMLfile(gsub("[a-zA-Z]+","",pathway_name),database_dir=database_dir)
xml <-parse_XMLfile(xml_path)
if(is.null(xml)){
  stop("pathway xml has no rectangle")
}

xml[,1]<-gsub("ko:","",xml[,1])

inter_KOs <- intersect(xml[,1],as.character(path_KOs[,1]))
if(length(inter_KOs)==0){
  stop("KOs and pathway xml has no intersections")
}

choose_xml <- xml[which(xml[,1] %in% as.character(inter_KOs)),c(1:5)]

choose_xml<-as.data.frame(choose_xml,stringsAsFactors=FALSE,drop=FALSE)
if(ncol(choose_xml) == 1){
  choose_xml <- as.data.frame(t(choose_xml),stringsAsFactors=FALSE,drop=FALSE)
}

colnames(choose_xml)<-c("name","x","y","width","height")
choose_xml<-transform(choose_xml, x = as.numeric(x), y = as.numeric(y),width = as.numeric(width),height = as.numeric(height))
choose_xml<-choose_xml[!duplicated(paste(choose_xml$x,choose_xml$y,sep="-")),]

require(png)
# img  <-  readPNG(paste(database_dir,"/",gsub("[a-zA-Z]+","map",pathway_name),".png",sep=""))
img  <-  readPNG(paste(png_path,sep=""))
width<-ncol(img)
height<-nrow(img)
err_x_location<-0
err_y_location<-0



png(paste(pdf_path,sep=""),width=width,height=height,type="cairo",units = "px")
par(yaxs="i")
par(xaxs="i")
par(mar=c(0,0,0,0))
plot(c(0,width),c(0,height),  type='n',xlab="",ylab="",axes=FALSE)
rasterImage(img,  0,  0,  width,  height,interpolate=F)


for(i in 1:length(choose_xml[,1])){
  #bg_col = unlist(strsplit(as.character(path_KOs[choose_xml[i,]$name,]$bg),","))
  #fg_col = unlist(strsplit(as.character(path_KOs[choose_xml[i,]$name,]$fg),","))
  #plot_KO(choose_xml[i,],bg_col,fg_col,height,err_x_location,err_y_location)
  plot_KO(choose_xml[i,],height,err_x_location,err_y_location)
}

dev.off()
