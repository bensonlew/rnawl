scaleNorm <- function(data,scalenm){
	data = data
	if(scalenm=='Ctr'){
		data<-apply(data, 2, MeanCenter);
	}else if(scalenm=='UV'){
		data<-apply(data, 2, AutoNorm);
	}else if(scalenm =='Par'){
		data<-apply(data, 2, ParetoNorm);
	}else if(scalenm=='RangeNorm'){
		data<-apply(data, 2, RangeNorm);
	}else{
		data <- data
	}
}
  

# normalize to zero mean and unit variance
AutoNorm<-function(x){
  (x - mean(x))/sd(x, na.rm=T);
}

# normalize to zero mean but variance/SE
ParetoNorm<-function(x){
  (x - mean(x))/sqrt(sd(x, na.rm=T));
}

# normalize to zero mean but variance/SE
MeanCenter<-function(x){
  x - mean(x);
}

RangeNorm<-function(x){
  if(max(x) == min(x)){
    x;
  }else{
    (x - mean(x))/(max(x)-min(x));
  }
}

ScaleUnit<-function(x){
    x /sqrt(sum(x^2));
}