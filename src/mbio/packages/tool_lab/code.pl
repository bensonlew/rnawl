#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my ($fasta,$out);

use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	#"fasta:s"=>\$fasta,
	"out:s"=>\$out,
	) or &USAGE;
&USAGE unless ($out);

my (@blk,@line,@new);
#if(-f "$out/$fasta.bulkoutfile"){`rm -rf $out\/$fasta.bulkoutfile`;`rm -rf $out/$fasta.result`;}
#`/mnt/ilustre/users/mdna/Softwares/codonW/codonw $fasta $out/$fasta.result $out/$fasta.bulkoutfile -all_indices -nomenu`;

my $clean1='sed -i \'s/\s\s\s\s\s/\ NAN\ /g\'';
`$clean1 $out/fasta.bulkoutfile`;
my $clean2='sed -i \'s/\s\s\s\s/NAN\ /g\'';
`$clean2  $out/fasta.bulkoutfile`;

open IN , "<$out/fasta.bulkoutfile";
open OUT,">$out/result.txt";
print OUT "AA\tCodon\tSample\n";
while (<IN>) {
	chomp;next if ($_ eq ""||/^$/||/^\d/);
	@blk = split (/\s+/);
	for (my $i=0;$i < 15 ;$i=$i+4) {
		if ($blk[$i] =~m/NAN/ ) {
			$blk[$i] = $line[$i];
		}
		my $w = join ("\t",@blk[$i,$i+1,$i+3]);
		push @new,$w;
	}
	@line = @blk;
	
}
print OUT join ("\n",@new),"\n";
close IN;
close OUT;

open RCMD, ">$out/codon.r";
print RCMD " def.par <- par(no.readonly=TRUE)

library(ggplot2)
library(ggsci)
library(gridExtra)
options(bitmapType='cairo')
newdata<-NULL
data<-read.csv(file=\"$out/result.txt\",header=T, sep=\"\t\")
code<-levels(as.factor(data[,1]))
for (i in 1:length(code)) {
        temp <- data[data\$AA %in% code[i],]
        temp\$text <- c(1)
        temp\$text_label <- 1:nrow(temp)
        temp\$factor <- paste(\"lable\",\"_\",temp\$text_label,sep=\"\")
		newdata <- rbind(newdata,temp)
}

for (i in 1:ncol(data)){
	name <- colnames(data)[i]
	if(grepl(\'AA\',name) || grepl(\'Codon\',name) || grepl(\'text\',name) || grepl(\'factor\',name) || grepl(\'text_label\',name)){
		next
	}else{
		Yaxis <- colnames(data)[i]
		colnames(newdata)[i] <- c(\"Name\")
		p1 <-	ggplot()+
			geom_bar(data=newdata,aes(x=AA,y=Name,fill=factor),stat=\"identity\",show.legend = FALSE,width=0.6)+ylab(paste(\"Frequency \(per thousand\)\"))+
			scale_y_continuous(expand = c(0,0),breaks=seq(0, 7, 0.5),limits=c(0,7)) +
			theme(panel.background=element_blank(),panel.grid.major.y=element_line(color=\"gray\"),
		axis.title.x=element_blank(),axis.title.y=element_text(size=15,face=\"bold\",vjust=-5,margin = unit(c(0,15, 0, 0), \"mm\")),
		axis.text.y=element_text(size=12,face=\"bold\"),axis.text.x=element_text(size=12,face=\"bold\"),axis.ticks.x=element_blank())
		p<-ggplot()+
			geom_bar(data=newdata,aes(x=AA,y=text,fill=factor),stat=\"identity\",colour=\"white\",show.legend = FALSE,width=0.6)+ylab(\"Codon\")+
			scale_y_continuous(trans = \"reverse\",expand = c(0,0),breaks=seq(0, 7, 0.5))+
			geom_text(data=newdata,aes(x=AA,y=text_label,label=Codon),size = 3,vjust=-1)+
			theme(panel.background=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(),
                axis.title.y=element_text(color=\"white\",size=15,face=\"bold\",vjust=-5,margin = unit(c(0,15, 0, 0), \"mm\")),
                axis.ticks.y=element_line(color=\"white\"),axis.text.y=element_text(color=\"white\",size=12,face=\"bold\"))

		p1_png = p1+scale_fill_npg()
		p_png = p+scale_fill_npg()
		codon_plot <- grid.arrange(p1_png, p_png, nrow = 2,heights = c(4,1.5))
		ggsave(codon_plot, file=paste(c(Yaxis),c(\"_\"),c(\"codon.png\"),sep=\"\"), width=10, height=6, device=\"png\")
		colnames(newdata)[i] <- Yaxis
	}
}
";
close RCMD;

#`R --restore --no-save < $out/codon.r`;
#`rm $out/codon.r`;


