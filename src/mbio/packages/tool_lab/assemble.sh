#!/bin/bash

phrepPhrap=$1
echo "第一个参数为：$phrepPhrap";

if [ $# -ne 1 ]
then
	PAIRDIR=`dirname $2`/`basename $2`
else
	PAIRDIR=`pwd`
fi

if [ ! -d $PAIRDIR ]
then
	echo "输入的目录不存在！" >&2
	exit 1
elif [ ! -w $PAIRDIR ]; then
	echo "输入的目录无权限写入！" >&2
	exit 1
#elif [ "`ls $PAIRDIR/*.ab1 $PAIRDIR/*.seq`" = "" ];then
#	exit $?
fi


echo "读取seq文件..."
: > $PAIRDIR/files.txt
#ls $PAIRDIR/*.seq | while read line; do
ls $PAIRDIR/*.fa | while read line; do
	echo "`basename $line`"
	echo "`basename $line`" >> $PAIRDIR/files.txt
done

#awk -F'-' '{print $1}' $PAIRDIR/files.txt | tee $PAIRDIR/files1.txt
awk -F'-' '{print $1}' $PAIRDIR/files.txt | sort | uniq -c | awk '{if($1==2){print $2"\n"$2 > "'$PAIRDIR'/files1.txt"}else{print $2": seq file != 2"}}'
echo "读取seq文件数大于1的样品数据..."
echo "      数量 样品"
sort $PAIRDIR/files1.txt | uniq -c -d


echo "开始拼接处理"

if [ ! -d $PAIRDIR/contigs ]; then
	mkdir $PAIRDIR/contigs
fi

if [ ! -d $PAIRDIR/tmp ]; then
	mkdir $PAIRDIR/tmp
fi


: > $PAIRDIR/contigs/readme.txt
sort $PAIRDIR/files1.txt | uniq -d | while read line; do
	echo -e "\n\n开始拼接$line样品"
	if [ ! -d $PAIRDIR/tmp/$line ]; then
		mkdir $PAIRDIR/tmp/$line
	fi
	echo -n "$line : " >> $PAIRDIR/contigs/readme.txt
	cd $PAIRDIR/tmp/$line
	#for i in `ls $PAIRDIR/$line*.seq`; do file=`basename $i`; echo ">"$file; cat $i; echo -e "\n"; done > $line.fasta.screen
	for i in `ls $PAIRDIR/$line*.fa`; do file=`basename $i`; echo ">"$file; cat $i; echo -e "\n"; done > $line.fasta.screen
	$phrepPhrap/bin/phrap $line.fasta.screen -new_ace
	$phrepPhrap/bin/tagRepeats.perl $line.fasta.screen.ace
	$phrepPhrap/bin/cross_match $line.contigs $phrepPhrap/lib/screenLibs/repeats.fasta -tags -minmatch 10
	if [ -f $PAIRDIR/tmp/$line/$line.contigs -a -s $PAIRDIR/tmp/$line/$line.contigs ]; then
		cp $PAIRDIR/tmp/$line/$line.contigs $PAIRDIR/contigs/
		echo  " 成功 " >> $PAIRDIR/contigs/readme.txt
	else
		echo  " 未成功 " >> $PAIRDIR/contigs/readme.txt
	fi
done
echo "contig文件为Linux下文本文件，windows下可以使用写字板打开，推荐安装使用gedit windows版本：http://projects.gnome.org/gedit/" >> $PAIRDIR/contigs/readme.txt

#cd  $PAIRDIR/contigs
#for shname in `ls *.contigs`; do name=`echo "$shname" | awk -F. '{print $1}'` ; more $shname|sed "s/>Contig1/>$name/g" >$name.seq; done
#cat *.seq |sed  '/^$/d' >all.contigs
#cd ../

echo "正在删除临时文件..."
rm -rf $PAIRDIR/tmp
rm -f $PAIRDIR/files.txt
rm -f $PAIRDIR/files1.txt
echo "拼接全部完成，拼接结果在$PAIRDIR/contigs目录下"
exit 0
