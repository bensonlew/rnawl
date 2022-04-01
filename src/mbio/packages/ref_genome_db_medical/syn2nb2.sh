#!/bin/bash
FILE=${1}

FILE=`readlink -f $FILE`
FILETO=$FILE
DIR=$FILE
if [ -d "$FILE" ]; then
    FILETO=`echo $FILE | sed 's/\/$//g' |sed 's/[^\/]*$//'`
    DIR=$FILETO
else
    DIR=`echo $FILE | sed 's/[^\/]*$//'`
fi
FILETO_SG=`echo $FILETO|sed 's/ilustre/lustre/;s/sanger-dev/sanger/;s/biocluster/sanger_bioinfo/'`
DIR_SG=`echo $DIR|sed 's/ilustre/lustre/;s/sanger-dev/sanger/;s/biocluster/sanger_bioinfo/'`

FILETO_ISG=`echo $FILETO|sed 's/sanger-dev/isanger/;s/biocluster/sanger_bioinfo/'`
DIR_ISG=`echo $DIR|sed 's/sanger-dev/isanger/;s/biocluster/sanger_bioinfo/'`

# echo $DIR
# create dir
HOME='/mnt/ilustre/users/sanger-dev'
ssh -i $HOME/sg-users/liubinxu/script/.key/nb_id_rsa sanger@10.2.0.110 '
if [ -d '$DIR_SG' ];then
echo "文件夹存在"
else
echo "文件夹不存在, 创建"
mkdir -p '$DIR_SG'
fi
'
# echo $FILETO
scp -r -i $HOME/sg-users/liubinxu/script/.key/nb_id_rsa $FILE sanger@10.2.0.110:$DIR_SG

# create dir
ssh -i $HOME/sg-users/liubinxu/script/.key/nb2_rsa isanger@10.2.0.115 '
if [ -d '$DIR_ISG' ];then
echo "文件夹存在"
else
echo "文件夹不存在, 创建"
mkdir -p '$DIR_ISG'
fi
'

# echo $FILETO
scp -r -i $HOME/sg-users/liubinxu/script/.key/nb2_rsa $FILE isanger@10.2.0.115:$DIR_ISG
