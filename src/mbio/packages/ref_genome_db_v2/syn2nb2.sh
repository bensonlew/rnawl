#!/bin/bash -x
FILE=${1}
task_status=${2}
if [ $task_status == "tsg_finish" ]; then
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

    FILETO_ISG=`echo $FILETO|sed 's/lustre/ilustre/g;s/sanger-dev/isanger/;s/biocluster/sanger_bioinfo/'`
    DIR_ISG=`echo $DIR|sed 's/lustre/ilustre/g;s/sanger-dev/isanger/;s/biocluster/sanger_bioinfo/'`

    # echo $DIR
    # create dir
    HOME='/mnt/lustre/users/sanger-dev'
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
elif [ $task_status == "sanger_finish" ]; then
    if [ -d "$FILE" ] ; then
        echo "sanger for sanger"
        FILE=`readlink -f $FILE`
        FILETO=$FILE
        DIR=$FILE
        FILETO=`echo $FILE | sed 's/\/$//g' |sed 's/[^\/]*$//'`
        DIR=$FILETO
        FILETO_ISG=`echo $FILETO|sed 's/lustre/ilustre/;s/sanger/isanger/'`
        DIR_ISG=`echo $DIR|sed 's/lustre/ilustre/;s/sanger/isanger/'`
        FILETO_TSG=`echo $FILETO|sed 's/sanger/sanger-dev/'`
        DIR_TSG=`echo $DIR|sed 's/sanger/sanger-dev/'`

        HOME='/mnt/lustre/users/sanger'
        ssh -i $HOME/app/database/Genome_DB_finish/.key/nb2_rsa isanger@10.2.0.115 '
        if [ -d '$DIR_ISG' ];then
            echo "文件夹存在"
        else
            echo "文件夹不存在, 创建"
            mkdir -p '$DIR_ISG'
        fi
        '
        # echo $FILETO
        scp -r -i $HOME/app/database/Genome_DB_finish/.key/nb2_rsa $FILE isanger@10.2.0.115:$DIR_ISG
        ssh -i $HOME/app/database/Genome_DB_finish/.key/new_tsg_id_rsa sanger-dev@10.2.3.173 '
        if [ -d '$DIR_TSG' ];then
            echo "文件夹存在"
        else
            echo "文件夹不存在, 创建"
            mkdir -p '$DIR_TSG'
        fi
        '
        # echo $FILETO
        scp -r -i $HOME/app/database/Genome_DB_finish/.key/new_tsg_id_rsa $FILE sanger-dev@10.2.3.173:$DIR_TSG
    else
        echo "sanger for isanger"
        FILETO=$FILE
        DIR=`echo $FILE | sed 's/[^\/]*$//'`
        FILETO_ISG=`echo $FILETO|sed 's/lustre/ilustre/;s/sanger/isanger/'`
        DIR_ISG=`echo $DIR|sed 's/lustre/ilustre/;s/sanger/isanger/'`
        FILETO_TSG=`echo $FILETO|sed 's/sanger/sanger-dev/'`
        DIR_TSG=`echo $DIR|sed 's/sanger/sanger-dev/'`
        HOME='/mnt/ilustre/users/isanger'
        if [ -d '$DIR_ISG' ];then
            echo "文件夹存在"
        else
            echo "文件夹不存在, 创建"
            mkdir -p '$DIR_ISG'
        fi
        scp -r -i $HOME/app/database/Genome_DB_finish/.key/nb_id_rsa sanger@10.2.0.110:$FILE $DIR_ISG
        ssh -i $HOME/app/database/Genome_DB_finish/.key/new_tsg_id_rsa sanger-dev@10.2.3.173 '
        if [ -d '$DIR_TSG' ];then
            echo "文件夹存在"
        else
            echo "文件夹不存在, 创建"
            mkdir -p '$DIR_TSG'
        fi
        '
        # echo $FILETO
        scp -r -i $HOME/app/database/Genome_DB_finish/.key/new_tsg_id_rsa $FILETO_ISG sanger-dev@10.2.3.173:$DIR_TSG
    fi
else
    if [ -d "$FILE" ] ; then
        echo "isanger for isanger"
        FILE=`readlink -f $FILE`
        FILETO=$FILE
        DIR=$FILE
        FILETO=`echo $FILE | sed 's/\/$//g' |sed 's/[^\/]*$//'`
        DIR=$FILETO
        FILETO_SG=`echo $FILETO|sed 's/ilustre/lustre/;s/isanger/sanger/'`
        DIR_SG=`echo $DIR|sed 's/ilustre/lustre/;s/isanger/sanger/'`
        FILETO_TSG=`echo $FILETO|sed 's/ilustre/lustre/;s/isanger/sanger-dev/'`
        DIR_TSG=`echo $DIR|sed 's/ilustre/lustre/;s/isanger/sanger-dev/'`

        HOME='/mnt/ilustre/users/isanger'
        ssh -i $HOME/app/database/Genome_DB_finish/.key/nb_id_rsa sanger@10.2.0.110 '
        if [ -d '$DIR_SG' ];then
            echo "文件夹存在"
        else
            echo "文件夹不存在, 创建"
            mkdir -p '$DIR_SG'
        fi
        '
        # echo $FILETO
        scp -r -i $HOME/app/database/Genome_DB_finish/.key/nb_id_rsa $FILE sanger@10.2.0.110:$DIR_SG
        ssh -i $HOME/app/database/Genome_DB_finish/.key/new_tsg_id_rsa sanger-dev@10.2.3.173 '
        if [ -d '$DIR_TSG' ];then
            echo "文件夹存在"
        else
            echo "文件夹不存在, 创建"
            mkdir -p '$DIR_TSG'
        fi
        '
        # echo $FILETO
        scp -r -i $HOME/app/database/Genome_DB_finish/.key/new_tsg_id_rsa $FILE sanger-dev@10.2.3.173:$DIR_TSG
    else
        echo "isanger for sanger"
        FILETO=$FILE
        DIR=`echo $FILE | sed 's/[^\/]*$//'`
        FILETO_SG=`echo $FILETO|sed 's/ilustre/lustre/;s/isanger/sanger/'`
        DIR_SG=`echo $DIR|sed 's/ilustre/lustre/;s/isanger/sanger/'`
        FILETO_TSG=`echo $FILETO|sed 's/ilustre/lustre/;s/isanger/sanger-dev/'`
        DIR_TSG=`echo $DIR|sed 's/ilustre/lustre/;s/isanger/sanger-dev/'`
        HOME='/mnt/lustre/users/sanger'
        if [ -d '$DIR_SG' ];then
            echo "文件夹存在"
        else
            echo "文件夹不存在, 创建"
            mkdir -p '$DIR_SG'
        fi
        scp -r -i $HOME/app/database/Genome_DB_finish/.key/nb2_rsa isanger@10.2.0.115:$FILE $DIR_SG
        ssh -i $HOME/app/database/Genome_DB_finish/.key/new_tsg_id_rsa sanger-dev@10.2.3.173 '
        if [ -d '$DIR_TSG' ];then
            echo "文件夹存在"
        else
            echo "文件夹不存在, 创建"
            mkdir -p '$DIR_TSG'
        fi
        '
        # echo $FILETO
        scp -r -i $HOME/app/database/Genome_DB_finish/.key/new_tsg_id_rsa $FILETO_SG sanger-dev@10.2.3.173:$DIR_TSG
    fi
fi