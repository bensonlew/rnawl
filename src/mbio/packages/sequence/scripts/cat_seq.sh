#! /bin/bash
#cat numbers of seq(fq ,fa or other types) into one file
let a=$#
b=cat
for i in $@
do
	let count=$count+1
	if (("$count"<"$a"))
	then
		b=$b' '$i
	else
		b=$b' >> '$i
	fi
done
eval $b
#echo $b
