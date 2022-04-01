#! /bin/bash
for i in $@
do
	let count=$count+1
	if [ $count -lt $# ] && [ $count -gt 2 ]
	then
		$1 -jar $2 -i $i -q 30 >> ${@:$#:1}
	fi
done