#!/bin/sh


#lines=`grep "~~~~ Cut 1 ~~~~" output.log -n | cut -f1 -d:`
lines=( $(grep "~~~~ Cut 1 ~~~~" -n output.log | cut -f1 -d:) )

flines=`wc -l output.log | cut -f1 -d " "`
lout=( $(grep "~~~~ Cut 1 ~~~~" -n output.log | tail -n+2 | cut -f1 -d:) $flines)

echo "Looping"
for ((i=0;i<${#lines[@]};i++));do
	file=`sed -n "$(( ${lines[$i]}-2 )),$(( ${lout[$i]}-2 ))p" output.log | sed -e 's/\r/\n/g'`
	newline=`echo "$file" | grep -n "~~~~ Cut 1 ~~~~" | cut -f1 -d:`
	echo "$file" | tail -n+$(( $newline-3 )) > event$i
done
