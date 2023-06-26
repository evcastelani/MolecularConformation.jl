#!/bin/sh

echo -n "list_of_problems = [" > list.txt
for file in *.xyz; do			
	echo -n '"' >> list.txt
	echo -n ${file%.xyz} >> list.txt
	echo -n '", ' >> list.txt
done
echo "]" >> list.txt