#!/usr/bin/env bash


if [ $1 = '' ]
then
	echo "Please supply the dot file name (along with path) to convert. Program ending."

else
	echo "converting to a temporary tmp file"
	cat $1 | sed 's/L:\([^ ]\+\), R:/\1-/g'   | sed 's/label="[0-9][^"]*"/label=""/g' | awk -F'=' '{ if ($1 == "penwidth") {print $1 "=" ($2 ^ 6) ","} else {print $0 }}'   | tr '\n' ' '   | sed "s/;/\n/g"  > tmp

	echo "converting to svg file"
	cat tmp   | dot -Kneato -n -y -Tsvg -Efontname="Arial" -Nlabel="" -Nwidth=1 -Nheight=1 -Epenwidth=20 -Nstyle="filled"> edge_graph.svg #-Npenwidth=20 

	echo "converting to a pdf file"
	cat tmp   | dot -Kneato -n -y -Tpdf -Efontname="Arial" -Nlabel="" -Nwidth=1 -Nheight=1 -Epenwidth=20 -Nstyle="filled"> edge_graph.pdf #-Npenwidth=20 

	echo "Deleting the temporary temp file"
	rm tmp
  

fi


