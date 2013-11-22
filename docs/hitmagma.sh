#!/bin/sh

# So long as the server requests take > 90 seconds, you can run this script all night without getting kicked off their site from a single IP.

x=9

while [ $x -le "63" ]
do

y=36

while [ $y -le "63" ]
do

myvar=$(sed -e "s/FREEVAR/$x/g" -e "s/VARTWO/$y/g" magma_req.txt)

curl --referer http://magma.maths.usyd.edu.au/calc/ --data-urlencode "input=$myvar" http://magma.maths.usyd.edu.au/xml/calculator.xml

y=$(($y+6))

done 

x=$(($x+1))

done 

