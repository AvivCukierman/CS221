#/bin/bash

for x in 1 
do
  for e in 80
  do
    for i in 10
    do
      echo $i
      echo $e
      echo $x
      python SGD.py -e $e -i $i -x $x
    done
  done
done    
