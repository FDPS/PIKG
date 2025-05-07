#!/bin/bash

for i in $@;do
  grep Total $i > /dev/null
  if [ $? == 0 ]; then
    ninteract=$(grep Total $i | awk '{print $2}')
  else
    ninteract=0
  fi
  n=$(echo $i | sed -e s/ng/\ /g | awk '{print $2}' | cut -f1 -d'_')
  l=$(echo $i | sed -e s/_l/\ /g | awk '{print $2}' | cut -f1 -d'.')
  echo $n $l $ninteract
done | sort -n -k1 | sort -n -k2 --stable
echo
echo
