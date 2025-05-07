#!/bin/bash
set -x
for i in $@; do
  ng=$(echo $i | sed -e s/ng/\ / | cut -f2 -d' ' | cut -f1 -d'_')
  tot=$(grep Total $i | cut -f5 -d' ')
  nint=$(grep sec $i | head -n1 | awk '{print $2}') 
  echo $ng $nint $tot | awk '{print $1,$2,$3}'
done | sort -n -k1
