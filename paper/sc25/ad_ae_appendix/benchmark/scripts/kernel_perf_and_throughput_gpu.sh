#!/bin/bash

for i in $@; do
  ng=$(echo $i | sed -e s/ng/\ / | cut -f2 -d' ' | cut -f1 -d'_')
  tot=$(grep Total $i | cut -f5 -d' ')
  kernel=$(grep Kernel $i | awk '{print $2}')
  nint=$(grep Total $i | awk '{print $2}') 
  echo $ng $nint $tot $kernel | awk '{print $1,$2,$3,24*$2/$4*1e-9}'
done | sort -n -k1
