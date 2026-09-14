#!bin/bash

for file in "ree1.txt" "ree2.txt"; do
val=$(awk '!/^#/ {print $2}' $file)
echo $val > "p_$file"
done
