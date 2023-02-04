#!/bin/bash

python parboil clean $1

sources=("$(ls benchmarks/$1/src)")
datasets=("$(ls datasets/$1)")

for source in $sources; do
    echo "python parboil compile $1 $source"
    python parboil compile $1 $source
done

rm results/$1.txt
touch results/$1.txt

for source in $sources; do
    for dataset in $datasets; do
        echo "python parboil run $1 $source $dataset"
        result="$(python parboil run $1 $source $dataset)"
        if [ $? -eq 0 ]
        then
            echo "$2" >> results/$1.txt
            echo "$dataset $source" >> results/$1.txt
            echo "$2" >> results/$1.txt
            echo "$result" >> results/$1.txt
        fi
    done
done

python parboil clean $1
