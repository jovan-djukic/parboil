#!/bin/bash

benchmarks=("$(ls benchmarks)")

for benchmark in $benchmarks; do
    echo "python parboil clean $benchmark"
    python2 parboil clean $benchmark
done
