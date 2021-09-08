#!/bin/bash

for i in {1..10}; do
  for j in {1..10}; do
    echo $((i*j))
  done
done
