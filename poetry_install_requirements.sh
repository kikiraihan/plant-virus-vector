#!/bin/bash

while IFS= read -r line || [ -n "$line" ]; do
    poetry add "$line"
done < requirements_slim.txt
