#!/usr/bin/env bash

for i in *.fastq; do mv "$i" "${i/lane1-/}"; done
