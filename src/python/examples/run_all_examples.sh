#!/usr/bin/env bash
export FRIGUS_DATADIR_ROOT=../../../data
export PYTHONPATH=..:${PYTHONPATH} 
ls -1 *.py | xargs -L1 python

