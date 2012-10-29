#!/bin/bash
rm generated/*
python runpytests.py
make
