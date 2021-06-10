#!/bin/bash

awk '/Space group of crystal/{print FILENAME, $6, $7}' *.castep > spg.txt
