#!/usr/bin/env bash

gfortran -fpic -shared solve_tableau.f95 dgesv.f -o solve_tableau.so
