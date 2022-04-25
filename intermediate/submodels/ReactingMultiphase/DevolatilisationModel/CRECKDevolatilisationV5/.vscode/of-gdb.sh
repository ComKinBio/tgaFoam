#!/bin/bash
. ~/OpenFOAM/OpenFOAM-v2006/etc/bashrc WM_NCOMPROCS=2; export WM_COMPILE_OPTION=Debug
/usr/bin/gdb "$@"