#!/bin/sh

idl<<EOF
dir='/sdss4/data1/esheldon/REGRESS/'
matrixfile=dir+'matrix_752_756_g_N1.txt'
evectfile = dir+'evec_752_756_g_N1.txt'
rsumfile=dir+'rsum_752_756_g_N1.txt'

runsolve,matrixfile,evectfile,rsumfile
EOF

idl<<EOF
dir='/sdss4/data1/esheldon/REGRESS/'
matrixfile=dir+'matrix_752_756_r_N1.txt'
evectfile = dir+'evec_752_756_r_N1.txt'
rsumfile=dir+'rsum_752_756_r_N1.txt'

runsolve,matrixfile,evectfile,rsumfile
EOF

idl<<EOF
dir='/sdss4/data1/esheldon/REGRESS/'
matrixfile=dir+'matrix_752_756_i_N1.txt'
evectfile = dir+'evec_752_756_i_N1.txt'
rsumfile=dir+'rsum_752_756_i_N1.txt'

runsolve,matrixfile,evectfile,rsumfile
EOF
