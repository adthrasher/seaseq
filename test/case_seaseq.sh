#!/bin/bash
## Test WDL script using samplefiles provided.

# CROMWELL on local machine
local_cromwell="/Users/madetunj/Downloads/cromwell-61.jar"

# CROMWELL on St. Jude hpc user directory
sjhpc_cromwell="/rgs01/project_space/abrahgrp/Software_Dev_Sandbox/common/madetunj/software/cromwell-52.jar"
lsf_config="/home/madetunj/.commands/lsf.conf"
input="inputs-case.json"
option="options-case.json"

# STD OUT and ERR files
logout="wdlseaseq-case_out"
logerr="wdlseaseq-case_err"

# check if cromwell app exists
if [ -f "$local_cromwell" ]; then
    java -jar $local_cromwell \
        run ../seaseq-case.wdl \
        -i $input \
        -o $option \
        1>$logout 2>$logerr

elif [ -f "$sjhpc_cromwell" ]; then 
    #script syntax
    wdlscript="java -Dconfig.file=$lsf_config \
        -jar $sjhpc_cromwell \
        run ../seaseq-case.wdl \
        --inputs $input \
        --options $option"
    bsub -P watcher -q compbio \
        -R "rusage[mem=10000]" \
        -J wdlseaseq-case \
        -o $logout \
        -e $logerr \
        -N $wdlscript

else
    echo "cromwell executable jar doesN'T exist"

fi

