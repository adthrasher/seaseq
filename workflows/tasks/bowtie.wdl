version 1.0

task bowtie {
    input {
        File fastqfile
        File? metricsfile
        Array[File]+ index_files

        String outputfile = sub(basename(fastqfile),'\_R*[12]\_.+[0-9]\.f.*q\.gz','.sam')
        String outputfile_ = sub(basename(fastqfile),'\.f.*q\.gz','.sam')
        
        Int? read_length = 75
        Int limit_alignments = 2
        Int good_alignments = 2
        Boolean best_alignments = true

        Int memory_gb = 10
        Int max_retries = 1
        Int ncpu = 20
    }
    command <<<
        if [ -f "~{metricsfile}" ]; then
            readlength=$(tail -n 1 ~{metricsfile} | awk '{print $4}');
        else
            readlength=~{read_length}
        fi

        if [ "~{outputfile}" = "~{basename(fastqfile)}" ]; then
            samfile=~{outputfile_}
        else
            samfile=~{outputfile}
        fi

        bowtie \
            -l $readlength \
            -p ~{ncpu} \
            -k ~{good_alignments} \
            -m ~{limit_alignments} \
            ~{true="--best" false="" best_alignments} \
            -S \
            ~{sub(index_files[0], "(\.rev)?\.[0-9]\.ebwt$", "")} \
            ~{fastqfile} \
            > $samfile
    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'quay.io/biocontainers/bowtie:1.2.3--py36hf1ae8f4_2'
        cpu: ncpu
    }
    output {
        Array[File?] samfile = glob("*.sam")
    }
}

task index {
    input {
        File reference

        Int memory_gb = 10
        Int max_retries = 1
        Int ncpu = 20
    }
    command <<<
        bowtie-build --threads ~{ncpu} ~{reference} ~{basename(reference)}
    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'quay.io/biocontainers/bowtie:1.2.3--py36hf1ae8f4_2'
        cpu: ncpu
    }
    output {
        Array[File] bowtie_indexes = glob("~{basename(reference)}*")
    }
}
