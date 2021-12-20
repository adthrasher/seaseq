version 1.0
import "workflows/tasks/fastqc.wdl"
import "workflows/tasks/bedtools.wdl"
import "workflows/tasks/bowtie.wdl"
import "workflows/tasks/samtools.wdl"
import "workflows/tasks/macs.wdl"
import "workflows/workflows/bamtogff.wdl"
import "workflows/tasks/sicer.wdl"
import "workflows/workflows/motifs.wdl"
import "workflows/tasks/rose.wdl"
import "workflows/tasks/util.wdl"
import "workflows/workflows/visualization.wdl" as viz
import "workflows/workflows/mapping.wdl"
import "workflows/tasks/runspp.wdl"
import "workflows/tasks/sortbed.wdl"
import "workflows/tasks/sratoolkit.wdl" as sra

workflow peaseq {
    String pipeline_ver = 'v2.0.0'

    meta {
        title: 'PEAseq Analysis'
        summary: 'Paired-End Antibody Sequencing (PEAseq) Pipeline'
        description: 'A comprehensive automated computational pipeline for all ChIP-Seq/CUT&RUN data analysis.'
        version: '1.0.0'
        details: {
            citation: 'pending',
            contactEmail: 'modupeore.adetunji@stjude.org',
            contactOrg: "St Jude Children's Research Hospital",
            contactUrl: "",
            upstreamLicenses: "MIT",
            upstreamUrl: 'https://github.com/stjude/seaseq',
            whatsNew: [
                {
                    version: "1.0",
                    changes: "Initial release"
                }
            ]
        }
        parameter_group: {
            reference_genome: {
                title: 'Reference genome',
                description: 'Genome specific files. e.g. reference FASTA, GTF, blacklist, motif databases, FASTA index, bowtie index .',
                help: 'Input reference genome files as defined. If some genome data are missing then analyses using such data will be skipped.'
            },
            input_genomic_data: {
                title: 'Input FASTQ data',
                description: 'Genomic input files for experiment.',
                help: 'Input one or more sample data and/or SRA identifiers.'
            },
            analysis_parameter: {
                title: 'Analysis parameter',
                description: 'Analysis settings needed for experiment.',
                help: 'Analysis settings; such output analysis file name.'
            }
        }
    }
    input {
        # group: reference_genome
        File reference
        File? blacklist
        File gtf
        Array[File]? bowtie_index
        Array[File]? motif_databases

        # group: input_genomic_data
        Array[String]? sample_sraid
        Array[File]? sample_R1_fastq
        Array[File]? sample_R2_fastq

        # group: analysis_parameter
        String? results_name
        Boolean run_motifs=true

    }

    parameter_meta {
        reference: {
            description: 'Reference FASTA file',
            group: 'reference_genome',
            patterns: ["*.fa", "*.fasta", "*.fa.gz", "*.fasta.gz"]
        }
        blacklist: {
            description: 'Blacklist file in BED format',
            group: 'reference_genome',
            help: 'If defined, blacklist regions listed are excluded after reference alignment.',
            patterns: ["*.bed", "*.bed.gz"]
        }
        gtf: {
            description: 'gene annotation file (.gtf)',
            group: 'reference_genome',
            help: 'Input gene annotation file from RefSeq or GENCODE (.gtf).',
            patterns: ["*.gtf", "*.gtf.gz", "*.gff", "*.gff.gz", "*.gff3", "*.gff3.gz"]
        }
        bowtie_index: {
            description: 'bowtie v1 index files (*.ebwt)',
            group: 'reference_genome',
            help: 'If not defined, bowtie v1 index files are generated, will take a longer compute time.',
            patterns: ["*.ebwt"]
        }
        motif_databases: {
            description: 'One or more of the MEME suite motif databases (*.meme)',
            group: 'reference_genome',
            help: 'Input one or more motif databases available from the MEME suite (https://meme-suite.org/meme/db/motifs).',
            patterns: ["*.meme"]
        }
        sample_sraid: {
            description: 'One or more sample SRA (Sequence Read Archive) run identifiers',
            group: 'input_genomic_data',
            help: 'Input publicly available FASTQs (SRRs). Multiple SRRs are separated by commas (,).',
            example: 'SRR12345678'
        }
        sample_R1_fastq: {
            description: 'One or more sample R1 FASTQs',
            group: 'input_genomic_data',
            help: 'Upload zipped FASTQ files.',
            patterns: ["*.fq.gz", "*.fastq.gz"]
        }
        sample_R2_fastq: {
            description: 'One or more sample R2 FASTQs',
            group: 'input_genomic_data',
            help: 'Upload zipped FASTQ files.',
            patterns: ["*.fq.gz", "*.fastq.gz"]
        }
        results_name: {
            description: 'Experiment results custom name',
            group: 'analysis_parameter',
            help: 'Input preferred analysis results name (recommended if multiple FASTQs are provided).',
            example: 'AllMerge_mapped'
        }
        run_motifs: {
            description: 'Perform Motif Analysis',
            group: 'analysis_parameter',
            help: 'Setting this means Motif Discovery and Enrichment analysis will be performed.',
            example: true
        }
    }

### ---------------------------------------- ###
### ------------ S E C T I O N 1 ----------- ###
### ------ pre-process analysis files ------ ###
### ---------------------------------------- ###

    # Process SRRs
    if ( defined(sample_sraid) ) {
        # Download sample file(s) from SRA database
        # outputs:
        #    fastqdump.fastqfile : downloaded sample files in fastq.gz format
        Array[String] string_sra = [1] #buffer to allow for sra_id optionality

        Array[String] s_sraid = select_first([sample_sraid, string_sra])
        scatter (eachsra in s_sraid) {
            call sra.fastqdump {
                input :
                    sra_id=eachsra,
                    cloud=false
            }
            File R1end = select_first([fastqdump.R1end, string_sra[0]])
            File R2end = select_first([fastqdump.R2end, string_sra[0]])
        } # end scatter each sra

        Array[File] sample_R1_srafile_ = R1end 
        Array[File] sample_R2_srafile_ = R2end 
    } # end if sample_sraid

    # Generating INDEX files
    #1. Bowtie INDEX files if not provided
    if ( !defined(bowtie_index) ) {
        # create bowtie index when not provided
        call bowtie.index as bowtie_idx {
            input :
                reference=reference
        }
    }
    #2. Make sure indexes are six else build indexes
    if ( defined(bowtie_index) ) {
        # check total number of bowtie indexes provided
        Array[String] string_bowtie_index = [1] #buffer to allow for bowtie_index optionality
        Array[File] int_bowtie_index = select_first([bowtie_index, string_bowtie_index])
        if ( length(int_bowtie_index) != 6 ) {
            # create bowtie index if 6 index files aren't provided
            call bowtie.index as bowtie_idx_2 {
                input :
                    reference=reference
            }
        }
    }
    Array[File] bowtie_index_ = select_first([bowtie_idx_2.bowtie_indexes, bowtie_idx.bowtie_indexes, bowtie_index])

    # FASTA faidx and chromsizes and effective genome size
    call samtools.faidx as samtools_faidx {
        # create FASTA index and chrom sizes files
        input :
            reference=reference
    }
    call util.effective_genome_size as egs {
        # effective genome size for FASTA
        input :
            reference=reference
    }

    # Process FASTQs
    if ( defined(sample_R1_fastq) ) {
        Array[String] string_fastq = [1] #buffer to allow for fastq optionality
        Array[File] s_R1_fastq = select_first([sample_R1_fastq, string_fastq])
        Array[File] sample_R1_fastqfile_ = s_R1_fastq
        Array[File] s_R2_fastq = select_first([sample_R2_fastq, string_fastq])
        Array[File] sample_R2_fastqfile_ = s_R2_fastq
    } 
    # collate all fastqfiles
    Array[File] sample_R1 = flatten(select_all([sample_R1_srafile_, sample_R1_fastqfile_]))
    Array[File] sample_R2 = flatten(select_all([sample_R2_srafile_,sample_R2_fastqfile_]))

    # transpose to paired-end tuples
    Array[Array[File]] sample_fastqfiles = transpose([sample_R1, sample_R2])
    
### ------------------------------------------------- ###
### ---------------- S E C T I O N 2 ---------------- ###
### ---- A: analysis if multiple FASTQs provided ---- ###
### ------------------------------------------------- ###

    # if multiple fastqfiles are provided
    Boolean multi_fastq = if length(sample_fastqfiles) > 2 then true else false
    Boolean one_fastq = if length(sample_fastqfiles) == 2 then true else false

    if ( multi_fastq ) {
        scatter (eachfastq in fastqfiles) {
            # Execute analysis on each fastq file provided
            # Analysis executed:
            #   FastQC
            #   FASTQ read length distribution
            #   Reference Alignment using Bowtie (-k2 -m2)
            #   Convert SAM to BAM
            #   FastQC on BAM files
            #   Remove Blacklists (if provided)
            #   Remove read duplicates
            #   Summary statistics on FASTQs
            #   Combine html files into one for easy viewing
            
            call fastqc.fastqc as indv_fastqc {
                input :
                    inputfile=eachfastq,
                    default_location='SAMPLE/' + sub(basename(eachfastq),'\.fastq\.gz|\.fq\.gz','') + '/QC/FastQC'
            }

            call util.basicfastqstats as indv_bfs {
                input :
                    fastqfile=eachfastq,
                    default_location='SAMPLE/' + sub(basename(eachfastq),'\.fastq\.gz|\.fq\.gz','') + '/QC/SummaryStats'
            }

            call mapping.mapping as indv_mapping {
                input :
                    fastqfile=eachfastq,
                    index_files=bowtie_index_,
                    metricsfile=indv_bfs.metrics_out,
                    blacklist=blacklist,
                    default_location='SAMPLE/' + sub(basename(eachfastq),'\.fastq\.gz|\.fq\.gz','') + '/BAM_files'
            }

            call fastqc.fastqc as indv_bamfqc {
                input :
                    inputfile=indv_mapping.sorted_bam,
                    default_location='SAMPLE/' + sub(basename(eachfastq),'\.fastq\.gz|\.fq\.gz','') + '/QC/FastQC'
            }

            call runspp.runspp as indv_runspp {
                input:
                    bamfile=select_first([indv_mapping.bklist_bam, indv_mapping.sorted_bam])
            }

            call bedtools.bamtobed as indv_bamtobed {
                input:
                    bamfile=select_first([indv_mapping.bklist_bam, indv_mapping.sorted_bam])
            }

            call util.evalstats as indv_summarystats {
                input:
                    fastq_type="Sample FASTQ",
                    bambed=indv_bamtobed.bedfile,
                    sppfile=indv_runspp.spp_out,
                    fastqczip=indv_fastqc.zipfile,
                    bamflag=indv_mapping.bam_stats,
                    rmdupflag=indv_mapping.mkdup_stats,
                    bkflag=indv_mapping.bklist_stats,
                    fastqmetrics=indv_bfs.metrics_out,
                    default_location='SAMPLE/' + sub(basename(eachfastq),'\.fastq\.gz|\.fq\.gz','') + '/QC/SummaryStats'
            }
        } # end scatter (for each sample fastq)

        # MERGE BAM FILES
        # Execute analysis on merge bam file
        # Analysis executed:
        #   Merge BAM (if more than 1 fastq is provided)
        #   FastQC on Merge BAM (AllMerge_<number>_mapped)

        # merge bam files and perform fasTQC if more than one is provided
        call util.mergehtml {
            input:
                htmlfiles=indv_summarystats.xhtml,
                txtfiles=indv_summarystats.textfile,
                default_location='SAMPLE',
                outputfile = 'AllMapped_' + length(fastqfiles) + '_seaseq-summary-stats.html'
        }

        call samtools.mergebam {
            input:
                bamfiles=indv_mapping.sorted_bam,
                default_location = if defined(results_name) then results_name + '/BAM_files' else 'AllMerge_' + length(indv_mapping.sorted_bam) + '_mapped' + '/BAM_files',
                outputfile = if defined(results_name) then results_name + '.sorted.bam' else 'AllMerge_' + length(fastqfiles) + '_mapped.sorted.bam'
        }

        call fastqc.fastqc as mergebamfqc {
            input:
	        inputfile=mergebam.mergebam,
                default_location=sub(basename(mergebam.mergebam),'\.sorted\.b.*$','') + '/QC/FastQC'
        }

        call samtools.indexstats as mergeindexstats {
            input:
                bamfile=mergebam.mergebam,
                default_location=sub(basename(mergebam.mergebam),'\.sorted\.b.*$','') + '/BAM_files'
        }

        if ( defined(blacklist) ) {
            # remove blacklist regions
            String string_blacklist = "" #buffer to allow for blacklist optionality
            File blacklist_ = select_first([blacklist, string_blacklist])
            call bedtools.intersect as merge_rmblklist {
                input :
                    fileA=mergebam.mergebam,
                    fileB=blacklist_,
                    default_location=sub(basename(mergebam.mergebam),'\.sorted\.b.*$','') + '/BAM_files',
                    nooverlap=true
            }
            call samtools.indexstats as merge_bklist {
                input :
                    bamfile=merge_rmblklist.intersect_out,
                    default_location=sub(basename(mergebam.mergebam),'\.sorted\.b.*$','') + '/BAM_files'
            }
        } # end if blacklist provided

        File mergebam_afterbklist = select_first([merge_rmblklist.intersect_out, mergebam.mergebam])

        call samtools.markdup as merge_markdup {
            input :
                bamfile=mergebam_afterbklist,
                default_location=sub(basename(mergebam_afterbklist),'\.sorted\.b.*$','') + '/BAM_files'
        }

        call samtools.indexstats as merge_mkdup {
            input :
                bamfile=merge_markdup.mkdupbam,
                default_location=sub(basename(mergebam_afterbklist),'\.sorted\.b.*$','') + '/BAM_files'
        }
    } # end if length(fastqfiles) > 1: multi_fastq

### ---------------------------------------- ###
### ------------ S E C T I O N 2 ----------- ###
### -- B: analysis if one FASTQ provided --- ###
### ---------------------------------------- ###

    # if only one fastqfile is provided
    if ( one_fastq ) {
        # Execute analysis on each fastq file provided
        # Analysis executed:
        #   FastQC
        #   FASTQ read length distribution
        #   Reference Alignment using Bowtie (-k2 -m2)
        #   Convert SAM to BAM
        #   FastQC on BAM files
        #   Remove Blacklists (if provided)
        #   Remove read duplicates
        #   Summary statistics on FASTQs
        #   Combine html files into one for easy viewing

        call fastqc.fastqc as uno_fastqc {
            input :
                inputfile=fastqfiles[0],
                default_location=sub(basename(fastqfiles[0]),'\.fastq\.gz|\.fq\.gz','') + '/QC/FastQC'
        }

        call util.basicfastqstats as uno_bfs {
            input :
                fastqfile=fastqfiles[0],
                default_location=sub(basename(fastqfiles[0]),'\.fastq\.gz|\.fq\.gz','') + '/QC/SummaryStats'
        }

        call mapping.mapping {
            input :
                fastqfile=fastqfiles[0],
                index_files=bowtie_index_,
                metricsfile=uno_bfs.metrics_out,
                blacklist=blacklist,
                default_location=sub(basename(fastqfiles[0]),'\.fastq\.gz|\.fq\.gz','') + '/BAM_files'
        }

        call fastqc.fastqc as uno_bamfqc {
            input :
                inputfile=mapping.sorted_bam,
                default_location=sub(basename(fastqfiles[0]),'\.fastq\.gz|\.fq\.gz','') + '/QC/FastQC'
        }

        call runspp.runspp as uno_runspp {
            input:
                bamfile=select_first([mapping.bklist_bam, mapping.sorted_bam])
        }

        call bedtools.bamtobed as uno_bamtobed {
            input:
                bamfile=select_first([mapping.bklist_bam, mapping.sorted_bam])
        }
    } # end if length(fastqfiles) == 1: one_fastq

### ---------------------------------------- ###
### ------------ S E C T I O N 3 ----------- ###
### ----------- ChIP-seq analysis ---------- ###
### ---------------------------------------- ###

    # ChIP-seq and downstream analysis
    # Execute analysis on merge bam file
    # Analysis executed:
    #   FIRST: Check if reads are mapped
    #   Peaks identification (SICER, MACS, ROSE)
    #   Motif analysis
    #   Complete Summary statistics

    #collate correct files for downstream analysis
    File sample_bam = select_first([mergebam_afterbklist, mapping.bklist_bam, mapping.sorted_bam])

    call macs.macs {
        input :
            bamfile=sample_bam,
            pvalue="1e-9",
            keep_dup="auto",
            egs=egs.genomesize,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/PEAKS/NARROW_peaks' + '/' + basename(sample_bam,'\.bam') + '-p9_kd-auto'
    }

    call util.addreadme {
        input :
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/PEAKS'
    }

    call macs.macs as all {
        input :
            bamfile=sample_bam,
            pvalue="1e-9",
            keep_dup="all",
            egs=egs.genomesize,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/PEAKS/NARROW_peaks' + '/' + basename(sample_bam,'\.bam') + '-p9_kd-all'
    }

    call macs.macs as nomodel {
        input :
            bamfile=sample_bam,
            nomodel=true,
            egs=egs.genomesize,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/PEAKS/NARROW_peaks' + '/' + basename(sample_bam,'\.bam') + '-nm'
    }

    call bamtogff.bamtogff {
        input :
            gtffile=gtf,
            chromsizes=samtools_faidx.chromsizes,
            bamfile=select_first([merge_markdup.mkdupbam, mapping.mkdup_bam]),
            bamindex=select_first([merge_mkdup.indexbam, mapping.mkdup_index]),
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/BAM_Density'
    }

    call bedtools.bamtobed as forsicerbed {
        input :
            bamfile=select_first([merge_markdup.mkdupbam, mapping.mkdup_bam])
    }
    
    call sicer.sicer {
        input :
            bedfile=forsicerbed.bedfile,
            chromsizes=samtools_faidx.chromsizes,
            genome_fraction=egs.genomefraction,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/PEAKS/BROAD_peaks'
    }

    call rose.rose {
        input :
            gtffile=gtf,
            bamfile=sample_bam,
            bamindex=select_first([merge_bklist.indexbam, mergeindexstats.indexbam, mapping.bklist_index, mapping.bam_index]),
            bedfile_auto=macs.peakbedfile,
            bedfile_all=all.peakbedfile,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/PEAKS/STITCHED_peaks'
    }

    call runspp.runspp {
        input:
            bamfile=sample_bam
    }

    call util.peaksanno {
        input :
            gtffile=gtf,
            bedfile=macs.peakbedfile,
            chromsizes=samtools_faidx.chromsizes,
            summitfile=macs.summitsfile,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(macs.peakbedfile),'\_peaks.bed','')
    }

    call util.peaksanno as all_peaksanno {
        input :
            gtffile=gtf,
            bedfile=all.peakbedfile,
            chromsizes=samtools_faidx.chromsizes,
            summitfile=all.summitsfile,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(all.peakbedfile),'\_peaks.bed','')
    }

    call util.peaksanno as nomodel_peaksanno {
        input :
            gtffile=gtf,
            bedfile=nomodel.peakbedfile,
            chromsizes=samtools_faidx.chromsizes,
            summitfile=nomodel.summitsfile,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(nomodel.peakbedfile),'\_peaks.bed','')
    }

    call util.peaksanno as sicer_peaksanno {
        input :
            gtffile=gtf,
            bedfile=sicer.scoreisland,
            chromsizes=samtools_faidx.chromsizes,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/PEAKS_Annotation/BROAD_peaks'
    }

    # Motif Analysis
    if (run_motifs) { 
        call motifs.motifs {
            input:
                reference=reference,
                reference_index=samtools_faidx.faidx_file,
                bedfile=macs.peakbedfile,
                motif_databases=motif_databases,
                default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/MOTIFS'
        }

        call util.flankbed {
            input :
                bedfile=macs.summitsfile,
                default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/MOTIFS'
        }

        call motifs.motifs as flank {
            input:
                reference=reference,
                reference_index=samtools_faidx.faidx_file,
                bedfile=flankbed.flankbedfile,
                motif_databases=motif_databases,
                default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/MOTIFS'
        }
    }

    call viz.visualization {
        input:
            wigfile=macs.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            xlsfile=macs.peakxlsfile,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(macs.peakbedfile),'\_peaks.bed','')
    }

    call viz.visualization as vizall {
        input:
            wigfile=all.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            xlsfile=all.peakxlsfile,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(all.peakbedfile),'\_peaks.bed','')
    }

    call viz.visualization as viznomodel {
        input:
            wigfile=nomodel.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            xlsfile=nomodel.peakxlsfile,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(nomodel.peakbedfile),'\_peaks.bed','')
    }

    call viz.visualization as vizsicer {
        input:
            wigfile=sicer.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/COVERAGE_files/BROAD_peaks'
    }

    call bedtools.bamtobed as finalbed {
        input:
            bamfile=sample_bam
    }

    call sortbed.sortbed {
        input:
            bedfile=finalbed.bedfile
    }

    call bedtools.intersect {
        input:
            fileA=macs.peakbedfile,
            fileB=sortbed.sortbed_out,
            countoverlap=true,
            sorted=true
    }

### ---------------------------------------- ###
### ------------ S E C T I O N 4 ----------- ###
### ---------- Summary Statistics ---------- ###
### ---------------------------------------- ###

    String string_qual = "" #buffer to allow for optionality in if statement

    #SUMMARY STATISTICS
    if ( one_fastq ) {
        call util.evalstats as uno_summarystats {
            # SUMMARY STATISTICS of sample file (only 1 sample file provided)
            input:
                fastq_type="Sample FASTQ",
                bambed=finalbed.bedfile,
                sppfile=runspp.spp_out,
                fastqczip=select_first([uno_bamfqc.zipfile, string_qual]),
                bamflag=mapping.bam_stats,
                rmdupflag=mapping.mkdup_stats,
                bkflag=mapping.bklist_stats,
                fastqmetrics=uno_bfs.metrics_out,
                countsfile=intersect.intersect_out,
                peaksxls=macs.peakxlsfile,
                enhancers=rose.enhancers,
                superenhancers=rose.super_enhancers,
                default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/QC/SummaryStats'
        }

        call util.summaryreport as uno_overallsummary {
            # Presenting all quality stats for the analysis
            input:
                overallqc_html=uno_summarystats.xhtml,
                overallqc_txt=uno_summarystats.textfile
        }
    } # end if one_fastq

    if ( multi_fastq ) {
        call util.evalstats as merge_summarystats {
            # SUMMARY STATISTICS of all samples files (more than 1 sample file provided)
            input:
                fastq_type="Comprehensive",
                bambed=finalbed.bedfile,
                sppfile=runspp.spp_out,
                fastqczip=select_first([mergebamfqc.zipfile, string_qual]),
                bamflag=mergeindexstats.flagstats,
                rmdupflag=merge_mkdup.flagstats,
                bkflag=merge_bklist.flagstats,
                countsfile=intersect.intersect_out,
                peaksxls=macs.peakxlsfile,
                enhancers=rose.enhancers,
                superenhancers=rose.super_enhancers,
                default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/QC/SummaryStats'
        }
        
        call util.summaryreport as merge_overallsummary {
            # Presenting all quality stats for the analysis
            input:
                sampleqc_html=mergehtml.xhtml,
                overallqc_html=merge_summarystats.xhtml,
                sampleqc_txt=mergehtml.mergetxt,
                overallqc_txt=merge_summarystats.textfile
        }
    } # end if multi_fastq

    output {
        #FASTQC
        Array[File?]? indv_s_htmlfile = indv_fastqc.htmlfile
        Array[File?]? indv_s_zipfile = indv_fastqc.zipfile
        Array[File?]? indv_s_bam_htmlfile = indv_bamfqc.htmlfile
        Array[File?]? indv_s_bam_zipfile = indv_bamfqc.zipfile

        File? s_mergebam_htmlfile = mergebamfqc.htmlfile
        File? s_mergebam_zipfile = mergebamfqc.zipfile

        File? uno_s_htmlfile = uno_fastqc.htmlfile
        File? uno_s_zipfile = uno_fastqc.zipfile
        File? uno_s_bam_htmlfile = uno_bamfqc.htmlfile
        File? uno_s_bam_zipfile = uno_bamfqc.zipfile

        #BASICMETRICS
        Array[File?]? s_metrics_out = indv_bfs.metrics_out
        File? uno_s_metrics_out = uno_bfs.metrics_out

        #BAMFILES
        Array[File?]? indv_s_sortedbam = indv_mapping.sorted_bam
        Array[File?]? indv_s_indexbam = indv_mapping.bam_index
        Array[File?]? indv_s_bkbam = indv_mapping.bklist_bam
        Array[File?]? indv_s_bkindexbam = indv_mapping.bklist_index
        Array[File?]? indv_s_rmbam = indv_mapping.mkdup_bam
        Array[File?]? indv_s_rmindexbam = indv_mapping.mkdup_index

        File? uno_s_sortedbam = mapping.sorted_bam
        File? uno_s_indexstatsbam = mapping.bam_index
        File? uno_s_bkbam = mapping.bklist_bam
        File? uno_s_bkindexbam = mapping.bklist_index
        File? uno_s_rmbam = mapping.mkdup_bam
        File? uno_s_rmindexbam = mapping.mkdup_index

        File? s_mergebamfile = mergebam.mergebam
        File? s_mergebamindex = mergeindexstats.indexbam
        File? s_bkbam = merge_rmblklist.intersect_out
        File? s_bkindexbam = merge_bklist.indexbam
        File? s_rmbam = merge_markdup.mkdupbam
        File? s_rmindexbam = merge_mkdup.indexbam

        #MACS
        File? peakbedfile = macs.peakbedfile
        File? peakxlsfile = macs.peakxlsfile
        File? summitsfile = macs.summitsfile
        File? negativexlsfile = macs.negativepeaks
        File? wigfile = macs.wigfile
        File? all_peakbedfile = all.peakbedfile
        File? all_peakxlsfile = all.peakxlsfile
        File? all_summitsfile = all.summitsfile
        File? all_negativexlsfile = all.negativepeaks
        File? all_wigfile = all.wigfile
        File? nm_peakbedfile = nomodel.peakbedfile
        File? nm_peakxlsfile = nomodel.peakxlsfile
        File? nm_summitsfile = nomodel.summitsfile
        File? nm_negativexlsfile = nomodel.negativepeaks
        File? nm_wigfile = nomodel.wigfile
        File? readme_peaks = addreadme.readme_peaks

        #SICER
        File? scoreisland = sicer.scoreisland
        File? sicer_wigfile = sicer.wigfile

        #ROSE
        File? pngfile = rose.pngfile
        File? mapped_union = rose.mapped_union
        File? mapped_stitch = rose.mapped_stitch
        File? enhancers = rose.enhancers
        File? super_enhancers = rose.super_enhancers
        File? gff_file = rose.gff_file
        File? gff_union = rose.gff_union
        File? union_enhancers = rose.union_enhancers
        File? stitch_enhancers = rose.stitch_enhancers
        File? e_to_g_enhancers = rose.e_to_g_enhancers
        File? g_to_e_enhancers = rose.g_to_e_enhancers
        File? e_to_g_super_enhancers = rose.e_to_g_super_enhancers
        File? g_to_e_super_enhancers = rose.g_to_e_super_enhancers

        #MOTIFS
        File? flankbedfile = flankbed.flankbedfile

        File? ame_tsv = motifs.ame_tsv
        File? ame_html = motifs.ame_html
        File? ame_seq = motifs.ame_seq
        File? meme = motifs.meme_out
        File? meme_summary = motifs.meme_summary

        File? summit_ame_tsv = flank.ame_tsv
        File? summit_ame_html = flank.ame_html
        File? summit_ame_seq = flank.ame_seq
        File? summit_meme = flank.meme_out
        File? summit_meme_summary = flank.meme_summary

        #BAM2GFF
        File? s_matrices = bamtogff.s_matrices
        File? densityplot = bamtogff.densityplot
        File? pdf_gene = bamtogff.pdf_gene
        File? pdf_h_gene = bamtogff.pdf_h_gene
        File? png_h_gene = bamtogff.png_h_gene
        File? pdf_promoters = bamtogff.pdf_promoters
        File? pdf_h_promoters = bamtogff.pdf_h_promoters
        File? png_h_promoters = bamtogff.png_h_promoters

        #PEAKS-ANNOTATION
        File? peak_promoters = peaksanno.peak_promoters
        File? peak_genebody = peaksanno.peak_genebody
        File? peak_window = peaksanno.peak_window
        File? peak_closest = peaksanno.peak_closest
        File? peak_comparison = peaksanno.peak_comparison
        File? gene_comparison = peaksanno.gene_comparison
        File? pdf_comparison = peaksanno.pdf_comparison

        File? all_peak_promoters = all_peaksanno.peak_promoters
        File? all_peak_genebody = all_peaksanno.peak_genebody
        File? all_peak_window = all_peaksanno.peak_window
        File? all_peak_closest = all_peaksanno.peak_closest
        File? all_peak_comparison = all_peaksanno.peak_comparison
        File? all_gene_comparison = all_peaksanno.gene_comparison
        File? all_pdf_comparison = all_peaksanno.pdf_comparison

        File? nomodel_peak_promoters = nomodel_peaksanno.peak_promoters
        File? nomodel_peak_genebody = nomodel_peaksanno.peak_genebody
        File? nomodel_peak_window = nomodel_peaksanno.peak_window
        File? nomodel_peak_closest = nomodel_peaksanno.peak_closest
        File? nomodel_peak_comparison = nomodel_peaksanno.peak_comparison
        File? nomodel_gene_comparison = nomodel_peaksanno.gene_comparison
        File? nomodel_pdf_comparison = nomodel_peaksanno.pdf_comparison

        File? sicer_peak_promoters = sicer_peaksanno.peak_promoters
        File? sicer_peak_genebody = sicer_peaksanno.peak_genebody
        File? sicer_peak_window = sicer_peaksanno.peak_window
        File? sicer_peak_closest = sicer_peaksanno.peak_closest
        File? sicer_peak_comparison = sicer_peaksanno.peak_comparison
        File? sicer_gene_comparison = sicer_peaksanno.gene_comparison
        File? sicer_pdf_comparison = sicer_peaksanno.pdf_comparison

        #VISUALIZATION
        File? bigwig = visualization.bigwig
        File? norm_wig = visualization.norm_wig
        File? tdffile = visualization.tdffile
        File? n_bigwig = viznomodel.bigwig
        File? n_norm_wig = viznomodel.norm_wig
        File? n_tdffile = viznomodel.tdffile
        File? a_bigwig = vizall.bigwig
        File? a_norm_wig = vizall.norm_wig
        File? a_tdffile = vizall.tdffile

        File? s_bigwig = vizsicer.bigwig
        File? s_norm_wig = vizsicer.norm_wig
        File? s_tdffile = vizsicer.tdffile

        #QC-STATS
        Array[File?]? s_qc_statsfile = indv_summarystats.statsfile
        Array[File?]? s_qc_htmlfile = indv_summarystats.htmlfile
        Array[File?]? s_qc_textfile = indv_summarystats.textfile
        File? s_qc_mergehtml = mergehtml.mergefile

        File? s_uno_statsfile = uno_summarystats.statsfile
        File? s_uno_htmlfile = uno_summarystats.htmlfile
        File? s_uno_textfile = uno_summarystats.textfile

        File? statsfile = merge_summarystats.statsfile
        File? htmlfile = merge_summarystats.htmlfile
        File? textfile = merge_summarystats.textfile

        File? summaryhtml = uno_overallsummary.summaryhtml
        File? summarytxt = uno_overallsummary.summarytxt
        File? m_summaryhtml = merge_overallsummary.summaryhtml
        File? m_summarytxt = merge_overallsummary.summarytxt
    }
}
