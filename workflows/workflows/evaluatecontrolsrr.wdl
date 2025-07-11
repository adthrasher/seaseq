version 1.0

import "../../peaseq-control.wdl" as pc
import "../../seaseq-control.wdl" as sc
import "../../workflows/tasks/sratoolkit.wdl" as sra

workflow evaluatesrr {
    input {
        File reference
        File? spikein_reference
        File? blacklist
        File gtf
        Array[File]? bowtie_index
        Array[File]? spikein_bowtie_index
        Array[File]? motif_databases
        Array[String]? sample_sraid
        Array[File]? sample_R1_fastq
        Array[File]? sample_R2_fastq
        Array[String]? control_sraid
        Array[File]? control_R1_fastq
        Array[File]? control_R2_fastq
        String? results_name
        Boolean run_motifs = true
        Int insertsize = 600
        String strandedness = "fr"
        Boolean paired = true
    }

    if (defined(sample_sraid)) {
        Array[String] string_sra = [
            "1",
        ]  # buffer to allow for sra_id optionality
        Array[String] s_sraid = select_first([
            sample_sraid,
            string_sra,
        ])
        scatter (eachsra in s_sraid) {
            call sra.srameta { input: sra_id = eachsra }
            if (!srameta.paired_end) {
                Boolean? paired_sample = srameta.paired_end
            }
        }  # end scatter each sra
    }
    Boolean paired_sample_m = select_first([
        paired_sample[0],
        paired,
    ])

    if (defined(control_sraid)) {
        Array[String] c_sra = [
            "1",
        ]  # buffer to allow for sra_id optionality
        Array[String] c_sraid = select_first([
            control_sraid,
            c_sra,
        ])
        scatter (eachsra in c_sraid) {
            call sra.srameta as c_srameta { input: sra_id = eachsra }
            if (!c_srameta.paired_end) {
                Boolean? paired_control = c_srameta.paired_end
            }
        }  # end scatter each sra
    }  # end if control_sra_id
    Boolean paired_control_m = select_first([
        paired_control[0],
        paired,
    ])

    if (!(paired_control_m && paired_sample_m)) {
        call sc.seaseq as sc_seaseq { input:
            reference = reference,
            spikein_reference = spikein_reference,
            blacklist = blacklist,
            gtf = gtf,
            bowtie_index = bowtie_index,
            spikein_bowtie_index = spikein_bowtie_index,
            motif_databases = motif_databases,
            sample_fastq = sample_R1_fastq,
            control_fastq = control_R1_fastq,
            sample_sraid = sample_sraid,
            control_sraid = control_sraid,
            results_name = results_name,
            run_motifs = run_motifs,
        }
    }
    if (paired_sample_m && paired_control_m) {
        call pc.peaseq as pc_peaseq { input:
            reference = reference,
            spikein_reference = spikein_reference,
            blacklist = blacklist,
            gtf = gtf,
            bowtie_index = bowtie_index,
            spikein_bowtie_index = spikein_bowtie_index,
            motif_databases = motif_databases,
            sample_R1_fastq = sample_R1_fastq,
            sample_R2_fastq = sample_R2_fastq,
            control_R1_fastq = control_R1_fastq,
            control_R2_fastq = control_R2_fastq,
            sample_sraid = sample_sraid,
            control_sraid = control_sraid,
            insertsize = insertsize,
            strandedness = strandedness,
            results_name = results_name,
            run_motifs = run_motifs,
        }
    }

    # Processing OUTPUTs
    output {
        Array[File?]? spikein_indv_s_htmlfile = if (paired_sample_m && paired_control_m)
            then pc_peaseq.spikein_indv_s_htmlfile else sc_seaseq.spikein_indv_s_htmlfile
        Array[File?]? spikein_indv_s_zipfile = if (paired_sample_m && paired_control_m)
            then pc_peaseq.spikein_indv_s_zipfile else sc_seaseq.spikein_indv_s_zipfile
        Array[File?]? spikein_s_metrics_out = if (paired_sample_m && paired_control_m)
            then pc_peaseq.spikein_s_metrics_out else sc_seaseq.spikein_s_metrics_out
        Array[File?]? spikein_indv_c_htmlfile = if (paired_sample_m && paired_control_m)
            then pc_peaseq.spikein_indv_c_htmlfile else sc_seaseq.spikein_indv_c_htmlfile
        Array[File?]? spikein_indv_c_zipfile = if (paired_sample_m && paired_control_m)
            then pc_peaseq.spikein_indv_c_zipfile else sc_seaseq.spikein_indv_c_zipfile
        Array[File?]? spikein_c_metrics_out = if (paired_sample_m && paired_control_m)
            then pc_peaseq.spikein_c_metrics_out else sc_seaseq.spikein_c_metrics_out
        Array[File?]? indv_s_htmlfile = if (paired_sample_m && paired_control_m) then pc_peaseq.indv_s_htmlfile
            else sc_seaseq.indv_s_htmlfile
        Array[File?]? indv_s_zipfile = if (paired_sample_m && paired_control_m) then pc_peaseq.indv_s_zipfile
            else sc_seaseq.indv_s_zipfile
        Array[File?]? indv_s_bam_htmlfile = if (paired_sample_m && paired_control_m) then pc_peaseq.indv_s_bam_htmlfile
            else sc_seaseq.indv_s_bam_htmlfile
        Array[File?]? indv_s_bam_zipfile = if (paired_sample_m && paired_control_m) then pc_peaseq.indv_s_bam_zipfile
            else sc_seaseq.indv_s_bam_zipfile
        Array[File?]? indv_c_htmlfile = if (paired_sample_m && paired_control_m) then pc_peaseq.indv_c_htmlfile
            else sc_seaseq.indv_c_htmlfile
        Array[File?]? indv_c_zipfile = if (paired_sample_m && paired_control_m) then pc_peaseq.indv_c_zipfile
            else sc_seaseq.indv_c_zipfile
        Array[File?]? indv_c_bam_htmlfile = if (paired_sample_m && paired_control_m) then pc_peaseq.indv_c_bam_htmlfile
            else sc_seaseq.indv_c_bam_htmlfile
        Array[File?]? indv_c_bam_zipfile = if (paired_sample_m && paired_control_m) then pc_peaseq.indv_c_bam_zipfile
            else sc_seaseq.indv_c_bam_zipfile
        File? s_mergebam_htmlfile = if (paired_sample_m && paired_control_m) then pc_peaseq.s_mergebam_htmlfile
            else sc_seaseq.s_mergebam_htmlfile
        File? s_mergebam_zipfile = if (paired_sample_m && paired_control_m) then pc_peaseq.s_mergebam_zipfile
            else sc_seaseq.s_mergebam_zipfile
        File? c_mergebam_htmlfile = if (paired_sample_m && paired_control_m) then pc_peaseq.c_mergebam_htmlfile
            else sc_seaseq.c_mergebam_htmlfile
        File? c_mergebam_zipfile = if (paired_sample_m && paired_control_m) then pc_peaseq.c_mergebam_zipfile
            else sc_seaseq.c_mergebam_zipfile
        Array[File?]? indv_sp_bam_htmlfile = pc_peaseq.indv_sp_bam_htmlfile
        Array[File?]? indv_sp_bam_zipfile = pc_peaseq.indv_sp_bam_zipfile
        File? sp_mergebam_htmlfile = pc_peaseq.sp_mergebam_htmlfile
        File? sp_mergebam_zipfile = pc_peaseq.sp_mergebam_zipfile
        Array[File?]? indv_cp_bam_htmlfile = pc_peaseq.indv_cp_bam_htmlfile
        Array[File?]? indv_cp_bam_zipfile = pc_peaseq.indv_cp_bam_zipfile
        File? cp_mergebam_htmlfile = pc_peaseq.cp_mergebam_htmlfile
        File? cp_mergebam_zipfile = pc_peaseq.cp_mergebam_zipfile
        File? uno_s_htmlfile = sc_seaseq.uno_s_htmlfile
        File? uno_s_zipfile = sc_seaseq.uno_s_zipfile
        File? uno_s_bam_htmlfile = if (paired_sample_m && paired_control_m) then pc_peaseq.uno_s_bam_htmlfile
            else sc_seaseq.uno_s_bam_htmlfile
        File? uno_s_bam_zipfile = if (paired_sample_m && paired_control_m) then pc_peaseq.uno_s_bam_zipfile
            else sc_seaseq.uno_s_bam_zipfile
        File? uno_c_htmlfile = sc_seaseq.uno_c_htmlfile
        File? uno_c_zipfile = sc_seaseq.uno_c_zipfile
        File? uno_c_bam_htmlfile = if (paired_sample_m && paired_control_m) then pc_peaseq.uno_c_bam_htmlfile
            else sc_seaseq.uno_c_bam_htmlfile
        File? uno_c_bam_zipfile = if (paired_sample_m && paired_control_m) then pc_peaseq.uno_c_bam_zipfile
            else sc_seaseq.uno_c_bam_zipfile
        Array[File?]? s_metrics_out = if (paired_sample_m && paired_control_m) then pc_peaseq.s_metrics_out
            else sc_seaseq.s_metrics_out
        File? uno_s_metrics_out = sc_seaseq.uno_s_metrics_out
        Array[File?]? c_metrics_out = if (paired_sample_m && paired_control_m) then pc_peaseq.c_metrics_out
            else sc_seaseq.c_metrics_out
        File? uno_c_metrics_out = sc_seaseq.uno_c_metrics_out
        Array[File?]? indv_s_sortedbam = if (paired_sample_m && paired_control_m) then pc_peaseq.indv_s_sortedbam
            else sc_seaseq.indv_s_sortedbam
        Array[File?]? indv_s_indexbam = if (paired_sample_m && paired_control_m) then pc_peaseq.indv_s_indexbam
            else sc_seaseq.indv_s_indexbam
        Array[File?]? indv_s_bkbam = if (paired_sample_m && paired_control_m) then pc_peaseq.indv_s_bkbam
            else sc_seaseq.indv_s_bkbam
        Array[File?]? indv_s_bkindexbam = if (paired_sample_m && paired_control_m) then pc_peaseq.indv_s_bkindexbam
            else sc_seaseq.indv_s_bkindexbam
        Array[File?]? indv_s_rmbam = if (paired_sample_m && paired_control_m) then pc_peaseq.indv_s_rmbam
            else sc_seaseq.indv_s_rmbam
        Array[File?]? indv_s_rmindexbam = if (paired_sample_m && paired_control_m) then pc_peaseq.indv_s_rmindexbam
            else sc_seaseq.indv_s_rmindexbam
        Array[File?]? indv_c_sortedbam = if (paired_sample_m && paired_control_m) then pc_peaseq.indv_c_sortedbam
            else sc_seaseq.indv_c_sortedbam
        Array[File?]? indv_c_indexbam = if (paired_sample_m && paired_control_m) then pc_peaseq.indv_c_indexbam
            else sc_seaseq.indv_c_indexbam
        Array[File?]? indv_c_bkbam = if (paired_sample_m && paired_control_m) then pc_peaseq.indv_c_bkbam
            else sc_seaseq.indv_c_bkbam
        Array[File?]? indv_c_bkindexbam = if (paired_sample_m && paired_control_m) then pc_peaseq.indv_c_bkindexbam
            else sc_seaseq.indv_c_bkindexbam
        Array[File?]? indv_c_rmbam = if (paired_sample_m && paired_control_m) then pc_peaseq.indv_c_rmbam
            else sc_seaseq.indv_c_rmbam
        Array[File?]? indv_c_rmindexbam = if (paired_sample_m && paired_control_m) then pc_peaseq.indv_c_rmindexbam
            else sc_seaseq.indv_c_rmindexbam
        Array[File?]? indv_sp_sortedbam = pc_peaseq.indv_sp_sortedbam
        Array[File?]? indv_sp_indexbam = pc_peaseq.indv_sp_indexbam
        Array[File?]? indv_sp_bkbam = pc_peaseq.indv_sp_bkbam
        Array[File?]? indv_sp_bkindexbam = pc_peaseq.indv_sp_bkindexbam
        Array[File?]? indv_sp_rmbam = pc_peaseq.indv_sp_rmbam
        Array[File?]? indv_sp_rmindexbam = pc_peaseq.indv_sp_rmindexbam
        Array[File?]? indv_cp_sortedbam = pc_peaseq.indv_cp_sortedbam
        Array[File?]? indv_cp_indexbam = pc_peaseq.indv_cp_indexbam
        Array[File?]? indv_cp_bkbam = pc_peaseq.indv_cp_bkbam
        Array[File?]? indv_cp_bkindexbam = pc_peaseq.indv_cp_bkindexbam
        Array[File?]? indv_cp_rmbam = pc_peaseq.indv_cp_rmbam
        Array[File?]? indv_cp_rmindexbam = pc_peaseq.indv_cp_rmindexbam
        File? uno_s_sortedbam = if (paired_sample_m && paired_control_m) then pc_peaseq.uno_s_sortedbam
            else sc_seaseq.uno_s_sortedbam
        File? uno_s_indexstatsbam = if (paired_sample_m && paired_control_m) then pc_peaseq.uno_s_indexstatsbam
            else sc_seaseq.uno_s_indexstatsbam
        File? uno_s_bkbam = if (paired_sample_m && paired_control_m) then pc_peaseq.uno_s_bkbam
            else sc_seaseq.uno_s_bkbam
        File? uno_s_bkindexbam = if (paired_sample_m && paired_control_m) then pc_peaseq.uno_s_bkindexbam
            else sc_seaseq.uno_s_bkindexbam
        File? uno_s_rmbam = if (paired_sample_m && paired_control_m) then pc_peaseq.uno_s_rmbam
            else sc_seaseq.uno_s_rmbam
        File? uno_s_rmindexbam = if (paired_sample_m && paired_control_m) then pc_peaseq.uno_s_rmindexbam
            else sc_seaseq.uno_s_rmindexbam
        File? uno_c_sortedbam = if (paired_sample_m && paired_control_m) then pc_peaseq.uno_c_sortedbam
            else sc_seaseq.uno_c_sortedbam
        File? uno_c_indexstatsbam = if (paired_sample_m && paired_control_m) then pc_peaseq.uno_c_indexstatsbam
            else sc_seaseq.uno_c_indexstatsbam
        File? uno_c_bkbam = if (paired_sample_m && paired_control_m) then pc_peaseq.uno_c_bkbam
            else sc_seaseq.uno_c_bkbam
        File? uno_c_bkindexbam = if (paired_sample_m && paired_control_m) then pc_peaseq.uno_c_bkindexbam
            else sc_seaseq.uno_c_bkindexbam
        File? uno_c_rmbam = if (paired_sample_m && paired_control_m) then pc_peaseq.uno_c_rmbam
            else sc_seaseq.uno_c_rmbam
        File? uno_c_rmindexbam = if (paired_sample_m && paired_control_m) then pc_peaseq.uno_c_rmindexbam
            else sc_seaseq.uno_c_rmindexbam
        File? s_mergebamfile = if (paired_sample_m && paired_control_m) then pc_peaseq.s_mergebamfile
            else sc_seaseq.s_mergebamfile
        File? s_mergebamindex = if (paired_sample_m && paired_control_m) then pc_peaseq.s_mergebamindex
            else sc_seaseq.s_mergebamindex
        File? s_bkbam = if (paired_sample_m && paired_control_m) then pc_peaseq.s_bkbam else sc_seaseq.s_bkbam
        File? s_bkindexbam = if (paired_sample_m && paired_control_m) then pc_peaseq.s_bkindexbam
            else sc_seaseq.s_bkindexbam
        File? s_rmbam = if (paired_sample_m && paired_control_m) then pc_peaseq.s_rmbam else sc_seaseq.s_rmbam
        File? s_rmindexbam = if (paired_sample_m && paired_control_m) then pc_peaseq.s_rmindexbam
            else sc_seaseq.s_rmindexbam
        File? c_mergebamfile = if (paired_sample_m && paired_control_m) then pc_peaseq.c_mergebamfile
            else sc_seaseq.c_mergebamfile
        File? c_mergebamindex = if (paired_sample_m && paired_control_m) then pc_peaseq.c_mergebamindex
            else sc_seaseq.c_mergebamindex
        File? c_bkbam = if (paired_sample_m && paired_control_m) then pc_peaseq.c_bkbam else sc_seaseq.c_bkbam
        File? c_bkindexbam = if (paired_sample_m && paired_control_m) then pc_peaseq.c_bkindexbam
            else sc_seaseq.c_bkindexbam
        File? c_rmbam = if (paired_sample_m && paired_control_m) then pc_peaseq.c_rmbam else sc_seaseq.c_rmbam
        File? c_rmindexbam = if (paired_sample_m && paired_control_m) then pc_peaseq.c_rmindexbam
            else sc_seaseq.c_rmindexbam
        File? sp_mergebamfile = pc_peaseq.sp_mergebamfile
        File? sp_mergebamindex = pc_peaseq.sp_mergebamindex
        File? sp_bkbam = pc_peaseq.sp_bkbam
        File? sp_bkindexbam = pc_peaseq.sp_bkindexbam
        File? sp_rmbam = pc_peaseq.sp_rmbam
        File? sp_rmindexbam = pc_peaseq.sp_rmindexbam
        File? cp_mergebamfile = pc_peaseq.cp_mergebamfile
        File? cp_mergebamindex = pc_peaseq.cp_mergebamindex
        File? cp_bkbam = pc_peaseq.cp_bkbam
        File? cp_bkindexbam = pc_peaseq.cp_bkindexbam
        File? cp_rmbam = pc_peaseq.cp_rmbam
        File? cp_rmindexbam = pc_peaseq.cp_rmindexbam
        File? s_fragments_bam = pc_peaseq.s_fragments_bam
        File? s_fragments_indexbam = pc_peaseq.s_fragments_indexbam
        File? c_fragments_bam = pc_peaseq.c_fragments_bam
        File? c_fragments_indexbam = pc_peaseq.c_fragments_indexbam
        File? peakbedfile = if (paired_sample_m && paired_control_m) then pc_peaseq.peakbedfile
            else sc_seaseq.peakbedfile
        File? peakxlsfile = if (paired_sample_m && paired_control_m) then pc_peaseq.peakxlsfile
            else sc_seaseq.peakxlsfile
        File? negativexlsfile = if (paired_sample_m && paired_control_m) then pc_peaseq.negativexlsfile
            else sc_seaseq.negativexlsfile
        File? summitsfile = if (paired_sample_m && paired_control_m) then pc_peaseq.summitsfile
            else sc_seaseq.summitsfile
        File? wigfile = if (paired_sample_m && paired_control_m) then pc_peaseq.wigfile else sc_seaseq.wigfile
        File? ctrlwigfile = if (paired_sample_m && paired_control_m) then pc_peaseq.ctrlwigfile
            else sc_seaseq.ctrlwigfile
        File? all_peakbedfile = if (paired_sample_m && paired_control_m) then pc_peaseq.all_peakbedfile
            else sc_seaseq.all_peakbedfile
        File? all_peakxlsfile = if (paired_sample_m && paired_control_m) then pc_peaseq.all_peakxlsfile
            else sc_seaseq.all_peakxlsfile
        File? all_negativexlsfile = if (paired_sample_m && paired_control_m) then pc_peaseq.all_negativexlsfile
            else sc_seaseq.all_negativexlsfile
        File? all_summitsfile = if (paired_sample_m && paired_control_m) then pc_peaseq.all_summitsfile
            else sc_seaseq.all_summitsfile
        File? all_wigfile = if (paired_sample_m && paired_control_m) then pc_peaseq.all_wigfile
            else sc_seaseq.all_wigfile
        File? all_ctrlwigfile = if (paired_sample_m && paired_control_m) then pc_peaseq.all_ctrlwigfile
            else sc_seaseq.all_ctrlwigfile
        File? nm_peakbedfile = if (paired_sample_m && paired_control_m) then pc_peaseq.nm_peakbedfile
            else sc_seaseq.nm_peakbedfile
        File? nm_peakxlsfile = if (paired_sample_m && paired_control_m) then pc_peaseq.nm_peakxlsfile
            else sc_seaseq.nm_peakxlsfile
        File? nm_negativexlsfile = if (paired_sample_m && paired_control_m) then pc_peaseq.nm_negativexlsfile
            else sc_seaseq.nm_negativexlsfile
        File? nm_summitsfile = if (paired_sample_m && paired_control_m) then pc_peaseq.nm_summitsfile
            else sc_seaseq.nm_summitsfile
        File? nm_wigfile = if (paired_sample_m && paired_control_m) then pc_peaseq.nm_wigfile
            else sc_seaseq.nm_wigfile
        File? nm_ctrlwigfile = if (paired_sample_m && paired_control_m) then pc_peaseq.nm_ctrlwigfile
            else sc_seaseq.nm_ctrlwigfile
        File? readme_peaks = if (paired_sample_m && paired_control_m) then pc_peaseq.readme_peaks
            else sc_seaseq.readme_peaks
        File? only_s_peakbedfile = if (paired_sample_m && paired_control_m) then pc_peaseq.only_s_peakbedfile
            else sc_seaseq.only_s_peakbedfile
        File? only_s_peakxlsfile = if (paired_sample_m && paired_control_m) then pc_peaseq.only_s_peakxlsfile
            else sc_seaseq.only_s_peakxlsfile
        File? only_s_summitsfile = if (paired_sample_m && paired_control_m) then pc_peaseq.only_s_summitsfile
            else sc_seaseq.only_s_summitsfile
        File? only_s_wigfile = if (paired_sample_m && paired_control_m) then pc_peaseq.only_s_wigfile
            else sc_seaseq.only_s_wigfile
        File? only_c_peakbedfile = if (paired_sample_m && paired_control_m) then pc_peaseq.only_c_peakbedfile
            else sc_seaseq.only_c_peakbedfile
        File? only_c_peakxlsfile = if (paired_sample_m && paired_control_m) then pc_peaseq.only_c_peakxlsfile
            else sc_seaseq.only_c_peakxlsfile
        File? only_c_summitsfile = if (paired_sample_m && paired_control_m) then pc_peaseq.only_c_summitsfile
            else sc_seaseq.only_c_summitsfile
        File? only_c_wigfile = if (paired_sample_m && paired_control_m) then pc_peaseq.only_c_wigfile
            else sc_seaseq.only_c_wigfile
        File? only_sp_peakbedfile = pc_peaseq.only_sp_peakbedfile
        File? only_sp_peakxlsfile = pc_peaseq.only_sp_peakxlsfile
        File? only_sp_summitsfile = pc_peaseq.only_sp_summitsfile
        File? only_sp_wigfile = pc_peaseq.only_sp_wigfile
        File? only_cp_peakbedfile = pc_peaseq.only_cp_peakbedfile
        File? only_cp_peakxlsfile = pc_peaseq.only_cp_peakxlsfile
        File? only_cp_summitsfile = pc_peaseq.only_cp_summitsfile
        File? only_cp_wigfile = pc_peaseq.only_cp_wigfile
        File? sp_peakbedfile = pc_peaseq.sp_peakbedfile
        File? sp_peakxlsfile = pc_peaseq.sp_peakxlsfile
        File? sp_negativexlsfile = pc_peaseq.sp_negativexlsfile
        File? sp_summitsfile = pc_peaseq.sp_summitsfile
        File? sp_wigfile = pc_peaseq.sp_wigfile
        File? sp_ctrlwigfile = pc_peaseq.sp_ctrlwigfile
        File? sp_all_peakbedfile = pc_peaseq.sp_all_peakbedfile
        File? sp_all_peakxlsfile = pc_peaseq.sp_all_peakxlsfile
        File? sp_all_negativexlsfile = pc_peaseq.sp_all_negativexlsfile
        File? sp_all_summitsfile = pc_peaseq.sp_all_summitsfile
        File? sp_all_wigfile = pc_peaseq.sp_all_wigfile
        File? sp_all_ctrlwigfile = pc_peaseq.sp_all_ctrlwigfile
        File? sp_nm_peakbedfile = pc_peaseq.sp_nm_peakbedfile
        File? sp_nm_peakxlsfile = pc_peaseq.sp_nm_peakxlsfile
        File? sp_nm_negativexlsfile = pc_peaseq.sp_nm_negativexlsfile
        File? sp_nm_summitsfile = pc_peaseq.sp_nm_summitsfile
        File? sp_nm_wigfile = pc_peaseq.sp_nm_wigfile
        File? sp_nm_ctrlwigfile = pc_peaseq.sp_nm_ctrlwigfile
        File? sp_readme_peaks = pc_peaseq.sp_readme_peaks
        File? sc_scoreisland = if (paired_sample_m && paired_control_m) then pc_peaseq.scoreisland
            else sc_seaseq.scoreisland
        File? sicer_wigfile = if (paired_sample_m && paired_control_m) then pc_peaseq.sicer_wigfile
            else sc_seaseq.sicer_wigfile
        File? sicer_summary = if (paired_sample_m && paired_control_m) then pc_peaseq.sicer_summary
            else sc_seaseq.sicer_summary
        File? sicer_fdrisland = if (paired_sample_m && paired_control_m) then pc_peaseq.sicer_fdrisland
            else sc_seaseq.sicer_fdrisland
        File? sp_scoreisland = pc_peaseq.sp_scoreisland
        File? sp_sicer_wigfile = pc_peaseq.sp_sicer_wigfile
        File? sp_sicer_summary = pc_peaseq.sp_sicer_summary
        File? sp_sicer_fdrisland = pc_peaseq.sp_sicer_fdrisland
        File? pngfile = if (paired_sample_m && paired_control_m) then pc_peaseq.pngfile else sc_seaseq.pngfile
        File? mapped_union = if (paired_sample_m && paired_control_m) then pc_peaseq.mapped_union
            else sc_seaseq.mapped_union
        File? mapped_stitch = if (paired_sample_m && paired_control_m) then pc_peaseq.mapped_stitch
            else sc_seaseq.mapped_stitch
        File? enhancers = if (paired_sample_m && paired_control_m) then pc_peaseq.enhancers else sc_seaseq.enhancers
        File? super_enhancers = if (paired_sample_m && paired_control_m) then pc_peaseq.super_enhancers
            else sc_seaseq.super_enhancers
        File? gff_file = if (paired_sample_m && paired_control_m) then pc_peaseq.gff_file else sc_seaseq.gff_file
        File? gff_union = if (paired_sample_m && paired_control_m) then pc_peaseq.gff_union else sc_seaseq.gff_union
        File? union_enhancers = if (paired_sample_m && paired_control_m) then pc_peaseq.union_enhancers
            else sc_seaseq.union_enhancers
        File? stitch_enhancers = if (paired_sample_m && paired_control_m) then pc_peaseq.stitch_enhancers
            else sc_seaseq.stitch_enhancers
        File? e_to_g_enhancers = if (paired_sample_m && paired_control_m) then pc_peaseq.e_to_g_enhancers
            else sc_seaseq.e_to_g_enhancers
        File? g_to_e_enhancers = if (paired_sample_m && paired_control_m) then pc_peaseq.g_to_e_enhancers
            else sc_seaseq.g_to_e_enhancers
        File? e_to_g_super_enhancers = if (paired_sample_m && paired_control_m) then pc_peaseq.e_to_g_super_enhancers
            else sc_seaseq.e_to_g_super_enhancers
        File? g_to_e_super_enhancers = if (paired_sample_m && paired_control_m) then pc_peaseq.g_to_e_super_enhancers
            else sc_seaseq.g_to_e_super_enhancers
        File? sp_pngfile = pc_peaseq.sp_pngfile
        File? sp_mapped_union = pc_peaseq.sp_mapped_union
        File? sp_mapped_stitch = pc_peaseq.sp_mapped_stitch
        File? sp_enhancers = pc_peaseq.sp_enhancers
        File? sp_super_enhancers = pc_peaseq.sp_super_enhancers
        File? sp_gff_file = pc_peaseq.sp_gff_file
        File? sp_gff_union = pc_peaseq.sp_gff_union
        File? sp_union_enhancers = pc_peaseq.sp_union_enhancers
        File? sp_stitch_enhancers = pc_peaseq.sp_stitch_enhancers
        File? sp_e_to_g_enhancers = pc_peaseq.sp_e_to_g_enhancers
        File? sp_g_to_e_enhancers = pc_peaseq.sp_g_to_e_enhancers
        File? sp_e_to_g_super_enhancers = pc_peaseq.sp_e_to_g_super_enhancers
        File? sp_g_to_e_super_enhancers = pc_peaseq.sp_g_to_e_super_enhancers
        File? flankbedfile = if (paired_sample_m && paired_control_m) then pc_peaseq.flankbedfile
            else sc_seaseq.flankbedfile
        File? ame_tsv = if (paired_sample_m && paired_control_m) then pc_peaseq.ame_tsv else sc_seaseq.ame_tsv
        File? ame_html = if (paired_sample_m && paired_control_m) then pc_peaseq.ame_html else sc_seaseq.ame_html
        File? ame_seq = if (paired_sample_m && paired_control_m) then pc_peaseq.ame_seq else sc_seaseq.ame_seq
        File? meme = if (paired_sample_m && paired_control_m) then pc_peaseq.meme else sc_seaseq.meme
        File? meme_summary = if (paired_sample_m && paired_control_m) then pc_peaseq.meme_summary
            else sc_seaseq.meme_summary
        File? summit_ame_tsv = if (paired_sample_m && paired_control_m) then pc_peaseq.summit_ame_tsv
            else sc_seaseq.summit_ame_tsv
        File? summit_ame_html = if (paired_sample_m && paired_control_m) then pc_peaseq.summit_ame_html
            else sc_seaseq.summit_ame_html
        File? summit_ame_seq = if (paired_sample_m && paired_control_m) then pc_peaseq.summit_ame_seq
            else sc_seaseq.summit_ame_seq
        File? summit_meme = if (paired_sample_m && paired_control_m) then pc_peaseq.summit_meme
            else sc_seaseq.summit_meme
        File? summit_meme_summary = if (paired_sample_m && paired_control_m) then pc_peaseq.summit_meme_summary
            else sc_seaseq.summit_meme_summary
        File? sp_flankbedfile = pc_peaseq.sp_flankbedfile
        File? sp_ame_tsv = pc_peaseq.sp_ame_tsv
        File? sp_ame_html = pc_peaseq.sp_ame_html
        File? sp_ame_seq = pc_peaseq.sp_ame_seq
        File? sp_meme = pc_peaseq.sp_meme
        File? sp_meme_summary = pc_peaseq.sp_meme_summary
        File? sp_summit_ame_tsv = pc_peaseq.sp_summit_ame_tsv
        File? sp_summit_ame_html = pc_peaseq.sp_summit_ame_html
        File? sp_summit_ame_seq = pc_peaseq.sp_summit_ame_seq
        File? sp_summit_meme = pc_peaseq.sp_summit_meme
        File? sp_summit_meme_summary = pc_peaseq.sp_summit_meme_summary
        File? s_matrices = if (paired_sample_m && paired_control_m) then pc_peaseq.s_matrices
            else sc_seaseq.s_matrices
        File? c_matrices = if (paired_sample_m && paired_control_m) then pc_peaseq.c_matrices
            else sc_seaseq.c_matrices
        File? densityplot = if (paired_sample_m && paired_control_m) then pc_peaseq.densityplot
            else sc_seaseq.densityplot
        File? pdf_gene = if (paired_sample_m && paired_control_m) then pc_peaseq.pdf_gene else sc_seaseq.pdf_gene
        File? pdf_h_gene = if (paired_sample_m && paired_control_m) then pc_peaseq.pdf_h_gene
            else sc_seaseq.pdf_h_gene
        File? png_h_gene = if (paired_sample_m && paired_control_m) then pc_peaseq.png_h_gene
            else sc_seaseq.png_h_gene
        File? jpg_h_gene = if (paired_sample_m && paired_control_m) then pc_peaseq.jpg_h_gene
            else sc_seaseq.jpg_h_gene
        File? pdf_promoters = if (paired_sample_m && paired_control_m) then pc_peaseq.pdf_promoters
            else sc_seaseq.pdf_promoters
        File? pdf_h_promoters = if (paired_sample_m && paired_control_m) then pc_peaseq.pdf_h_promoters
            else sc_seaseq.pdf_h_promoters
        File? png_h_promoters = if (paired_sample_m && paired_control_m) then pc_peaseq.png_h_promoters
            else sc_seaseq.png_h_promoters
        File? jpg_h_promoters = if (paired_sample_m && paired_control_m) then pc_peaseq.jpg_h_promoters
            else sc_seaseq.jpg_h_promoters
        File? sp_s_matrices = pc_peaseq.sp_s_matrices
        File? sp_c_matrices = pc_peaseq.sp_c_matrices
        File? sp_densityplot = pc_peaseq.sp_densityplot
        File? sp_pdf_gene = pc_peaseq.sp_pdf_gene
        File? sp_pdf_h_gene = pc_peaseq.sp_pdf_h_gene
        File? sp_png_h_gene = pc_peaseq.sp_png_h_gene
        File? sp_jpg_h_gene = pc_peaseq.sp_jpg_h_gene
        File? sp_pdf_promoters = pc_peaseq.sp_pdf_promoters
        File? sp_pdf_h_promoters = pc_peaseq.sp_pdf_h_promoters
        File? sp_png_h_promoters = pc_peaseq.sp_png_h_promoters
        File? sp_jpg_h_promoters = pc_peaseq.sp_jpg_h_promoters
        File? peak_promoters = if (paired_sample_m && paired_control_m) then pc_peaseq.peak_promoters
            else sc_seaseq.peak_promoters
        File? peak_genebody = if (paired_sample_m && paired_control_m) then pc_peaseq.peak_genebody
            else sc_seaseq.peak_genebody
        File? peak_window = if (paired_sample_m && paired_control_m) then pc_peaseq.peak_window
            else sc_seaseq.peak_window
        File? peak_closest = if (paired_sample_m && paired_control_m) then pc_peaseq.peak_closest
            else sc_seaseq.peak_closest
        File? peak_comparison = if (paired_sample_m && paired_control_m) then pc_peaseq.peak_comparison
            else sc_seaseq.peak_comparison
        File? gene_comparison = if (paired_sample_m && paired_control_m) then pc_peaseq.gene_comparison
            else sc_seaseq.gene_comparison
        File? pdf_comparison = if (paired_sample_m && paired_control_m) then pc_peaseq.pdf_comparison
            else sc_seaseq.pdf_comparison
        File? all_peak_promoters = if (paired_sample_m && paired_control_m) then pc_peaseq.all_peak_promoters
            else sc_seaseq.all_peak_promoters
        File? all_peak_genebody = if (paired_sample_m && paired_control_m) then pc_peaseq.all_peak_genebody
            else sc_seaseq.all_peak_genebody
        File? all_peak_window = if (paired_sample_m && paired_control_m) then pc_peaseq.all_peak_window
            else sc_seaseq.all_peak_window
        File? all_peak_closest = if (paired_sample_m && paired_control_m) then pc_peaseq.all_peak_closest
            else sc_seaseq.all_peak_closest
        File? all_peak_comparison = if (paired_sample_m && paired_control_m) then pc_peaseq.all_peak_comparison
            else sc_seaseq.all_peak_comparison
        File? all_gene_comparison = if (paired_sample_m && paired_control_m) then pc_peaseq.all_gene_comparison
            else sc_seaseq.all_gene_comparison
        File? all_pdf_comparison = if (paired_sample_m && paired_control_m) then pc_peaseq.all_pdf_comparison
            else sc_seaseq.all_pdf_comparison
        File? nomodel_peak_promoters = if (paired_sample_m && paired_control_m) then pc_peaseq.nomodel_peak_promoters
            else sc_seaseq.nomodel_peak_promoters
        File? nomodel_peak_genebody = if (paired_sample_m && paired_control_m) then pc_peaseq.nomodel_peak_genebody
            else sc_seaseq.nomodel_peak_genebody
        File? nomodel_peak_window = if (paired_sample_m && paired_control_m) then pc_peaseq.nomodel_peak_window
            else sc_seaseq.nomodel_peak_window
        File? nomodel_peak_closest = if (paired_sample_m && paired_control_m) then pc_peaseq.nomodel_peak_closest
            else sc_seaseq.nomodel_peak_closest
        File? nomodel_peak_comparison = if (paired_sample_m && paired_control_m) then pc_peaseq.nomodel_peak_comparison
            else sc_seaseq.nomodel_peak_comparison
        File? nomodel_gene_comparison = if (paired_sample_m && paired_control_m) then pc_peaseq.nomodel_gene_comparison
            else sc_seaseq.nomodel_gene_comparison
        File? nomodel_pdf_comparison = if (paired_sample_m && paired_control_m) then pc_peaseq.nomodel_pdf_comparison
            else sc_seaseq.nomodel_pdf_comparison
        File? sicer_peak_promoters = if (paired_sample_m && paired_control_m) then pc_peaseq.sicer_peak_promoters
            else sc_seaseq.sicer_peak_promoters
        File? sicer_peak_genebody = if (paired_sample_m && paired_control_m) then pc_peaseq.sicer_peak_genebody
            else sc_seaseq.sicer_peak_genebody
        File? sicer_peak_window = if (paired_sample_m && paired_control_m) then pc_peaseq.sicer_peak_window
            else sc_seaseq.sicer_peak_window
        File? sicer_peak_closest = if (paired_sample_m && paired_control_m) then pc_peaseq.sicer_peak_closest
            else sc_seaseq.sicer_peak_closest
        File? sicer_peak_comparison = if (paired_sample_m && paired_control_m) then pc_peaseq.sicer_peak_comparison
            else sc_seaseq.sicer_peak_comparison
        File? sicer_gene_comparison = if (paired_sample_m && paired_control_m) then pc_peaseq.sicer_gene_comparison
            else sc_seaseq.sicer_gene_comparison
        File? sicer_pdf_comparison = if (paired_sample_m && paired_control_m) then pc_peaseq.sicer_pdf_comparison
            else sc_seaseq.sicer_pdf_comparison
        File? sp_peak_promoters = pc_peaseq.sp_peak_promoters
        File? sp_peak_genebody = pc_peaseq.sp_peak_genebody
        File? sp_peak_window = pc_peaseq.sp_peak_window
        File? sp_peak_closest = pc_peaseq.sp_peak_closest
        File? sp_peak_comparison = pc_peaseq.sp_peak_comparison
        File? sp_gene_comparison = pc_peaseq.sp_gene_comparison
        File? sp_pdf_comparison = pc_peaseq.sp_pdf_comparison
        File? sp_all_peak_promoters = pc_peaseq.sp_all_peak_promoters
        File? sp_all_peak_genebody = pc_peaseq.sp_all_peak_genebody
        File? sp_all_peak_window = pc_peaseq.sp_all_peak_window
        File? sp_all_peak_closest = pc_peaseq.sp_all_peak_closest
        File? sp_all_peak_comparison = pc_peaseq.sp_all_peak_comparison
        File? sp_all_gene_comparison = pc_peaseq.sp_all_gene_comparison
        File? sp_all_pdf_comparison = pc_peaseq.sp_all_pdf_comparison
        File? sp_nomodel_peak_promoters = pc_peaseq.sp_nomodel_peak_promoters
        File? sp_nomodel_peak_genebody = pc_peaseq.sp_nomodel_peak_genebody
        File? sp_nomodel_peak_window = pc_peaseq.sp_nomodel_peak_window
        File? sp_nomodel_peak_closest = pc_peaseq.sp_nomodel_peak_closest
        File? sp_nomodel_peak_comparison = pc_peaseq.sp_nomodel_peak_comparison
        File? sp_nomodel_gene_comparison = pc_peaseq.sp_nomodel_gene_comparison
        File? sp_nomodel_pdf_comparison = pc_peaseq.sp_nomodel_pdf_comparison
        File? sp_sicer_peak_promoters = pc_peaseq.sp_sicer_peak_promoters
        File? sp_sicer_peak_genebody = pc_peaseq.sp_sicer_peak_genebody
        File? sp_sicer_peak_window = pc_peaseq.sp_sicer_peak_window
        File? sp_sicer_peak_closest = pc_peaseq.sp_sicer_peak_closest
        File? sp_sicer_peak_comparison = pc_peaseq.sp_sicer_peak_comparison
        File? sp_sicer_gene_comparison = pc_peaseq.sp_sicer_gene_comparison
        File? sp_sicer_pdf_comparison = pc_peaseq.sp_sicer_pdf_comparison
        File? bigwig = if (paired_sample_m && paired_control_m) then pc_peaseq.bigwig else sc_seaseq.bigwig
        File? norm_wig = if (paired_sample_m && paired_control_m) then pc_peaseq.norm_wig else sc_seaseq.norm_wig
        File? tdffile = if (paired_sample_m && paired_control_m) then pc_peaseq.tdffile else sc_seaseq.tdffile
        File? n_bigwig = if (paired_sample_m && paired_control_m) then pc_peaseq.n_bigwig else sc_seaseq.n_bigwig
        File? n_norm_wig = if (paired_sample_m && paired_control_m) then pc_peaseq.n_norm_wig
            else sc_seaseq.n_norm_wig
        File? n_tdffile = if (paired_sample_m && paired_control_m) then pc_peaseq.n_tdffile else sc_seaseq.n_tdffile
        File? a_bigwig = if (paired_sample_m && paired_control_m) then pc_peaseq.a_bigwig else sc_seaseq.a_bigwig
        File? a_norm_wig = if (paired_sample_m && paired_control_m) then pc_peaseq.a_norm_wig
            else sc_seaseq.a_norm_wig
        File? a_tdffile = if (paired_sample_m && paired_control_m) then pc_peaseq.a_tdffile else sc_seaseq.a_tdffile
        File? c_bigwig = if (paired_sample_m && paired_control_m) then pc_peaseq.c_bigwig else sc_seaseq.c_bigwig
        File? c_norm_wig = if (paired_sample_m && paired_control_m) then pc_peaseq.c_norm_wig
            else sc_seaseq.c_norm_wig
        File? c_tdffile = if (paired_sample_m && paired_control_m) then pc_peaseq.c_tdffile else sc_seaseq.c_tdffile
        File? c_n_bigwig = if (paired_sample_m && paired_control_m) then pc_peaseq.c_n_bigwig
            else sc_seaseq.c_n_bigwig
        File? c_n_norm_wig = if (paired_sample_m && paired_control_m) then pc_peaseq.c_n_norm_wig
            else sc_seaseq.c_n_norm_wig
        File? c_n_tdffile = if (paired_sample_m && paired_control_m) then pc_peaseq.c_n_tdffile
            else sc_seaseq.c_n_tdffile
        File? c_a_bigwig = if (paired_sample_m && paired_control_m) then pc_peaseq.c_a_bigwig
            else sc_seaseq.c_a_bigwig
        File? c_a_norm_wig = if (paired_sample_m && paired_control_m) then pc_peaseq.c_a_norm_wig
            else sc_seaseq.c_a_norm_wig
        File? c_a_tdffile = if (paired_sample_m && paired_control_m) then pc_peaseq.c_a_tdffile
            else sc_seaseq.c_a_tdffile
        File? s_bigwig = if (paired_sample_m && paired_control_m) then pc_peaseq.s_bigwig else sc_seaseq.s_bigwig
        File? s_norm_wig = if (paired_sample_m && paired_control_m) then pc_peaseq.s_norm_wig
            else sc_seaseq.s_norm_wig
        File? s_tdffile = if (paired_sample_m && paired_control_m) then pc_peaseq.s_tdffile else sc_seaseq.s_tdffile
        File? sp_bigwig = pc_peaseq.sp_bigwig
        File? sp_norm_wig = pc_peaseq.sp_norm_wig
        File? sp_tdffile = pc_peaseq.sp_tdffile
        File? sp_n_bigwig = pc_peaseq.sp_n_bigwig
        File? sp_n_norm_wig = pc_peaseq.sp_n_norm_wig
        File? sp_n_tdffile = pc_peaseq.sp_n_tdffile
        File? sp_a_bigwig = pc_peaseq.sp_a_bigwig
        File? sp_a_norm_wig = pc_peaseq.sp_a_norm_wig
        File? sp_a_tdffile = pc_peaseq.sp_a_tdffile
        File? cp_bigwig = pc_peaseq.cp_bigwig
        File? cp_norm_wig = pc_peaseq.cp_norm_wig
        File? cp_tdffile = pc_peaseq.cp_tdffile
        File? cp_n_bigwig = pc_peaseq.cp_n_bigwig
        File? cp_n_norm_wig = pc_peaseq.cp_n_norm_wig
        File? cp_n_tdffile = pc_peaseq.cp_n_tdffile
        File? cp_a_bigwig = pc_peaseq.cp_a_bigwig
        File? cp_a_norm_wig = pc_peaseq.cp_a_norm_wig
        File? cp_a_tdffile = pc_peaseq.cp_a_tdffile
        File? sp_s_bigwig = pc_peaseq.sp_s_bigwig
        File? sp_s_norm_wig = pc_peaseq.sp_s_norm_wig
        File? sp_s_tdffile = pc_peaseq.sp_s_tdffile
        File? sf_bigwig = pc_peaseq.sf_bigwig
        File? sf_tdffile = pc_peaseq.sf_tdffile
        File? sf_wigfile = pc_peaseq.sf_wigfile
        File? cf_bigwig = pc_peaseq.cf_bigwig
        File? cf_tdffile = pc_peaseq.cf_tdffile
        File? cf_wigfile = pc_peaseq.cf_wigfile
        Array[File?]? s_qc_statsfile = if (paired_sample_m && paired_control_m) then pc_peaseq.s_qc_statsfile
            else sc_seaseq.s_qc_statsfile
        Array[File?]? s_qc_htmlfile = if (paired_sample_m && paired_control_m) then pc_peaseq.s_qc_htmlfile
            else sc_seaseq.s_qc_htmlfile
        Array[File?]? s_qc_textfile = if (paired_sample_m && paired_control_m) then pc_peaseq.s_qc_textfile
            else sc_seaseq.s_qc_textfile
        File? s_qc_mergehtml = if (paired_sample_m && paired_control_m) then pc_peaseq.s_qc_mergehtml
            else sc_seaseq.s_qc_mergehtml
        Array[File?]? c_qc_statsfile = if (paired_sample_m && paired_control_m) then pc_peaseq.c_qc_statsfile
            else sc_seaseq.c_qc_statsfile
        Array[File?]? c_qc_htmlfile = if (paired_sample_m && paired_control_m) then pc_peaseq.c_qc_htmlfile
            else sc_seaseq.c_qc_htmlfile
        Array[File?]? c_qc_textfile = if (paired_sample_m && paired_control_m) then pc_peaseq.c_qc_textfile
            else sc_seaseq.c_qc_textfile
        File? c_qc_mergehtml = if (paired_sample_m && paired_control_m) then pc_peaseq.c_qc_mergehtml
            else sc_seaseq.c_qc_mergehtml
        Array[File?]? sp_qc_statsfile = pc_peaseq.sp_qc_statsfile
        Array[File?]? sp_qc_htmlfile = pc_peaseq.sp_qc_htmlfile
        Array[File?]? sp_qc_textfile = pc_peaseq.sp_qc_textfile
        Array[File?]? cp_qc_statsfile = pc_peaseq.cp_qc_statsfile
        Array[File?]? cp_qc_htmlfile = pc_peaseq.cp_qc_htmlfile
        Array[File?]? cp_qc_textfile = pc_peaseq.cp_qc_textfile
        File? s_uno_statsfile = if (paired_sample_m && paired_control_m) then pc_peaseq.s_uno_statsfile
            else sc_seaseq.s_uno_statsfile
        File? s_uno_htmlfile = if (paired_sample_m && paired_control_m) then pc_peaseq.s_uno_htmlfile
            else sc_seaseq.s_uno_htmlfile
        File? s_uno_textfile = if (paired_sample_m && paired_control_m) then pc_peaseq.s_uno_textfile
            else sc_seaseq.s_uno_textfile
        File? c_uno_statsfile = if (paired_sample_m && paired_control_m) then pc_peaseq.c_uno_statsfile
            else sc_seaseq.c_uno_statsfile
        File? c_uno_htmlfile = if (paired_sample_m && paired_control_m) then pc_peaseq.c_uno_htmlfile
            else sc_seaseq.c_uno_htmlfile
        File? c_uno_textfile = if (paired_sample_m && paired_control_m) then pc_peaseq.c_uno_textfile
            else sc_seaseq.c_uno_textfile
        File? statsfile = if (paired_sample_m && paired_control_m) then pc_peaseq.statsfile else sc_seaseq.statsfile
        File? htmlfile = if (paired_sample_m && paired_control_m) then pc_peaseq.htmlfile else sc_seaseq.htmlfile
        File? textfile = if (paired_sample_m && paired_control_m) then pc_peaseq.textfile else sc_seaseq.textfile
        File? s_statsfile = pc_peaseq.s_statsfile
        File? s_htmlfile = pc_peaseq.s_htmlfile
        File? s_textfile = pc_peaseq.s_textfile
        File? c_statsfile = pc_peaseq.c_statsfile
        File? c_htmlfile = pc_peaseq.c_htmlfile
        File? c_textfile = pc_peaseq.c_textfile
        File? sp_statsfile = pc_peaseq.sp_statsfile
        File? sp_htmlfile = pc_peaseq.sp_htmlfile
        File? sp_textfile = pc_peaseq.sp_textfile
        File? s_summaryhtml = pc_peaseq.s_summaryhtml
        File? s_summarystats = pc_peaseq.s_summarystats
        File? s_summarytxt = pc_peaseq.s_summarytxt
        File? c_summaryhtml = pc_peaseq.c_summaryhtml
        File? c_summarystats = pc_peaseq.c_summarystats
        File? c_summarytxt = pc_peaseq.c_summarytxt
        File? s_fragsize = pc_peaseq.s_fragsize
        File? c_fragsize = pc_peaseq.c_fragsize
        File? summaryhtml = if (paired_sample_m && paired_control_m) then pc_peaseq.summaryhtml
            else sc_seaseq.summaryhtml
        File? summarytxt = if (paired_sample_m && paired_control_m) then pc_peaseq.summarytxt
            else sc_seaseq.summarytxt
        File? sp_summaryhtml = pc_peaseq.sp_summaryhtml
        File? sp_summarytxt = pc_peaseq.sp_summarytxt
    }
}
