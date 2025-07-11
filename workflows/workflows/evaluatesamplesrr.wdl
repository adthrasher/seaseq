version 1.0

import "../../peaseq-case.wdl" as ps
import "../../seaseq-case.wdl" as ss
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
        String? results_name
        Boolean run_motifs = true
        Int insertsize = 600
        String strandedness = "fr"
        Boolean paired = true
    }

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
    Boolean paired_sample_m = select_first([
        paired_sample[0],
        paired,
    ])

    if (!paired_sample_m) {
        call ss.seaseq as ss_seaseq { input:
            reference = reference,
            spikein_reference = spikein_reference,
            blacklist = blacklist,
            gtf = gtf,
            bowtie_index = bowtie_index,
            spikein_bowtie_index = spikein_bowtie_index,
            motif_databases = motif_databases,
            sample_fastq = sample_R1_fastq,
            sample_sraid = sample_sraid,
            results_name = results_name,
            run_motifs = run_motifs,
        }
    }

    if (paired_sample_m) {
        call ps.peaseq as ps_peaseq { input:
            reference = reference,
            spikein_reference = spikein_reference,
            blacklist = blacklist,
            gtf = gtf,
            bowtie_index = bowtie_index,
            spikein_bowtie_index = spikein_bowtie_index,
            motif_databases = motif_databases,
            sample_R1_fastq = sample_R1_fastq,
            sample_R2_fastq = sample_R2_fastq,
            sample_sraid = sample_sraid,
            insertsize = insertsize,
            strandedness = strandedness,
            results_name = results_name,
            run_motifs = run_motifs,
        }
    }

    # Processing OUTPUTs
    output {
        Array[File?]? spikein_indv_s_htmlfile = if paired_sample_m then ps_peaseq.spikein_indv_s_htmlfile
            else ss_seaseq.spikein_indv_s_htmlfile
        Array[File?]? spikein_indv_s_zipfile = if paired_sample_m then ps_peaseq.spikein_indv_s_zipfile
            else ss_seaseq.spikein_indv_s_zipfile
        Array[File?]? spikein_s_metrics_out = if paired_sample_m then ps_peaseq.spikein_s_metrics_out
            else ss_seaseq.spikein_s_metrics_out
        Array[File?]? indv_s_htmlfile = if paired_sample_m then ps_peaseq.indv_s_htmlfile else ss_seaseq.indv_s_htmlfile
        Array[File?]? indv_s_zipfile = if paired_sample_m then ps_peaseq.indv_s_zipfile else ss_seaseq.indv_s_zipfile
        Array[File?]? indv_s_bam_htmlfile = if paired_sample_m then ps_peaseq.indv_s_bam_htmlfile
            else ss_seaseq.indv_s_bam_htmlfile
        Array[File?]? indv_s_bam_zipfile = if paired_sample_m then ps_peaseq.indv_s_bam_zipfile
            else ss_seaseq.indv_s_bam_zipfile
        File? s_mergebam_htmlfile = if paired_sample_m then ps_peaseq.s_mergebam_htmlfile else ss_seaseq.s_mergebam_htmlfile
        File? s_mergebam_zipfile = if paired_sample_m then ps_peaseq.s_mergebam_zipfile else ss_seaseq.s_mergebam_zipfile
        Array[File?]? indv_sp_bam_htmlfile = ps_peaseq.indv_sp_bam_htmlfile
        Array[File?]? indv_sp_bam_zipfile = ps_peaseq.indv_sp_bam_zipfile
        File? sp_mergebam_htmlfile = ps_peaseq.sp_mergebam_htmlfile
        File? sp_mergebam_zipfile = ps_peaseq.sp_mergebam_zipfile
        File? uno_s_htmlfile = ss_seaseq.uno_s_htmlfile
        File? uno_s_zipfile = ss_seaseq.uno_s_zipfile
        File? uno_s_bam_htmlfile = if paired_sample_m then ps_peaseq.uno_s_bam_htmlfile else ss_seaseq.uno_s_bam_htmlfile
        File? uno_s_bam_zipfile = if paired_sample_m then ps_peaseq.uno_s_bam_zipfile else ss_seaseq.uno_s_bam_zipfile
        Array[File?]? s_metrics_out = if paired_sample_m then ps_peaseq.s_metrics_out else ss_seaseq.s_metrics_out
        File? uno_s_metrics_out = ss_seaseq.uno_s_metrics_out
        Array[File?]? indv_s_sortedbam = if paired_sample_m then ps_peaseq.indv_s_sortedbam else ss_seaseq.indv_s_sortedbam
        Array[File?]? indv_s_indexbam = if paired_sample_m then ps_peaseq.indv_s_indexbam else ss_seaseq.indv_s_indexbam
        Array[File?]? indv_s_bkbam = if paired_sample_m then ps_peaseq.indv_s_bkbam else ss_seaseq.indv_s_bkbam
        Array[File?]? indv_s_bkindexbam = if paired_sample_m then ps_peaseq.indv_s_bkindexbam
            else ss_seaseq.indv_s_bkindexbam
        Array[File?]? indv_s_rmbam = if paired_sample_m then ps_peaseq.indv_s_rmbam else ss_seaseq.indv_s_rmbam
        Array[File?]? indv_s_rmindexbam = if paired_sample_m then ps_peaseq.indv_s_rmindexbam
            else ss_seaseq.indv_s_rmindexbam
        Array[File?]? indv_sp_sortedbam = ps_peaseq.indv_sp_sortedbam
        Array[File?]? indv_sp_indexbam = ps_peaseq.indv_sp_indexbam
        Array[File?]? indv_sp_bkbam = ps_peaseq.indv_sp_bkbam
        Array[File?]? indv_sp_bkindexbam = ps_peaseq.indv_sp_bkindexbam
        Array[File?]? indv_sp_rmbam = ps_peaseq.indv_sp_rmbam
        Array[File?]? indv_sp_rmindexbam = ps_peaseq.indv_sp_rmindexbam
        File? uno_s_sortedbam = if paired_sample_m then ps_peaseq.uno_s_sortedbam else ss_seaseq.uno_s_sortedbam
        File? uno_s_indexstatsbam = if paired_sample_m then ps_peaseq.uno_s_indexstatsbam else ss_seaseq.uno_s_indexstatsbam
        File? uno_s_bkbam = if paired_sample_m then ps_peaseq.uno_s_bkbam else ss_seaseq.uno_s_bkbam
        File? uno_s_bkindexbam = if paired_sample_m then ps_peaseq.uno_s_bkindexbam else ss_seaseq.uno_s_bkindexbam
        File? uno_s_rmbam = if paired_sample_m then ps_peaseq.uno_s_rmbam else ss_seaseq.uno_s_rmbam
        File? uno_s_rmindexbam = if paired_sample_m then ps_peaseq.uno_s_rmindexbam else ss_seaseq.uno_s_rmindexbam
        File? s_mergebamfile = if paired_sample_m then ps_peaseq.s_mergebamfile else ss_seaseq.s_mergebamfile
        File? s_mergebamindex = if paired_sample_m then ps_peaseq.s_mergebamindex else ss_seaseq.s_mergebamindex
        File? s_bkbam = if paired_sample_m then ps_peaseq.s_bkbam else ss_seaseq.s_bkbam
        File? s_bkindexbam = if paired_sample_m then ps_peaseq.s_bkindexbam else ss_seaseq.s_bkindexbam
        File? s_rmbam = if paired_sample_m then ps_peaseq.s_rmbam else ss_seaseq.s_rmbam
        File? s_rmindexbam = if paired_sample_m then ps_peaseq.s_rmindexbam else ss_seaseq.s_rmindexbam
        File? sp_mergebamfile = ps_peaseq.sp_mergebamfile
        File? sp_mergebamindex = ps_peaseq.sp_mergebamindex
        File? sp_bkbam = ps_peaseq.sp_bkbam
        File? sp_bkindexbam = ps_peaseq.sp_bkindexbam
        File? sp_rmbam = ps_peaseq.sp_rmbam
        File? sp_rmindexbam = ps_peaseq.sp_rmindexbam
        File? s_fragments_bam = ps_peaseq.s_fragments_bam
        File? s_fragments_indexbam = ps_peaseq.s_fragments_indexbam
        File? peakbedfile = if paired_sample_m then ps_peaseq.peakbedfile else ss_seaseq.peakbedfile
        File? peakxlsfile = if paired_sample_m then ps_peaseq.peakxlsfile else ss_seaseq.peakxlsfile
        File? negativexlsfile = if paired_sample_m then ps_peaseq.negativexlsfile else ss_seaseq.negativexlsfile
        File? summitsfile = if paired_sample_m then ps_peaseq.summitsfile else ss_seaseq.summitsfile
        File? wigfile = if paired_sample_m then ps_peaseq.wigfile else ss_seaseq.wigfile
        File? all_peakbedfile = if paired_sample_m then ps_peaseq.all_peakbedfile else ss_seaseq.all_peakbedfile
        File? all_peakxlsfile = if paired_sample_m then ps_peaseq.all_peakxlsfile else ss_seaseq.all_peakxlsfile
        File? all_negativexlsfile = if paired_sample_m then ps_peaseq.all_negativexlsfile else ss_seaseq.all_negativexlsfile
        File? all_summitsfile = if paired_sample_m then ps_peaseq.all_summitsfile else ss_seaseq.all_summitsfile
        File? all_wigfile = if paired_sample_m then ps_peaseq.all_wigfile else ss_seaseq.all_wigfile
        File? nm_peakbedfile = if paired_sample_m then ps_peaseq.nm_peakbedfile else ss_seaseq.nm_peakbedfile
        File? nm_peakxlsfile = if paired_sample_m then ps_peaseq.nm_peakxlsfile else ss_seaseq.nm_peakxlsfile
        File? nm_negativexlsfile = if paired_sample_m then ps_peaseq.nm_negativexlsfile else ss_seaseq.nm_negativexlsfile
        File? nm_summitsfile = if paired_sample_m then ps_peaseq.nm_summitsfile else ss_seaseq.nm_summitsfile
        File? nm_wigfile = if paired_sample_m then ps_peaseq.nm_wigfile else ss_seaseq.nm_wigfile
        File? readme_peaks = if paired_sample_m then ps_peaseq.readme_peaks else ss_seaseq.readme_peaks
        File? sp_peakbedfile = ps_peaseq.sp_peakbedfile
        File? sp_peakxlsfile = ps_peaseq.sp_peakxlsfile
        File? sp_negativexlsfile = ps_peaseq.sp_negativexlsfile
        File? sp_summitsfile = ps_peaseq.sp_summitsfile
        File? sp_wigfile = ps_peaseq.sp_wigfile
        File? sp_all_peakbedfile = ps_peaseq.sp_all_peakbedfile
        File? sp_all_peakxlsfile = ps_peaseq.sp_all_peakxlsfile
        File? sp_all_negativexlsfile = ps_peaseq.sp_all_negativexlsfile
        File? sp_all_summitsfile = ps_peaseq.sp_all_summitsfile
        File? sp_all_wigfile = ps_peaseq.sp_all_wigfile
        File? sp_nm_peakbedfile = ps_peaseq.sp_nm_peakbedfile
        File? sp_nm_peakxlsfile = ps_peaseq.sp_nm_peakxlsfile
        File? sp_nm_negativexlsfile = ps_peaseq.sp_nm_negativexlsfile
        File? sp_nm_summitsfile = ps_peaseq.sp_nm_summitsfile
        File? sp_nm_wigfile = ps_peaseq.sp_nm_wigfile
        File? sp_readme_peaks = ps_peaseq.sp_readme_peaks
        File? scoreisland = if paired_sample_m then ps_peaseq.scoreisland else ss_seaseq.scoreisland
        File? sicer_wigfile = if paired_sample_m then ps_peaseq.sicer_wigfile else ss_seaseq.sicer_wigfile
        File? sp_scoreisland = ps_peaseq.sp_scoreisland
        File? sp_sicer_wigfile = ps_peaseq.sp_sicer_wigfile
        File? pngfile = if paired_sample_m then ps_peaseq.pngfile else ss_seaseq.pngfile
        File? mapped_union = if paired_sample_m then ps_peaseq.mapped_union else ss_seaseq.mapped_union
        File? mapped_stitch = if paired_sample_m then ps_peaseq.mapped_stitch else ss_seaseq.mapped_stitch
        File? enhancers = if paired_sample_m then ps_peaseq.enhancers else ss_seaseq.enhancers
        File? super_enhancers = if paired_sample_m then ps_peaseq.super_enhancers else ss_seaseq.super_enhancers
        File? gff_file = if paired_sample_m then ps_peaseq.gff_file else ss_seaseq.gff_file
        File? gff_union = if paired_sample_m then ps_peaseq.gff_union else ss_seaseq.gff_union
        File? union_enhancers = if paired_sample_m then ps_peaseq.union_enhancers else ss_seaseq.union_enhancers
        File? stitch_enhancers = if paired_sample_m then ps_peaseq.stitch_enhancers else ss_seaseq.stitch_enhancers
        File? e_to_g_enhancers = if paired_sample_m then ps_peaseq.e_to_g_enhancers else ss_seaseq.e_to_g_enhancers
        File? g_to_e_enhancers = if paired_sample_m then ps_peaseq.g_to_e_enhancers else ss_seaseq.g_to_e_enhancers
        File? e_to_g_super_enhancers = if paired_sample_m then ps_peaseq.e_to_g_super_enhancers
            else ss_seaseq.e_to_g_super_enhancers
        File? g_to_e_super_enhancers = if paired_sample_m then ps_peaseq.g_to_e_super_enhancers
            else ss_seaseq.g_to_e_super_enhancers
        File? sp_pngfile = ps_peaseq.sp_pngfile
        File? sp_mapped_union = ps_peaseq.sp_mapped_union
        File? sp_mapped_stitch = ps_peaseq.sp_mapped_stitch
        File? sp_enhancers = ps_peaseq.sp_enhancers
        File? sp_super_enhancers = ps_peaseq.sp_super_enhancers
        File? sp_gff_file = ps_peaseq.sp_gff_file
        File? sp_gff_union = ps_peaseq.sp_gff_union
        File? sp_union_enhancers = ps_peaseq.sp_union_enhancers
        File? sp_stitch_enhancers = ps_peaseq.sp_stitch_enhancers
        File? sp_e_to_g_enhancers = ps_peaseq.sp_e_to_g_enhancers
        File? sp_g_to_e_enhancers = ps_peaseq.sp_g_to_e_enhancers
        File? sp_e_to_g_super_enhancers = ps_peaseq.sp_e_to_g_super_enhancers
        File? sp_g_to_e_super_enhancers = ps_peaseq.sp_g_to_e_super_enhancers
        File? flankbedfile = if paired_sample_m then ps_peaseq.flankbedfile else ss_seaseq.flankbedfile
        File? ame_tsv = if paired_sample_m then ps_peaseq.ame_tsv else ss_seaseq.ame_tsv
        File? ame_html = if paired_sample_m then ps_peaseq.ame_html else ss_seaseq.ame_html
        File? ame_seq = if paired_sample_m then ps_peaseq.ame_seq else ss_seaseq.ame_seq
        File? meme = if paired_sample_m then ps_peaseq.meme else ss_seaseq.meme
        File? meme_summary = if paired_sample_m then ps_peaseq.meme_summary else ss_seaseq.meme_summary
        File? summit_ame_tsv = if paired_sample_m then ps_peaseq.summit_ame_tsv else ss_seaseq.summit_ame_tsv
        File? summit_ame_html = if paired_sample_m then ps_peaseq.summit_ame_html else ss_seaseq.summit_ame_html
        File? summit_ame_seq = if paired_sample_m then ps_peaseq.summit_ame_seq else ss_seaseq.summit_ame_seq
        File? summit_meme = if paired_sample_m then ps_peaseq.summit_meme else ss_seaseq.summit_meme
        File? summit_meme_summary = if paired_sample_m then ps_peaseq.summit_meme_summary else ss_seaseq.summit_meme_summary
        File? sp_flankbedfile = ps_peaseq.summit_meme_summary
        File? sp_ame_tsv = ps_peaseq.sp_ame_tsv
        File? sp_ame_html = ps_peaseq.sp_ame_html
        File? sp_ame_seq = ps_peaseq.sp_ame_seq
        File? sp_meme = ps_peaseq.sp_meme
        File? sp_meme_summary = ps_peaseq.sp_meme_summary
        File? sp_summit_ame_tsv = ps_peaseq.sp_summit_ame_tsv
        File? sp_summit_ame_html = ps_peaseq.sp_summit_ame_html
        File? sp_summit_ame_seq = ps_peaseq.sp_summit_ame_seq
        File? sp_summit_meme = ps_peaseq.sp_summit_meme
        File? sp_summit_meme_summary = ps_peaseq.sp_summit_meme_summary
        File? s_matrices = if paired_sample_m then ps_peaseq.s_matrices else ss_seaseq.s_matrices
        File? densityplot = if paired_sample_m then ps_peaseq.densityplot else ss_seaseq.densityplot
        File? pdf_gene = if paired_sample_m then ps_peaseq.pdf_gene else ss_seaseq.pdf_gene
        File? pdf_h_gene = if paired_sample_m then ps_peaseq.pdf_h_gene else ss_seaseq.pdf_h_gene
        File? png_h_gene = if paired_sample_m then ps_peaseq.png_h_gene else ss_seaseq.png_h_gene
        File? jpg_h_gene = if paired_sample_m then ps_peaseq.jpg_h_gene else ss_seaseq.jpg_h_gene
        File? pdf_promoters = if paired_sample_m then ps_peaseq.pdf_promoters else ss_seaseq.pdf_promoters
        File? pdf_h_promoters = if paired_sample_m then ps_peaseq.pdf_h_promoters else ss_seaseq.pdf_h_promoters
        File? png_h_promoters = if paired_sample_m then ps_peaseq.png_h_promoters else ss_seaseq.png_h_promoters
        File? jpg_h_promoters = if paired_sample_m then ps_peaseq.jpg_h_promoters else ss_seaseq.jpg_h_promoters
        File? sp_s_matrices = ps_peaseq.sp_s_matrices
        File? sp_densityplot = ps_peaseq.sp_densityplot
        File? sp_pdf_gene = ps_peaseq.sp_pdf_gene
        File? sp_pdf_h_gene = ps_peaseq.sp_pdf_h_gene
        File? sp_png_h_gene = ps_peaseq.sp_png_h_gene
        File? sp_jpg_h_gene = ps_peaseq.sp_jpg_h_gene
        File? sp_pdf_promoters = ps_peaseq.sp_pdf_promoters
        File? sp_pdf_h_promoters = ps_peaseq.sp_pdf_h_promoters
        File? sp_png_h_promoters = ps_peaseq.sp_png_h_promoters
        File? sp_jpg_h_promoters = ps_peaseq.sp_jpg_h_promoters
        File? peak_promoters = if paired_sample_m then ps_peaseq.peak_promoters else ss_seaseq.peak_promoters
        File? peak_genebody = if paired_sample_m then ps_peaseq.peak_genebody else ss_seaseq.peak_genebody
        File? peak_window = if paired_sample_m then ps_peaseq.peak_window else ss_seaseq.peak_window
        File? peak_closest = if paired_sample_m then ps_peaseq.peak_closest else ss_seaseq.peak_closest
        File? peak_comparison = if paired_sample_m then ps_peaseq.peak_comparison else ss_seaseq.peak_comparison
        File? gene_comparison = if paired_sample_m then ps_peaseq.gene_comparison else ss_seaseq.gene_comparison
        File? pdf_comparison = if paired_sample_m then ps_peaseq.pdf_comparison else ss_seaseq.pdf_comparison
        File? all_peak_promoters = if paired_sample_m then ps_peaseq.all_peak_promoters else ss_seaseq.all_peak_promoters
        File? all_peak_genebody = if paired_sample_m then ps_peaseq.all_peak_genebody else ss_seaseq.all_peak_genebody
        File? all_peak_window = if paired_sample_m then ps_peaseq.all_peak_window else ss_seaseq.all_peak_window
        File? all_peak_closest = if paired_sample_m then ps_peaseq.all_peak_closest else ss_seaseq.all_peak_closest
        File? all_peak_comparison = if paired_sample_m then ps_peaseq.all_peak_comparison else ss_seaseq.all_peak_comparison
        File? all_gene_comparison = if paired_sample_m then ps_peaseq.all_gene_comparison else ss_seaseq.all_gene_comparison
        File? all_pdf_comparison = if paired_sample_m then ps_peaseq.all_pdf_comparison else ss_seaseq.all_pdf_comparison
        File? nomodel_peak_promoters = if paired_sample_m then ps_peaseq.nomodel_peak_promoters
            else ss_seaseq.nomodel_peak_promoters
        File? nomodel_peak_genebody = if paired_sample_m then ps_peaseq.nomodel_peak_genebody
            else ss_seaseq.nomodel_peak_genebody
        File? nomodel_peak_window = if paired_sample_m then ps_peaseq.nomodel_peak_window else ss_seaseq.nomodel_peak_window
        File? nomodel_peak_closest = if paired_sample_m then ps_peaseq.nomodel_peak_closest else ss_seaseq.nomodel_peak_closest
        File? nomodel_peak_comparison = if paired_sample_m then ps_peaseq.nomodel_peak_comparison
            else ss_seaseq.nomodel_peak_comparison
        File? nomodel_gene_comparison = if paired_sample_m then ps_peaseq.nomodel_gene_comparison
            else ss_seaseq.nomodel_gene_comparison
        File? nomodel_pdf_comparison = if paired_sample_m then ps_peaseq.nomodel_pdf_comparison
            else ss_seaseq.nomodel_pdf_comparison
        File? sicer_peak_promoters = if paired_sample_m then ps_peaseq.sicer_peak_promoters else ss_seaseq.sicer_peak_promoters
        File? sicer_peak_genebody = if paired_sample_m then ps_peaseq.sicer_peak_genebody else ss_seaseq.sicer_peak_genebody
        File? sicer_peak_window = if paired_sample_m then ps_peaseq.sicer_peak_window else ss_seaseq.sicer_peak_window
        File? sicer_peak_closest = if paired_sample_m then ps_peaseq.sicer_peak_closest else ss_seaseq.sicer_peak_closest
        File? sicer_peak_comparison = if paired_sample_m then ps_peaseq.sicer_peak_comparison
            else ss_seaseq.sicer_peak_comparison
        File? sicer_gene_comparison = if paired_sample_m then ps_peaseq.sicer_gene_comparison
            else ss_seaseq.sicer_gene_comparison
        File? sicer_pdf_comparison = if paired_sample_m then ps_peaseq.sicer_pdf_comparison else ss_seaseq.sicer_pdf_comparison
        File? sp_peak_promoters = ps_peaseq.sp_peak_promoters
        File? sp_peak_genebody = ps_peaseq.sp_peak_genebody
        File? sp_peak_window = ps_peaseq.sp_peak_window
        File? sp_peak_closest = ps_peaseq.sp_peak_closest
        File? sp_peak_comparison = ps_peaseq.sp_peak_comparison
        File? sp_gene_comparison = ps_peaseq.sp_gene_comparison
        File? sp_pdf_comparison = ps_peaseq.sp_pdf_comparison
        File? sp_all_peak_promoters = ps_peaseq.sp_all_peak_promoters
        File? sp_all_peak_genebody = ps_peaseq.sp_all_peak_genebody
        File? sp_all_peak_window = ps_peaseq.sp_all_peak_window
        File? sp_all_peak_closest = ps_peaseq.sp_all_peak_closest
        File? sp_all_peak_comparison = ps_peaseq.sp_all_peak_comparison
        File? sp_all_gene_comparison = ps_peaseq.sp_all_gene_comparison
        File? sp_all_pdf_comparison = ps_peaseq.sp_all_pdf_comparison
        File? sp_nomodel_peak_promoters = ps_peaseq.sp_nomodel_peak_promoters
        File? sp_nomodel_peak_genebody = ps_peaseq.sp_nomodel_peak_genebody
        File? sp_nomodel_peak_window = ps_peaseq.sp_nomodel_peak_window
        File? sp_nomodel_peak_closest = ps_peaseq.sp_nomodel_peak_closest
        File? sp_nomodel_peak_comparison = ps_peaseq.sp_nomodel_peak_comparison
        File? sp_nomodel_gene_comparison = ps_peaseq.sp_nomodel_gene_comparison
        File? sp_nomodel_pdf_comparison = ps_peaseq.sp_nomodel_pdf_comparison
        File? sp_sicer_peak_promoters = ps_peaseq.sp_sicer_peak_promoters
        File? sp_sicer_peak_genebody = ps_peaseq.sp_sicer_peak_genebody
        File? sp_sicer_peak_window = ps_peaseq.sp_sicer_peak_window
        File? sp_sicer_peak_closest = ps_peaseq.sp_sicer_peak_closest
        File? sp_sicer_peak_comparison = ps_peaseq.sp_sicer_peak_comparison
        File? sp_sicer_gene_comparison = ps_peaseq.sp_sicer_gene_comparison
        File? sp_sicer_pdf_comparison = ps_peaseq.sp_sicer_pdf_comparison
        File? bigwig = if paired_sample_m then ps_peaseq.bigwig else ss_seaseq.bigwig
        File? norm_wig = if paired_sample_m then ps_peaseq.norm_wig else ss_seaseq.norm_wig
        File? tdffile = if paired_sample_m then ps_peaseq.tdffile else ss_seaseq.tdffile
        File? n_bigwig = if paired_sample_m then ps_peaseq.n_bigwig else ss_seaseq.n_bigwig
        File? n_norm_wig = if paired_sample_m then ps_peaseq.n_norm_wig else ss_seaseq.n_norm_wig
        File? n_tdffile = if paired_sample_m then ps_peaseq.n_tdffile else ss_seaseq.n_tdffile
        File? a_bigwig = if paired_sample_m then ps_peaseq.a_bigwig else ss_seaseq.a_bigwig
        File? a_norm_wig = if paired_sample_m then ps_peaseq.a_norm_wig else ss_seaseq.a_norm_wig
        File? a_tdffile = if paired_sample_m then ps_peaseq.a_tdffile else ss_seaseq.a_tdffile
        File? s_bigwig = if paired_sample_m then ps_peaseq.s_bigwig else ss_seaseq.s_bigwig
        File? s_norm_wig = if paired_sample_m then ps_peaseq.s_norm_wig else ss_seaseq.s_norm_wig
        File? s_tdffile = if paired_sample_m then ps_peaseq.s_tdffile else ss_seaseq.s_tdffile
        File? sp_bigwig = ps_peaseq.sp_bigwig
        File? sp_norm_wig = ps_peaseq.sp_norm_wig
        File? sp_tdffile = ps_peaseq.sp_tdffile
        File? sp_n_bigwig = ps_peaseq.sp_n_bigwig
        File? sp_n_norm_wig = ps_peaseq.sp_n_norm_wig
        File? sp_n_tdffile = ps_peaseq.sp_n_tdffile
        File? sp_a_bigwig = ps_peaseq.sp_a_bigwig
        File? sp_a_norm_wig = ps_peaseq.sp_a_norm_wig
        File? sp_a_tdffile = ps_peaseq.sp_a_tdffile
        File? sp_s_bigwig = ps_peaseq.sp_s_bigwig
        File? sp_s_norm_wig = ps_peaseq.sp_s_norm_wig
        File? sp_s_tdffile = ps_peaseq.sp_s_tdffile
        File? sf_bigwig = ps_peaseq.sf_bigwig
        File? sf_tdffile = ps_peaseq.sf_tdffile
        File? sf_wigfile = ps_peaseq.sf_wigfile
        Array[File?]? s_qc_statsfile = if paired_sample_m then ps_peaseq.s_qc_statsfile else ss_seaseq.s_qc_statsfile
        Array[File?]? s_qc_htmlfile = if paired_sample_m then ps_peaseq.s_qc_htmlfile else ss_seaseq.s_qc_htmlfile
        Array[File?]? s_qc_textfile = if paired_sample_m then ps_peaseq.s_qc_textfile else ss_seaseq.s_qc_textfile
        File? s_qc_mergehtml = if paired_sample_m then ps_peaseq.s_qc_mergehtml else ss_seaseq.s_qc_mergehtml
        Array[File?]? sp_qc_statsfile = ps_peaseq.sp_qc_statsfile
        Array[File?]? sp_qc_htmlfile = ps_peaseq.sp_qc_htmlfile
        Array[File?]? sp_qc_textfile = ps_peaseq.sp_qc_textfile
        File? s_uno_statsfile = if paired_sample_m then ps_peaseq.s_uno_statsfile else ss_seaseq.s_uno_statsfile
        File? s_uno_htmlfile = if paired_sample_m then ps_peaseq.s_uno_htmlfile else ss_seaseq.s_uno_htmlfile
        File? s_uno_textfile = if paired_sample_m then ps_peaseq.s_uno_textfile else ss_seaseq.s_uno_textfile
        File? statsfile = if paired_sample_m then ps_peaseq.statsfile else ss_seaseq.statsfile
        File? htmlfile = if paired_sample_m then ps_peaseq.htmlfile else ss_seaseq.htmlfile
        File? textfile = if paired_sample_m then ps_peaseq.textfile else ss_seaseq.textfile
        File? sp_statsfile = ps_peaseq.sp_statsfile
        File? sp_htmlfile = ps_peaseq.sp_htmlfile
        File? sp_textfile = ps_peaseq.sp_textfile
        File? s_summaryhtml = ps_peaseq.s_summaryhtml
        File? s_summarytxt = ps_peaseq.s_summarytxt
        File? s_fragsize = ps_peaseq.s_fragsize
        File? summaryhtml = if paired_sample_m then ps_peaseq.summaryhtml else ss_seaseq.summaryhtml
        File? summarytxt = if paired_sample_m then ps_peaseq.summarytxt else ss_seaseq.summarytxt
        File? sp_summaryhtml = ps_peaseq.sp_summaryhtml
        File? sp_summarytxt = ps_peaseq.sp_summarytxt
    }
}
