nextflow.enable.dsl = 2

params.samples = "${projectDir}/ascairn_samples.csv"
params.outdir = "${projectDir}/results"
params.resource_dir = "/scratch/hain/ascairn/ascairn_resource/resource/panel/ascairn_paper_2025"
params.reference = "/scratch/hain/centromere/revision_ascairn_comparison/GRCh38_full_analysis_set_plus_decoy_hla.fa"
params.reference_name = "hg38"
params.cpus = 8
params.ascairn_executable = "ascairn"

process ASCAIRN_TYPE_ALL {
    tag "${sample}"
    cpus params.cpus
    publishDir params.outdir, mode: 'copy', overwrite: false

    input:
    tuple val(sample), val(cram_url)

    output:
    path "${sample}.cen_type_all.txt", emit: results

    script:
    """
    ${params.ascairn_executable} type_all \\
        '${cram_url}' \\
        -o '${sample}' \\
        --resource_dir '${params.resource_dir}' \\
        --reference '${params.reference_name}' \\
        -r '${params.reference}' \\
        -t ${task.cpus}
    """
}

workflow {
    Channel
        .fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(row.SAMPLE_NAME, row.ENA_FILE_PATH) }
        .set { samples }

    ASCAIRN_TYPE_ALL(samples)
}