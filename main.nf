nextflow.enable.dsl=2

// Parameters
params.accession = 'M21012'
params.in        = null
params.out       = 'results'

process FETCH_REFERENCE {
    tag { accession }
    conda 'bioconda::entrez-direct=24.0'

    input:
    val accession

    output:
    path "reference.fasta"

    script:
    """
    esearch -db nucleotide -query "${accession}" | \
    efetch -format fasta > reference.fasta
    """
}

// Process_1
process COMBINE_SEQUENCES {
    tag "merging"

    input:
    path ref_fasta
    path "genomes/*" 

    output:
    path "combined.fasta"

    script:
    """
    cat ${ref_fasta} > combined.fasta
    cat genomes/* >> combined.fasta
    """
}

process MAFFT_ALIGN {
    tag "mafft"
    conda 'bioconda::mafft=7.525'

    input:
    path combined_fasta

    output:
    path "alignment.fasta"

    script:
    """
    mafft --auto --inputorder "${combined_fasta}" > alignment.fasta
    """
}

process TRIMAL_CLEAN {
    tag "trimal"
    conda 'bioconda::trimal=1.5.0'
    publishDir "${params.out}", mode: 'copy'

    input:
    path alignment_fasta

    output:
    path "alignment_trimmed.fasta"
    path "alignment_trimmed.html"

    script:
    """
    trimal -in "${alignment_fasta}" \
           -out alignment_trimmed.fasta \
           -automated1 \
           -htmlout alignment_trimmed.html
    """
}

workflow {
    if (!params.in) error 'Missing --in (e.g. --in "*.fasta")'

    // 1. Fetch Reference
    ch_ref = FETCH_REFERENCE(params.accession)

    // 2. Collect and Sort genomes into a single List
    ch_genomes_sorted = Channel
        .fromPath(params.in, checkIfExists: true)
        .toSortedList({ a, b -> a.name <=> b.name })

    // 3. Combine: Pass the single ref and the list of genomes
    ch_combined = COMBINE_SEQUENCES(ch_ref, ch_genomes_sorted)

    // 4. Align and Trim
    ch_align = MAFFT_ALIGN(ch_combined)
    TRIMAL_CLEAN(ch_align)
}


