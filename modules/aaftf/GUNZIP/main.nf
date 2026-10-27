// Generic decompress: takes a gzip-compressed file and emits the plain
// version, for stages that need to feed a durably-published .gz FASTA
// (reused assembly/vecscreen output, or batched contam_clean output) back
// into an AAFTF subcommand — those only accept plain (uncompressed) FASTA.
process GUNZIP {
    tag   { sample }
    label 'aaftf_lite'

    input:
    tuple val(sample), path(gzfile)

    output:
    // Name the output per sample. CONTAM_CLEAN_BATCH stages every sample's
    // FASTA into one work dir, so a shared name such as out.fasta collides.
    tuple val(sample), path("${sample}.fasta"), emit: plain

    script:
    """
    gunzip -c ${gzfile} > ${sample}.fasta
    """

    stub:
    """
    gunzip -c ${gzfile} > ${sample}.fasta
    """
}
