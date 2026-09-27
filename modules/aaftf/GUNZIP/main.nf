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
    tuple val(sample), path("out.fasta"), emit: plain

    script:
    """
    gunzip -c ${gzfile} > out.fasta
    """

    stub:
    """
    gunzip -c ${gzfile} > out.fasta
    """
}
