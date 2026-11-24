nextflow.enable.dsl = 2

params.vcf_file_list = ""
params.map_file = ""
params.nHCs = 150

// pbwt -paintSparse <prefix> <nchunks> <ploidy> <cutoff>:
// number of chunks per region, ploidy, and the minimum match length (in SNPs) for a match to be counted
params.paint_nchunks = 100
params.paint_ploidy = 2
params.paint_cutoff = 500

// combine --restrict: the maximum number of entries in the combined chunk length matrix (smallest values
// are dropped first). The default, 2^31-1, is the most R's Matrix package can hold
params.restrict = 2147483647
// rescale each row of the matrix to its sum before --restrict (and --remove_pairs) removed entries, before computing HCs
params.reweight = false
// optional: a file of sample ID pairs (e.g. relatives; first two whitespace-separated columns) whose entries are
// removed from the matrix before computing HCs. Use with --reweight
params.remove_pairs = ""
// optional: a file of sample IDs (one per line). Only these samples' rows and columns of the matrix are used for the
// SVD (their entries to other samples are removed) and only they get HCs, e.g. HCs within one group from a painting of
// everyone. Use with --reweight
params.fit_samples = ""

// optional: a file of sample IDs (one per line) to keep. The input files are subset to these samples before
// painting, so only these samples are painted (against each other). Output rows follow the order of this file
params.keep_samples = ""

params.outfile = "HCs.csv"
params.diagnostics = "HC_diagnostics.tsv"
// also publish the combined chunk length matrix, its row sums, and the sample list its indices refer to,
// to <resultsdir>/matrix, e.g. to re-run the HC step without repainting
params.publish_matrix = false
params.resultsdir = "results"

workflow {

    // a (channel of a) tuple of the chromosome index and the path to the bcf file for that chromosome,
    // parsed from a plain text file with one line per chromosome, each line containing the chromosome index and the path to the bcf file for that chromosome
    // with columns separated by whitespace. Relative paths are resolved against the directory containing vcf_file_list
    def vcf_file_list = file(params.vcf_file_list)
    c_vcf_input = channel.from(
        vcf_file_list
            .readLines()
        ) |\
        map { line ->
                def fields = line.split("\\h+")
                if (fields.size() != 2) {
                    throw new Exception("Each line in the vcf_file_list must contain exactly two fields: chromosome number and path to vcf/bcf file.")
                }
                return [String.format('%02d', Integer.parseInt(fields[0])), file(vcf_file_list.parent.resolve(fields[1]).toString())]

        }

    // optionally subset every input file to the samples in keep_samples, before anything else
    if (params.keep_samples) {
        c_vcf = c_vcf_input |\
            combine(channel.value(file(params.keep_samples, checkIfExists: true))) |\
            subset_samples
    } else {
        c_vcf = c_vcf_input
    }

    // a (channel of a) plain text file of the sample IDs, one per line, in the order of the input files.
    // check_samples fails if the input files don't all have the same samples in the same order, because
    // the chunk length matrix refers to samples by their position in this list
    c_samples_file = c_vcf |\
            toSortedList {a, b -> a[0] <=> b[0]} |\
            map { l -> l.collect { _chr, bcf -> bcf } } |\
            check_samples

    // a (channel of a) plain text file with one line per chromosome, each line containing the number of SNPs in that vcf, in the same order as the input vcf_file_listx
    c_nsnps_file = c_vcf |\
        count_snps |\
        toSortedList {a, b -> a[0] <=> b[0]} |\
        flatMap |\
        map{ _chr, nsnps -> nsnps } |\
        // collectFile(name: 'nsnps.txt', sort: 'index', storeDir: params.resultsdir + '/intermediate')
        collectFile(name: 'nsnps.txt', sort: 'index')

    // carry out sparse painting (per chromosome)
    c_chunks = c_vcf |\
        paint |\
        toSortedList {a, b -> a[0] <=> b[0]} |\
        flatMap

    // channel of the outputs of sparse painting (a chunklengths file), in the same order as the input vcf_file_list
    c_chunkfiles = c_chunks |\
        map{ _chr, chunks -> chunks }

    // a (channel of a) plain text file with one line per chromosome, each line containing the filename of a chunklengths file, in the same order as the input vcf_file_list
    c_chunkfilenames_file = c_chunks |\
        map{ _chr, chunks -> chunks.name } |\
        // collectFile(name: 'chunklength.files.txt', sort: 'index', newLine: true, storeDir: params.resultsdir + '/intermediate')
        collectFile(name: 'chunklength.files.txt', sort: 'index', newLine: true)

    // to the outputs of pbwt, we add:
    //      1) a list of those files's names in a plain text file
    //      2) a list of the number of SNPs in each chromosome in a plain text file
    //      3) a list of the genetic map lengths for each chromosome in a plain text file
    // then we combine the chunklengths files into a single matrix using combine (called by generate_chunk_length_matrix), 
    // add the sample IDs (which give the dimensions of the sparse matrix, and label the output),
    // and finally compute the haplotype components.
    c_chunkfiles |\
        collect |\
        map { 
            chunks -> [chunks]
        } |\
        combine(c_chunkfilenames_file) |\
        combine(c_nsnps_file) |\
        map { 
            t -> [
                t[0],
                t[1],
                t[2],
                file(params.map_file, checkIfExists: true),
            ]
        } |\
        generate_chunk_length_matrix |\
        combine(c_samples_file) |\
        map { matrix, sums, samples -> [matrix, sums, samples, params.nHCs, params.outfile, params.diagnostics, params.reweight.toString().toBoolean(),
                                      params.remove_pairs as boolean,
                                      params.remove_pairs ? file(params.remove_pairs, checkIfExists: true) : file("${projectDir}/assets/NO_FILE"),
                                      params.fit_samples as boolean,
                                      params.fit_samples ? file(params.fit_samples, checkIfExists: true) : file("${projectDir}/assets/NO_FILE")] } |\
        compute_HCs
}

process count_snps {
    // publishDir "${params.resultsdir}/intermediate", mode: 'copy'

    input:
    tuple val(chr), path(bcf)

    output:
    tuple val(chr), path("*_nsnps.txt")

    script:
    """
    bcftools view -H ${bcf} | wc -l > chr${chr}_nsnps.txt
    """
}

process subset_samples {
    // publishDir "${params.resultsdir}/intermediate", mode: 'copy'

    input:
    tuple val(chr), path(bcf), path(keep)

    output:
    tuple val(chr), path("chr${chr}.subset.bcf")

    script:
    // bcftools fails if a sample in the keep file is not in the input file, and outputs the samples in
    // the order of the keep file. Sites that become monomorphic are kept
    """
    bcftools view -S ${keep} -Ob -o chr${chr}.subset.bcf ${bcf}
    """
}

process check_samples {
    publishDir "${params.resultsdir}/matrix", mode: 'copy', enabled: params.publish_matrix.toString().toBoolean()

    input:
    path(bcfs)

    output:
    path("samples.txt")

    script:
    """
    first=""
    for FILE in ${bcfs}
    do
        bcftools query -l "\${FILE}" > "\${FILE}.samples"
        if [ -z "\${first}" ]; then
            first="\${FILE}"
        elif ! cmp -s "\${first}.samples" "\${FILE}.samples"; then
            echo "Error: \${FILE} does not have the same samples, in the same order, as \${first}." >&2
            exit 1
        fi
    done
    if [ ! -s "\${first}.samples" ]; then
        echo "Error: no samples found in \${first}." >&2
        exit 1
    fi
    dups=\$(sort "\${first}.samples" | uniq -d | head -5)
    if [ -n "\${dups}" ]; then
        echo "Error: duplicate sample IDs in the input files, e.g.: \${dups}" >&2
        exit 1
    fi
    cp "\${first}.samples" samples.txt
    """
}

process paint {
    // publishDir "${params.resultsdir}/intermediate", mode: 'copy'

    input:
    tuple val(chr), path(bcf)

    output:
    tuple val(chr), path("*chunklengths*")

    script:
    """
    pbwt -readVcfGT ${bcf} -paintSparse chr${chr} ${params.paint_nchunks} ${params.paint_ploidy} ${params.paint_cutoff}
    """
}

process generate_chunk_length_matrix {
    publishDir "${params.resultsdir}/matrix", mode: 'copy', enabled: params.publish_matrix.toString().toBoolean()

    input:
    tuple path(chunk_files), path(chunk_filenames_list), path(snp_counts_list), path(map_lengths_list)

    output:
    tuple path("combined.chunklengths.txt.gz"), path("combined.chunklengths.sums")

    script:
    """
    combine -c ${chunk_filenames_list} -s ${snp_counts_list} -m ${map_lengths_list} -r=${params.restrict} -o combined.chunklengths
    """
}

process compute_HCs {
    publishDir params.resultsdir, mode: 'copy'

    input:
    tuple path(chunks), path(sums), path(samples), val(nHCs), val(outfile), val(diagnostics), val(reweight), val(use_pairs), path(remove_pairs, stageAs: 'pairs/*'), val(use_fit), path(fit_samples, stageAs: 'fit/*')

    output:
    tuple path(outfile), path(diagnostics)

    script:
    """
    calculate_haplotype_components \\
        --matrix ${chunks} \\
        --samples ${samples} \\
        --sums ${sums} \\
        --nhcs ${nHCs} \\
        --out ${outfile} \\
        --diagnostics ${diagnostics} \\
        ${reweight ? '--reweight' : ''} \\
        ${use_pairs ? "--remove_pairs ${remove_pairs}" : ''} \\
        ${use_fit ? "--fit_samples ${fit_samples}" : ''}
    """
}
