nextflow.enable.dsl = 2

params.vcf_file_list = ""
params.map_file = ""
params.nHCs = 150

params.outfile = "HCs.csv"
params.resultsdir = "results"

workflow {

    // a (channel of a) tuple of the chromosome index and the path to the bcf file for that chromosome,
    // parsed from a plain text file with one line per chromosome, each line containing the chromosome index and the path to the bcf file for that chromosome
    // with columns separated by whitespace
    c_vcf = channel.from(
        file(params.vcf_file_list)
            .readLines()
        ) |\
        map { line -> 
                def fields = line.split("\\h+")
                if (fields.size() != 2) {
                    throw new Exception("Each line in the vcf_file_list must contain exactly two fields: chromosome number and path to vcf/bcf file.")
                }
                return [String.format('%02d', Integer.parseInt(fields[0])), file(fields[1])]
            
        }
    
    // a (channel of a) plain text file with one line containing the number of individuals in the dataset
    c_nindividuals_val = c_vcf |\
            map { _i, bcf -> bcf } |\
            collect |\
            count_individuals |\
            map { file -> file.text.trim() }

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
    // add the number of individuals in the dataset in order to make the sparse matrix with the correct dimensions, 
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
                params.map_file,
            ]
        } |\
        generate_chunk_length_matrix |\
        combine(c_nindividuals_val) |\
        map { file, _sums, nind -> [file, nind, params.nHCs, params.outfile] } |\
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

process count_individuals {
    // publishDir "${params.resultsdir}/intermediate", mode: 'copy'

    input:
    path(bcf)

    output:
    path("nind.txt")

    script:
    """
    for FILE in *bcf
    do
        bcftools query -l \${FILE} | wc -l >> ninds.txt
    done
    sort ninds.txt | uniq > nind.txt
    if [ \$(wc -l < nind.txt) -gt 1 ]; then
        echo "Error: More than one unique number of individuals found."
        exit 1
    fi
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
    pbwt -readVcfGT ${bcf} -paintSparse chr${chr} 100 2 500
    """
}

process generate_chunk_length_matrix {
    // publishDir "${params.resultsdir}/intermediate", mode: 'copy'

    input:
    tuple path(chunk_files), path(chunk_filenames_list), path(snp_counts_list), path(map_lengths_list)

    output:
    tuple path("combined.chunklengths.txt.gz"), path("combined.chunklengths.sums")

    script:
    """
    combine -c ${chunk_filenames_list} -s ${snp_counts_list} -m ${map_lengths_list} -r -o combined.chunklengths
    """
}

process compute_HCs {
    publishDir params.resultsdir, mode: 'copy'

    input:
    tuple path(chunks), val(nind), val(nHCs), val(outfile)

    output:
    path("*")

    script:
    """
    calculate_haplotype_components ${chunks} ${nind} ${nHCs} ${outfile}
    """
}
