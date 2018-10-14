configfile: "sampleconfig.yaml"

KMERSIZE      = config["KMERSIZE"]
ABUNDANCE_MIN = config["ABUNDANCE_MIN"]
TMP_DIR       = config["TMP_DIR"]
OUT_DIR       = config["OUT_DIR"]

rule all: 
    input:
        kmer_results   = expand("cvfiles/{x}.kmer_results", x=config["CONTIG_FILES"]),
        bwa_results    = expand("cvfiles/{x}.bwa_align_results", x=config["CONTIG_FILES"]),
        exact_results  = "cvfiles/all.exact_align_results"
    
    output:
        temp("cvfiles/temp0"),
        temp("cvfiles/temp1"),
        temp("cvfiles/temp2"),
        temp("cvfiles/temp3"),
        "cvfiles/results.txt"

    shell:
        """
        echo -e "kmer_recall%\tkmer_precision%\tbwa_align%\texact_align%" > cvfiles/temp0
        cat {input.kmer_results} > cvfiles/temp1
        cat {input.bwa_results}  | sed s/\%//g > cvfiles/temp2
        paste cvfiles/temp1 cvfiles/temp2 {input.exact_results} > cvfiles/temp3
        cat cvfiles/temp0 cvfiles/temp3 > cvfiles/results.txt
        """


rule concat_refs:
    input:
        {config["PRIMARY_REF"]},
        expand("{x}", x=config["SECONDARY_REFS"])

    output:
        "cvfiles/refs.fa"
    
    shell:
        "cat {input} > cvfiles/refs.fa"

rule kmc_ref:
    input:
        "cvfiles/refs.fa"
    output:
        "cvfiles/refs.fa.kmc_pre",
        "cvfiles/refs.fa.kmc_suf"
    threads: 
        4
    shell:
        "kmc -t{threads} -ci1 -k{KMERSIZE} -fm {input} cvfiles/refs.fa {TMP_DIR}"

rule kmc_contigs:
    input:
        "{contigs}"

    output:
        "cvfiles/all.{contigs}.kmc_pre",
        "cvfiles/all.{contigs}.kmc_suf"

    threads: 
        4

    shell:
        "kmc -t{threads} -ci{ABUNDANCE_MIN} -k{KMERSIZE} -fm {input} cvfiles/all.{wildcards.contigs} {TMP_DIR}"

rule intersect_kmers:
    input:
        "cvfiles/all.{contigs}.kmc_pre",
        "cvfiles/all.{contigs}.kmc_suf",
        "cvfiles/refs.fa.kmc_pre",
        "cvfiles/refs.fa.kmc_suf"

    output:
        temp("cvfiles/tp.{contigs}.kmc_pre"),
        temp("cvfiles/tp.{contigs}.kmc_suf"),
        temp("cvfiles/fn.{contigs}.kmc_pre"),
        temp("cvfiles/fn.{contigs}.kmc_suf"),
        temp("cvfiles/fp.{contigs}.kmc_pre"),
        temp("cvfiles/fp.{contigs}.kmc_suf"),
        "cvfiles/{contigs}.kmer_results"

    shell:
        """
        kmc_tools simple cvfiles/refs.fa -ci1 cvfiles/all.{wildcards.contigs} -ci1 intersect cvfiles/tp.{wildcards.contigs} -ci1
        kmc_tools simple cvfiles/refs.fa -ci1 cvfiles/all.{wildcards.contigs} -ci1 kmers_subtract cvfiles/fn.{wildcards.contigs} -ci1
        kmc_tools simple cvfiles/refs.fa -ci1 cvfiles/all.{wildcards.contigs} -ci1 reverse_kmers_subtract cvfiles/fp.{wildcards.contigs} -ci1
        tp=$(src/count_kmers_kmc cvfiles/tp.{wildcards.contigs})
        fn=$(src/count_kmers_kmc cvfiles/fn.{wildcards.contigs})
        fp=$(src/count_kmers_kmc cvfiles/fp.{wildcards.contigs})
        precision=$(bc -l <<< "scale=4; $tp * 100 / ($tp + $fp)")
        recall=$(bc -l <<< "scale=4; $tp * 100 / ($tp + $fn)")
        echo -e "$recall\t$precision" >> cvfiles/{wildcards.contigs}.kmer_results
        """


rule bwa_index_ref:
    input:
        {config["PRIMARY_REF"]}
    output:
        "cvfiles/refindex.bwt",
        "cvfiles/refindex.sa",
        "cvfiles/refindex.pac",
        "cvfiles/refindex.ann",
        "cvfiles/refindex.amb"
    shell:
        "bwa index -p cvfiles/refindex {input}"

rule bwa_align:
    input:
        "cvfiles/refindex.bwt",
        "{contigs}"

    output:
        "cvfiles/{contigs}.bwa_align_results",
        temp("cvfiles/{contigs}.bwa.bam"),
        temp("cvfiles/{contigs}.bwa.bam.bai")

    threads: 8

    shell:
        """
        bwa mem -t {threads} cvfiles/refindex {wildcards.contigs} | samtools sort > cvfiles/{wildcards.contigs}.bwa.bam
        samtools index cvfiles/{wildcards.contigs}.bwa.bam
        samtools flagstat cvfiles/{wildcards.contigs}.bwa.bam | grep mapped | cut -f 5 -d " "| cut -f 2 -d "(" | head -1 >> cvfiles/{wildcards.contigs}.bwa_align_results
        """

rule exact_align:
    input:
        refs = "cvfiles/refs.fa",
        contig_files = expand("{x}", x=config["CONTIG_FILES"])

    output:
        "cvfiles/all.exact_align_results",
        "cvfiles/all.exact_align_results.exact"

    shell:
        "src/exactalign cvfiles/refs.fa  cvfiles/all.exact_align_results {input.contig_files}"

        
