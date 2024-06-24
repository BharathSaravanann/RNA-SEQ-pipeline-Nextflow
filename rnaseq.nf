params.fastq="/home/bharath/pipelines/rna_pipeline/fastq/rna.fastq"
params.index="/home/bharath/pipelines/rna_pipeline/grch38/genome"
params.annotation="/home/bharath/pipelines/rna_pipeline/Homo_sapiens.GRCh38.106.gtf"
params.trimjar="/home/bharath/pipelines/Trimmomatic-0.39/trimmomatic-0.39.jar"
params.output="/home/bharath/pipelines/rna_pipeline/output"

workflow {

    process QC {
        publishDir("${params.output}", mode: 'copy')
        
        input:
        path fastq
        
        output:
        path "*"
        
        script:
        """
        fastqc ${fastq} 
        """
    }

    process trimmomatic {
        publishDir("${params.output}", mode: 'copy')
        
        input:
        val jar
        path fastq
        
        output:
        path "${fastq.baseName}.demo_trimmed.fastq"
        
        script:
        """
        java -jar ${params.trimjar} SE -threads 4 ${fastq} ${fastq.baseName}.demo_trimmed.fastq TRAILING:10 -phred33
        """
    }

    process hisat {
        publishDir("${params.output}", mode: 'copy')
        
        input:
        path index
        path fastq
        
        output:
        path "${fastq.baseName}.demo_trimmed.bam"
        
        script:
        """
        hisat2 -q --rna-strandness R -x ${params.index} -U ${params.fastq} | samtools sort -o ${fastq.baseName}.demo_trimmed.bam
        """
    }

    process featurecount {
        publishDir("${params.output}", mode: 'copy')
        
        input:
        path annotation
        path bam
        
        output:
        path "rna_featurecounts.txt"
        
        script:
        """
        featureCounts -S 2 -a ${params.annotation} -o rna_featurecounts.txt ${bam}
        """
    }

    index_ch = Channel.of(params.index)
    annotation_ch = Channel.of(params.annotation)
    fastq_ch = Channel.of(params.fastq)
    trim_ch = Channel.of(params.trimjar)
    quality_control = QC(fastq_ch)
    trimadapters = trimmomatic(trim_ch, fastq_ch)
    alignment = hisat(trimadapters, index_ch)
    featurecount(annotation_ch, alignment)

}

