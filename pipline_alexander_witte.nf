nextflow.enable.dsl = 2

params.accession = "null"
params.outdir = "SRA_data"


process sra_fetch{
    storeDir "${params.outdir}"
   
    container "https://depot.galaxyproject.org/singularity/sra-tools:2.11.0--pl5262h314213e_0"
   
    input: 
    val nr_acc

    output:
    path "sra/${nr_acc}.sra", emit: sradata

    script:
    """
    prefetch $nr_acc
    """
}

process fastq_conv{
    storeDir "${params.outdir}/sra/${nr_acc}_fastq/raw_fastq"
    container "https://depot.galaxyproject.org/singularity/sra-tools:2.11.0--pl5262h314213e_0"

    input:
    val nr_acc
    path sradata

    output:
    path "*.fastq", emit: fastqdata

    script:
    """
    fastq-dump --split-3 ${sradata}
    """
}

process fastq_report {
    
    publishDir "${params.outdir}/sra/${params.accession}_fastq/", mode: 'copy', overwrite: true

    container "https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0"
  
    input:
        path fastqdata
    output:
        path "fastqc_results/*.html", emit: results
        path "fastqc_results/*.zip", emit: fastqczip
    
    script:
        """
        mkdir fastqc_results
        fastqc $fastqdata --outdir fastqc_results
        """
}

process fastp {
  publishDir "${params.outdir}/sra/${nr_acc}_fastq/", mode: 'copy', overwrite: true
  
  container "https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0"

  input:
    path fastqdata
    val  nr_acc
  output:
    path "fastp_fastq/*.fastq", emit: fastqdata
    path "fastp_report", emit: fastpreport

  script:
      if(fastqdata instanceof List) {
      """
      mkdir fastp_fastq
      mkdir fastp_report
      fastp -i ${fastqdata[0]} -I ${fastqdata[1]} -o fastp_fastq/${fastqdata[0].getSimpleName()}_fastp.fastq -O fastp_fastq/${fastqdata[1].getSimpleName()}_fastp.fastq -h fastp_report/fastp.html -j fastp_report/fastp.json
      """
    } else {
      """
      mkdir fastp_fastq
      mkdir fastp_report
      fastp -i ${fastqdata} -o fastp_fastq/${fastqdata.getSimpleName()}_fastp.fastq -h fastp_report/fastp.html -j fastp_report/fastp.json
      """
    }
}

workflow sra_fetch_conv{
    take:
        acc_nr "($params.accession)"
        output "($params.outdir)"
    
    main:
        sra_get = sra_fetch(acc_nr)
        fastq_result = fastq_conv(acc_nr,sra_get).fastqdata
        trimmed_fastq = fastp(fastq_result, acc_nr).fastqdata 
        trimmed_fastq_flat = trimmed_fastq.flatten()
        report = fastq_report(trimmed_fastq.collect()).results
        
    emit:
        fastq_out=trimmed_fastq_flat
        fastqc_report = report
    }

workflow{
    sra_result=sra_fetch_conv(params.accession,params.outdir)
    sra_result.fastq_out.view()
    sra_result.fastqc_report.view()
}