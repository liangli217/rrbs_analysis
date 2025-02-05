task fastqc{
    #inputs from workflow config, or upstream task if preprocessing done
    File fastq1
    File fastq2
    String? fastq_suffix = '.fastq.gz'
    String? fastq1_prefix = basename(fastq1, '.fastq.gz')
    String? fastq2_prefix = basename(fastq2, fastq_suffix)

    String? fastqc_args = ""

    #runtime inputs
    String? docker = "gcr.io/broad-cga-bknisbac-wupo1/bisulfite_tools:0.3"
    Int? mem = "3"
    Int? threads = "1"
    Int? disk_size_buffer = "10"
    Int? disk_scaler = "1"
    Int? disk_size_gb = ceil( (size(fastq1, "G") + size(fastq2, "G")) * disk_scaler) + disk_size_buffer
    Int? num_preempt = "4"

    command {
        /src/monitor_script.sh > monitoring.log &

        mkdir fastqc_results

        fastqc \
        --outdir . \
          --threads ${threads} \
          ${fastqc_args} \
          ${fastq1} ${fastq2}

    }
    runtime {
        docker: "${docker}"
        memory: "${mem}GB"
        cpu: "${threads}"
        disks: "local-disk ${disk_size_gb} HDD"
        preemptible: "${num_preempt}"
    }
    output {
        File fq1_zip="${fastq1_prefix}_fastqc.zip"
        File fq2_zip="${fastq2_prefix}_fastqc.zip"
        File fq1_html="${fastq1_prefix}_fastqc.html"
        File fq2_html="${fastq2_prefix}_fastqc.html"
        File monitoring_log="monitoring.log"
    }
}