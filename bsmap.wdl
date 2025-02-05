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

task trim_fastqs{
    #inputs from workflow config
    File fastq1
    File fastq2
    String sample_id #for log output

    String? trim_args=""
    String? trim_cores_arg = "--cores 2" #--cores 2 actually uses 9 cores, 4 uses 15

    #runtime inputs
    String? docker="gcr.io/broad-cga-bknisbac-wupo1/bismark:0.1"
    Int? mem = "4"
    Int? threads = "9"
    Int? disk_size_buffer = "10"
    Int? disk_scaler = "3" #trim_galore makes 2 copies
    Int? disk_size_gb = ceil( (size(fastq1, "G") + size(fastq2, "G")) * disk_scaler) + disk_size_buffer
    Int? num_preempt = "4"

    command {
        /src/monitor_script.sh > monitoring.log &

        #cores in trim_galore is non-intuitive: --cores 2 uses 9 cores; --cores 4 would mean using 15 cores:
        #4 (read) + 4 (write) + 4 (Cutadapt) + 2 (extra Cutadapt) + 1 (Trim Galore) = 15
        #2 (read) + 2 (write) + 2 (Cutadapt) + 2 (extra Cutadapt) + 1 (Trim Galore) = 9
        #See: https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md
        trim_galore --paired --gzip ${trim_cores_arg} ${trim_args} ${fastq1} ${fastq2}
        mv *_val_1.*.gz ${sample_id}_1.trm.fastq.gz
        mv *_val_2.*.gz ${sample_id}_2.trm.fastq.gz

        cat *_trimming_report.txt > "${sample_id}.fastq.trimming.log"
    }

    runtime {
    	    docker: "${docker}"
    	    memory: "${mem}GB"
            cpu: "${threads}"
    	    disks: "local-disk ${disk_size_gb} HDD"
    	    preemptible: "${num_preempt}"
    }
    output {
	    File fastq1_trimmed = "${sample_id}_1.trm.fastq.gz"
	    File fastq2_trimmed = "${sample_id}_2.trm.fastq.gz"
        File trimming_log = "${sample_id}.fastq.trimming.log"
        File monitoring_log="monitoring.log"
    }
}

task trim_diversity_adapters {
  #inputs from workflow config, or upstream task if preprocessing done
  File fastq1
  File fastq2
  String? fastq_suffix = '.fastq'
  String? zip_suffix = '.gz'
  String? full_suffix = fastq_suffix + zip_suffix
  String? fastq1_prefix = basename(fastq1, zip_suffix)
  String? fastq2_prefix = basename(fastq2, zip_suffix)

  String? trim_suff = ".trmDivAdp"

  String? trim_div_args = ""

  #runtime inputs
  String? docker = "gcr.io/broad-cga-bknisbac-wupo1/divtrim:1.0"
  Int? mem = "3"
  Int? threads = "1"
  Int? disk_size_buffer = "10"
  Int? disk_scaler = "2"
  Int? disk_size_gb = ceil( (size(fastq1, "G") + size(fastq2, "G")) * disk_scaler) + disk_size_buffer
  Int? num_preempt = "4"

  command {
      /src/monitor_script.sh > monitoring.log &

      mkdir workdir
      new_fastq1=workdir/$(basename ${fastq1})
      new_fastq2=workdir/$(basename ${fastq2})
      mv ${fastq1} $new_fastq1
      mv ${fastq2} $new_fastq2

      python2 /src/trimRRBSdiversityAdaptCustomers.py -1 $new_fastq1 -2 $new_fastq2 ${trim_div_args}

      mv workdir/"${fastq1_prefix}_trimmed.fq${zip_suffix}" "${fastq1_prefix}${trim_suff}${full_suffix}"
      mv workdir/"${fastq2_prefix}_trimmed.fq${zip_suffix}" "${fastq2_prefix}${trim_suff}${full_suffix}"
  }
  runtime {
      docker: "${docker}"
      memory: "${mem}GB"
      cpu: "${threads}"
      disks: "local-disk ${disk_size_gb} HDD"
      preemptible: "${num_preempt}"
  }
  output {
      File fastq1_divtrim="${fastq1_prefix}${trim_suff}${full_suffix}"
      File fastq2_divtrim="${fastq2_prefix}${trim_suff}${full_suffix}"
      File monitoring_log="monitoring.log"
  }
}


task bamstat {
    File bam_file
    File bam_index
    String prefix = basename(bam_file, ".bam")

    #runtime inputs
    String? docker="gcr.io/broad-cga-bknisbac-wupo1/bisulfite_tools:0.3"
    Int? mem = "3"
    Int? threads = "4"
    Int? disk_size_buffer = "10"
    Int? disk_scaler = "1"
    Int? disk_size_gb = ceil( (size(bam_file, "G") + size(bam_index, "G")) * disk_scaler) + disk_size_buffer
    Int? num_preempt = "4"

    command {
            /src/monitor_script.sh > monitoring.log &
            /src/bamStat.pl --bam ${bam_file} --threads ${threads} --verbose > ${prefix}.stats.txt

        }

    runtime {
        docker: "${docker}"
        memory: "${mem}GB"
        cpu: "${threads}"
        disks: "local-disk ${disk_size_gb} HDD"
        preemptible: "${num_preempt}"
    }
    output {
        File bam_stats = "${prefix}.stats.txt"
        File monitoring_log="monitoring.log"
    }
}

task bsmap{
    File fastq1
    File fastq2
    String sample_id
    File reference_fa
    Int? seed_size="12" #default=12(RRBS mode), 16(WGBS mode). min=8, max=16.
    Int? max_insert_size="1000" # max insert size for PE mapping (-x)
    # -q is quality threshold.
    # -w<int>   maximum number of equal best hits to count, <=1000
    # -S seed for rng.  0 for system clock (not reproducible) otherwise produces reproducible results.
    # -u   report unmapped reads, default=off
    # -R          print corresponding reference sequences in SAM output, default=off
    # -r  [0,1]   how to report repeat hits, 0=none(unique hit/pair only); 1=random one, default=1.
    # -D C-CGG    Add if want to align specifically to mspI cut sites (consider for RRBS, not WGBS, and beware if trimming may discard the mspI cut site seq)

    # Note: the UMI pipeline using nudup requires unique alignments (then use: -w 1; CLL-map PE default was: "-q 20 -w 100 -S 1 -u -R")
    String? bsmap_args="-q 20 -w 100 -S 1 -u -R"

    #runtime inputs
    String? docker="gcr.io/broad-cga-bknisbac-wupo1/bisulfite_tools:0.3"
    #mem=10, threads=12 is sufficient for "samtools sort" default mem/thread usage [768 MiB]
    Int? mem = "10"
    Int? threads = "24"
    Int? disk_size_buffer = "10"
    Int? disk_scaler = "2"
    Int? disk_size_gb = ceil( (size(fastq1, "G") + size(fastq2, "G")) * disk_scaler) + disk_size_buffer
    Int? num_preempt = "4"

    command {
            /src/monitor_script.sh > monitoring.log &
            # -s = seed size, default=12(RRBS mode), 16(WGBS mode). min=8, max=16
            # -u = report unmapped reads, default=off
            # -R = print corresponding reference sequences in SAM output, default=off
            # -q =  quality threshold in trimming, 0-40, default=0 (no trim)
            # -w =  maximum number of equal best hits to count, <=1000
            # -S = 1 seed for random number generation used in selecting multiple hits, keep nonzero for reproucibility
            bsmap \
              -v 0.1 -s ${seed_size} ${bsmap_args} \
              -x ${max_insert_size} \
              -p ${threads} \
              -d ${reference_fa} \
              -a ${fastq1} \
              -b ${fastq2} \
              -o ${sample_id}.bsmap.bam 2> ${sample_id}.bsmap.stderr.log

            /src/parse_bsmap_report.py -s ${sample_id} -f ${sample_id}.bsmap.stderr.log -o ${sample_id}.bsmap.stats.tsv

            echo START stderr from bsmap:
            cat ${sample_id}.bsmap.stderr.log
            echo END stderr from bsmap.
 	    }
    runtime {
        docker: "${docker}"
        memory: "${mem}GB"
        cpu: "${threads}"
        disks: "local-disk ${disk_size_gb} HDD"
        preemptible: "${num_preempt}"
    }
    output {
        File bam = "${sample_id}.bsmap.bam"
        File stats_tsv = "${sample_id}.bsmap.stats.tsv"
        File monitoring_log="monitoring.log"
    }
}

task markdup_umi {
  File bam_file
  File umi_index_fastq
  String? dedup_umi_args = "--paired-end --start 6 --length 6 -T temp"

  String? file_pref = basename(bam_file, '.bam')
  String? bam_outsuff = ".sorted.md_umi.bam"

  #runtime inputs
  String? docker = "gcr.io/broad-cga-bknisbac-wupo1/nudup:1.0"
  Int? mem = "8"
  Int? threads = "4"
  Int? disk_size_buffer = "10"
  Int? disk_scaler = "3"
  Int? disk_size_gb = ceil( size(bam_file, "G") * disk_scaler) + disk_size_buffer
  Int? num_preempt = "4"

  command {
      /src/monitor_script.sh > monitoring.log &

      mkdir temp

      python2 /src/nudup.py ${dedup_umi_args} \
        -f ${umi_index_fastq} \
        --out ${file_pref} \
        ${bam_file}

        mv ${file_pref}.sorted.markdup.bam ${file_pref}${bam_outsuff}
        rm ${file_pref}.sorted.dedup.bam #I didn't need dedup file

        samtools index -@ ${threads} ${file_pref}${bam_outsuff}

  }
  runtime {
      docker: "${docker}"
      memory: "${mem}GB"
      cpu: "${threads}"
      disks: "local-disk ${disk_size_gb} HDD"
      preemptible: "${num_preempt}"
  }
  output {
      File markdup_bam = "${file_pref}${bam_outsuff}"
      File markdup_bai = "${file_pref}${bam_outsuff}.bai"
      File markdup_log = "${file_pref}_dup_log.txt"
      File monitoring_log="monitoring.log"
  }
}

task sort_bam {
    File bam_file
    String prefix = basename(bam_file, ".bam")

    #runtime inputs
    String? docker="gcr.io/broad-cga-bknisbac-wupo1/bisulfite_tools:0.3"
    #mem=10, threads=12 is sufficient for "samtools sort" default mem/thread usage [768 MiB]
    Int? mem = "10"
    Int? threads = "12"
    Int? disk_size_buffer = "10"
    Int? disk_scaler = "2"
    Int? disk_size_gb = ceil( size(bam_file, "G") * disk_scaler) + disk_size_buffer
    Int? num_preempt = "4"

    #sort args (here because some derived)
    String? sort_args = ""
    Int? mem_mb_scaling_factor = "900"
    Int mem_mb_per_sort_thread = floor(mem_mb_scaling_factor * mem / threads) #would have done 1000 but don't want to risk choking machine

    command {
            /src/monitor_script.sh > monitoring.log &

            samtools sort ${sort_args} --threads ${threads} -m ${mem_mb_per_sort_thread}M --output-fmt BAM -o ${prefix}.srt.bam ${bam_file}

            samtools index -@ ${threads} ${prefix}.srt.bam

 	    }
    runtime {
        docker: "${docker}"
        memory: "${mem}GB"
        cpu: "${threads}"
        disks: "local-disk ${disk_size_gb} HDD"
        preemptible: "${num_preempt}"
    }
    output {
        File bam = "${prefix}.srt.bam"
        File bam_index = "${prefix}.srt.bam.bai"
        File monitoring_log="monitoring.log"
    }
}

task markduplicates {
    File bam_file
    File bam_index
    String sample_id

    String? bam_prefix=basename(bam_file, ".bam")

    String? remove_dups="false" #WARNING: "true" LEADS TO LOSS OF READS IN FINAL BAM
    Int? max_records_in_ram="500000" #As in GTEx; Helene: "10000000"
    Int? max_open_files="8000" #GATK default
    Int? java_mem="3" # max heap size #As in GTEx; Helene: "16G"
    String? read_name_regex="null" #skips optical duplicate finding (this is setting for HiSeq 4000 and NovaSeq 6000 that shouldnt have optical duplicates)
    String? other_args=""

    #runtime inputs
    String? docker = "gcr.io/broad-cga-bknisbac-wupo1/bisulfite_tools:0.3"
    Int? mem = "3"
    Int? threads = "1"
    Int? disk_size_buffer = "20"
    Int? disk_scaler = "3" # 3 for potential intermediate files
    Int? disk_size_gb = ceil( size(bam_file, "G") + size(bam_index, "G") * disk_scaler) + disk_size_buffer
    Int? num_preempt = "4"

    command {
        /src/monitor_script.sh > monitoring.log &

        bam_dedup=${bam_prefix}.md.bam

        ## We skip optical duplicate finding by setting READ_NAME_REGEX=null.
        ## With HiSeq 4000 and NovaSeq 6000 there shouldn't be any "optical duplicates"
        ## due to their patterned flowcells. We will run only few WGBS libs on other machines.
        ## sklages, 2018-11-15
        temp_dir="tempdir"
        mkdir $temp_dir

        gatk \
          MarkDuplicates \
          --java-options "-XX:ParallelGCThreads=${threads} -Xmx${java_mem}g -Djava.io.tmpdir=$temp_dir" \
          --ASSUME_SORT_ORDER=coordinate \
          --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=${max_open_files} \
          --REMOVE_DUPLICATES=${remove_dups} \
          --MAX_RECORDS_IN_RAM=${max_records_in_ram} \
          --INPUT=${bam_file} \
          --OUTPUT=${bam_prefix}.md.bam \
          --METRICS_FILE=${bam_prefix}.md-metrics.txt \
          --TMP_DIR=$temp_dir \
          --VALIDATION_STRINGENCY=LENIENT \
          --READ_NAME_REGEX=${read_name_regex} \
          ${other_args}

        samtools index -@ ${threads} ${bam_prefix}.md.bam

    }

    runtime {
        docker: "${docker}"
        memory: "${mem}GB"
        cpu: "${threads}"
        disks: "local-disk ${disk_size_gb} HDD"
        preemptible: "${num_preempt}"
    }
    output {
        File bam_md="${bam_prefix}.md.bam"
        File bam_md_index="${bam_prefix}.md.bam.bai"
        File bam_md_metrics="${bam_prefix}.md-metrics.txt"
        File monitoring_log="monitoring.log"
    }
}

task mcall {
    File bam_file #either bsmap_bam or bam_md
    File bam_index
    String sample_id
    File reference_fa
    File reference_sizes

    #No trimming: --trimRRBSEndRepairSeq 0
    #For RRBS/WGBS (automatically detected): --trimRRBSEndRepairSeq 2
    #-F 256 for primary alignments only
    String? mcall_args="-F 256"

    String? bam_prefix=basename(bam_file, ".bam")

    #runtime inputs
    String? docker = "gcr.io/broad-cga-bknisbac-wupo1/bisulfite_tools:0.3"
    Int? mem = "8"
    Int? threads = "16"
    Int? disk_size_buffer = "10"
    Int? disk_scaler = "2"
    Int? disk_size_gb = ceil( size(bam_file, "G") + size(bam_index, "G") * disk_scaler) + disk_size_buffer
    Int? num_preempt = "4"

    command {
        /src/monitor_script.sh > monitoring.log &

        mcall \
          --threads ${threads} \
          --reference ${reference_fa} \
          --sampleName ${sample_id} \
          --mappedFiles ${bam_file} \
          ${mcall_args}

        mv ${bam_file}.G.bed ${sample_id}.CpG.bed

        #Stats are not unique to CpG (would be the same for e.g. CpA analysis), hence naming CpX.stats
        mv ${bam_file}_stat.txt ${sample_id}.mcall.stats.txt

        ## convert BED to BigBED using bedToBigBed
        echo "Converting BED (${sample_id}.CpG.bed) to BigBed"
        tmp_file=$(/src/mcall_bed_reformat.pl ${sample_id}.CpG.bed)
        bedToBigBed $tmp_file ${reference_sizes} ${sample_id}.CpG.bb
        rm -f $tmp_file

        gzip ${sample_id}.CpG.bed

    }
    runtime {
        docker: "${docker}"
        memory: "${mem}GB"
        cpu: "${threads}"
        disks: "local-disk ${disk_size_gb} HDD"
        preemptible: "${num_preempt}"
    }
    output {
        File mcall_cpg_bed_gz = "${sample_id}.CpG.bed.gz"
        File mcall_cpg_bigbed = "${sample_id}.CpG.bb"
        File mcall_stats = "${sample_id}.mcall.stats.txt"
        File monitoring_log="monitoring.log"
    }
}

#https://multiqc.info/docs/
### Should have access to:
#fastqc for raw fastqs
#fastqc for trimmed fastqs
#cutadapt
task multiqc {
    String sample_id
    File fastq1_fastqc
    File fastq2_fastqc
    File? fastq1_trimmed_fastqc
    File? fastq2_trimmed_fastqc
    File? trimming_log
    File? bam_markdups_report
    String? multiqc_args = ""

    #runtime inputs
    String? docker = "gcr.io/broad-cga-bknisbac-wupo1/bisulfite_tools:0.3"
    Int? mem = "3"
    Int? threads = "1"
    Int? disk_size_buffer = "10"
    Int? disk_scaler = "1"
    Int? disk_size_gb = ceil( (size(fastq1_fastqc, "G") + size(fastq2_fastqc, "G") + size(fastq1_trimmed_fastqc, "G") + size(fastq2_trimmed_fastqc, "G")) * disk_scaler) + disk_size_buffer
    Int? num_preempt = "4"

    command {
        /src/monitor_script.sh > monitoring.log &

        multiqc_dir="${sample_id}.multiqc"
        mkdir $multiqc_dir

        mv ${fastq1_fastqc} $multiqc_dir/
        mv ${fastq2_fastqc} $multiqc_dir/

        if [ -n "${trimming_log}" ]; then
            mv ${trimming_log} $multiqc_dir/
        fi

        if [ -n "${fastq1_trimmed_fastqc}" ]; then
            mv ${fastq1_trimmed_fastqc} $multiqc_dir/
        fi

        if [ -n "${fastq2_trimmed_fastqc}" ]; then
            mv ${fastq2_trimmed_fastqc} $multiqc_dir/
        fi

        if [ -n "${bam_markdups_report}" ]; then
            mv ${bam_markdups_report} $multiqc_dir/
        fi

        multiqc $multiqc_dir --filename ${sample_id}.multiqc_report ${multiqc_args} 2>multiqc.stderr

        echo creating tar with:
        tar -czvf $multiqc_dir.tar.gz $multiqc_dir

    }

    runtime {
        docker: "${docker}"
        memory: "${mem}GB"
        cpu: "${threads}"
        disks: "local-disk ${disk_size_gb} HDD"
        preemptible: "${num_preempt}"
    }
    output {
        File multiqc_report_html="${sample_id}.multiqc_report.html"
        File multiqc_tar_gz="${sample_id}.multiqc.tar.gz"
        File monitoring_log="monitoring.log"
    }
}

workflow bsmap_to_mcall_PE {
    String sample_id
    File fastq1
    File fastq2
    File reference_fa #for bsmap, mcall
    File reference_sizes #for bsmap, mcall

    ## Control variables
    Boolean? run_trim = true
    Boolean? run_trim_div_adapts = true #depends on run_trim
    Boolean? run_sort_bam = true
    Boolean? run_markdup_umi = false
    Boolean? run_markduplicates = true

    ### workflow-wide optional runtime variables
    String? docker="gcr.io/broad-cga-bknisbac-wupo1/bisulfite_tools:0.3"
    String? docker_trim="gcr.io/broad-cga-bknisbac-wupo1/bismark:0.1"
    String? docker_trim_div_adapts="gcr.io/broad-cga-bknisbac-wupo1/divtrim:1.0"
    String? docker_markdup_umi="gcr.io/broad-cga-bknisbac-wupo1/nudup:1.0"

    Int? num_preempt="4"

    #### Per tasks ####
    ## for fastqc
    String? fastq_suffix
    String? fastqc_args_raw
    String? fastqc_args_trimmed

    ## for trim_fastqs
    String? trim_fastqs_fastq_suffix
    String? trim_fastqs_adaptors
    Int? trim_fastqs_min_length
    Int? trim_fastqs_min_qual
    Int? trim_fastqs_swift_trim_off
    String? trim_fastqs_trimming_type

    ## for trim_diversity_adapters
    String? trim_div_adapts_fastq_suffix
    String? trim_div_adapts_zip_suffix
    String? trim_div_adapts_suff
    String? trim_div_adapts_args

    ## for bsmap
    Int? bsmap_seed_size
    Int? bsmap_max_insert_size
    String? bsmap_args

    ## for markdup_umi
    File? umi_index_fastq
    String? markdup_umi_args
    String? markdup_umi_bam_outsuff

    ## sort args
    Int? sort_bam_args
    Int? sort_bam_mem_mb_scaling_factor

    ## for markduplicates
    String? markdups_remove_dups #WARNING: "true" LEADS TO LOSS OF READS IN FINAL BAM
    Int? markdups_max_records_in_ram
    Int? markdups_max_open_files
    Int? markdups_java_mem
    String? markdups_read_name_regex #skips optical duplicate finding (this is setting for HiSeq 4000 and NovaSeq 6000 that shouldnt have optical duplicates)
    String? markdups_other_args

    ## for mcall
    String? mcall_args

    ## for multiqc
    String? multiqc_args

    ### Per-task optional runtime variables
    Int? trim_fastqs_mem
    Int? trim_fastqs_threads
    Int? trim_fastqs_disk_size_buffer
    Int? trim_fastqs_num_preempt

    Int? trim_div_adapts_mem
    Int? trim_div_adapts_threads
    Int? trim_div_adapts_disk_size_buffer
    Int? trim_div_adapts_disk_scaler
    Int? trim_div_adapts_num_preempt

    Int? fastqc_mem
    Int? fastqc_threads
    Int? fastqc_disk_size_buffer
    Int? fastqc_num_preempt

    Int? bsmap_mem
    Int? bsmap_threads
    Int? bsmap_disk_size_buffer
    Int? bsmap_num_preempt

    Int? markdup_umi_mem
    Int? markdup_umi_threads
    Int? markdup_umi_disk_size_buffer
    Int? markdup_umi_disk_scaler
    Int? markdup_umi_num_preempt

    Int? sort_bam_mem
    Int? sort_bam_threads
    Int? sort_bam_disk_size_buffer
    Int? sort_bam_num_preempt

    Int? bamstat_mem
    Int? bamstat_threads
    Int? bamstat_disk_size_buffer
    Int? bamstat_num_preempt

    Int? markdups_mem
    Int? markdups_threads
    Int? markdups_disk_size_buffer
    Int? markdups_num_preempt

    Int? mcall_mem
    Int? mcall_threads
    Int? mcall_disk_size_buffer
    Int? mcall_num_preempt

    Int? multiqc_mem
    Int? multiqc_threads
    Int? multiqc_disk_size_buffer
    Int? multiqc_num_preempt


  if(run_trim==true){
      call trim_fastqs {
          input:
              fastq1=fastq1,
              fastq2=fastq2,
              sample_id=sample_id,
              docker = select_first([docker_trim, docker]),
              mem = trim_fastqs_mem,
              threads = trim_fastqs_threads,
              disk_size_buffer = trim_fastqs_disk_size_buffer,
              num_preempt = select_first([trim_fastqs_num_preempt, num_preempt])
      }

      if(run_trim_div_adapts==true){
          call trim_diversity_adapters {
              input:
                  fastq1 = trim_fastqs.fastq1_trimmed,
                  fastq2 = trim_fastqs.fastq2_trimmed,
                  fastq_suffix = trim_div_adapts_fastq_suffix,
                  zip_suffix = trim_div_adapts_zip_suffix,
                  trim_suff = trim_div_adapts_suff,
                  trim_div_args = trim_div_adapts_args,
                  docker = select_first([docker_trim_div_adapts, docker]),
                  mem = trim_div_adapts_mem,
                  threads = trim_div_adapts_threads,
                  disk_size_buffer = trim_div_adapts_disk_size_buffer,
                  disk_scaler = trim_div_adapts_disk_scaler,
                  num_preempt = select_first([trim_div_adapts_num_preempt, num_preempt])
          }
      }
  }

  File? fastq1_trimmed = select_first([trim_diversity_adapters.fastq1_divtrim, trim_fastqs.fastq1_trimmed])
  File? fastq2_trimmed = select_first([trim_diversity_adapters.fastq2_divtrim, trim_fastqs.fastq2_trimmed])

    call fastqc as fastqc_raw {
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            fastq_suffix = fastq_suffix,
            fastqc_args = fastqc_args_raw,
            docker = docker,
            mem = fastqc_mem,
            threads = fastqc_threads,
            disk_size_buffer = fastqc_disk_size_buffer,
            num_preempt = select_first([fastqc_num_preempt, num_preempt])
    }

    call fastqc as fastqc_trimmed {
        input:
            fastq1 = fastq1_trimmed,
            fastq2 = fastq2_trimmed,
            docker = docker,
            mem = fastqc_mem,
            threads = fastqc_threads,
            disk_size_buffer = fastqc_disk_size_buffer,
            num_preempt = select_first([fastqc_num_preempt, num_preempt])
    }


    call bsmap{
        input:
            fastq1=select_first([fastq1_trimmed, fastq1]),
            fastq2=select_first([fastq2_trimmed, fastq2]),
            sample_id=sample_id,
            reference_fa=reference_fa,
            seed_size=bsmap_seed_size,
            max_insert_size=bsmap_max_insert_size,
            bsmap_args=bsmap_args,
            docker = docker,
            mem = bsmap_mem,
            threads = bsmap_threads,
            disk_size_buffer = bsmap_disk_size_buffer,
            num_preempt = select_first([bsmap_num_preempt, num_preempt])
    }

if(run_markdup_umi==true){ # markdup_umi sorts and marks dups
    call markdup_umi{
      input:
          bam_file=bsmap.bam,
          umi_index_fastq=umi_index_fastq,
          dedup_umi_args=markdup_umi_args,
          bam_outsuff=markdup_umi_bam_outsuff,

          docker=docker_markdup_umi,
          mem = markdup_umi_mem,
          threads = markdup_umi_threads,
          disk_size_buffer = markdup_umi_disk_size_buffer,
          disk_scaler=markdup_umi_disk_scaler,
          num_preempt = markdup_umi_num_preempt
    }
}

if(run_markdup_umi==false){ # markdup_umi sorts and marks dups
  if(run_sort_bam==true){
      call sort_bam {
          input:
              bam_file=bsmap.bam,
              sort_args=sort_bam_args,
              mem_mb_scaling_factor=sort_bam_mem_mb_scaling_factor,
              mem = sort_bam_mem,
              threads = sort_bam_threads,
              disk_size_buffer = sort_bam_disk_size_buffer,
              num_preempt = sort_bam_num_preempt
      }
  }

  call bamstat as bamstat_bsmap{
      input:
          bam_file = sort_bam.bam,
          bam_index = sort_bam.bam_index,
          docker = docker,
          mem = bamstat_mem,
          threads = bamstat_threads,
          disk_size_buffer = bamstat_disk_size_buffer,
          num_preempt = select_first([bamstat_num_preempt, num_preempt])
  }

  if(run_markduplicates==true){
      call markduplicates{
          input:
              bam_file = sort_bam.bam,
              bam_index = sort_bam.bam_index,
              sample_id = sample_id,
              remove_dups = markdups_remove_dups,
              max_records_in_ram = markdups_max_records_in_ram,
              max_open_files = markdups_max_open_files,
              java_mem = markdups_java_mem,
              read_name_regex = markdups_read_name_regex,
              other_args = markdups_other_args,
              docker = docker,
              mem = markdups_mem,
              threads = markdups_threads,
              disk_size_buffer = markdups_disk_size_buffer,
              num_preempt = select_first([markdups_num_preempt, num_preempt])
      }
      call bamstat as bamstat_md{
          input:
              bam_file = markduplicates.bam_md,
              bam_index = markduplicates.bam_md_index,
              docker = docker,
              mem = bamstat_mem,
              threads = bamstat_threads,
              disk_size_buffer = bamstat_disk_size_buffer,
              num_preempt = select_first([bamstat_num_preempt, num_preempt])
      }
  }
}

File? markdup_bam = select_first([markduplicates.bam_md, markdup_umi.markdup_bam])
File? markdup_bam_index = select_first([markduplicates.bam_md_index, markdup_umi.markdup_bai])

    call mcall {
        input:
            bam_file = select_first([markdup_bam, sort_bam.bam]),
            bam_index = select_first([markdup_bam_index, sort_bam.bam_index]),
            sample_id = sample_id,
            reference_fa=reference_fa,
            reference_sizes=reference_sizes,
            mcall_args = mcall_args,
            docker = docker,
            mem = mcall_mem,
            threads = mcall_threads,
            disk_size_buffer = mcall_disk_size_buffer,
            num_preempt = select_first([mcall_num_preempt, num_preempt])
    }

    call multiqc{
        ### Get .zip & .html for fastqc ; trimmed report; markdup report
        ### Outputs will include a tar.gz + multiqc output
        input:
            sample_id = sample_id,
            fastq1_fastqc = fastqc_raw.fq1_zip,
            fastq2_fastqc = fastqc_raw.fq2_zip,
            fastq1_trimmed_fastqc = fastq1_trimmed,
            fastq2_trimmed_fastqc = fastq2_trimmed,
            trimming_log = trim_fastqs.trimming_log,
            bam_markdups_report = markduplicates.bam_md_metrics,
            multiqc_args = multiqc_args,
            docker = docker,
            mem = multiqc_mem,
            threads = multiqc_threads,
            disk_size_buffer = multiqc_disk_size_buffer,
            num_preempt = select_first([multiqc_num_preempt, num_preempt])
    }


    File? bamstats_bsmap_file = bamstat_bsmap.bam_stats
    File? bamstats_md_file = bamstat_md.bam_stats

    meta {
        author: "Binyamin A. Knisbacher"
        email: "bknisbac@broadinstitute.org"
        author: "Andrew Dunford (wrote bsmap docker)"
        email: "adunford@broadinstitute.org"
    }

    output {
        #bsmap
        File bsmap_bam_final = select_first([markdup_umi.markdup_bam, markduplicates.bam_md, sort_bam.bam])
        File bsmap_bam_final_index = select_first([markdup_umi.markdup_bai, markduplicates.bam_md_index, sort_bam.bam_index])
        Array[File] bsmap_align_stats = select_all([bamstats_bsmap_file, bamstats_md_file])
        File bsmap_stats_tsv = bsmap.stats_tsv

        #multiqc
        File multiqc_report_html = multiqc.multiqc_report_html
        File multiqc_tar_gz = multiqc.multiqc_tar_gz
        #mcall
        File mcall_cpg_bed_gz = mcall.mcall_cpg_bed_gz
        File mcall_cpg_bigbed = mcall.mcall_cpg_bigbed
        File mcall_stats = mcall.mcall_stats
    }
}