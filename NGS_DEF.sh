#!/bin/bash

# 设置选项：命令执行失败时立即退出，未定义变量报错
set -eu

# 定义变量
export JAVA_HOME=/home/stereonote/software/jdk-17.0.14
export PATH="$JAVA_HOME/bin:$PATH"
BASE_DIR="/home/stereonote/work/work"  # 请根据实际情况修改
INPUT_DIR="$BASE_DIR/input"
OUTPUT_DIR="$BASE_DIR/output"
TRIMMOMATIC_PATH="/home/stereonote/miniconda3/bin/trimmomatic"  # 请根据实际情况修改
SAMTOOLS_PATH="/home/stereonote/miniconda3/envs/samtools/bin/samtools"
PICARD_JAR="/home/stereonote/software/picard/build/libs/picard.jar"  # Picard JAR 文件路径
LOG_FILE="$OUTPUT_DIR/pipeline.log"
MAX_JOBS=2  # 最大并行任务数
THREADS=64

# 创建基础输出目录
mkdir -p "$INPUT_DIR" "$OUTPUT_DIR" || exit 1

# 检查是否提供了样本 ID 参数
if [ $# -eq 0 ]; then
    echo "Usage: $0 <sample_id1> <sample_id2> ..." | tee -a "$LOG_FILE"
    exit 1
fi

# 初始化日志，记录所有样本 ID
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting pipeline with samples: $@" | tee -a "$LOG_FILE"

# 定义样本处理函数
process_sample() {
    local SAMPLE_ID="$1"

    # 定义样本特定的输入和输出路径
    SAMPLE_INPUT_DIR="$INPUT_DIR/"
    SAMPLE_OUTPUT_DIR="$OUTPUT_DIR/$SAMPLE_ID"
    READ1="$SAMPLE_INPUT_DIR/${SAMPLE_ID}_1.fq.gz"
    READ2="$SAMPLE_INPUT_DIR/${SAMPLE_ID}_2.fq.gz"
    REF_FA="$SAMPLE_INPUT_DIR/${SAMPLE_ID}.fasta"
    TRIMMED_READ1="$SAMPLE_OUTPUT_DIR/trimmomatic/${SAMPLE_ID}_1.trimmed.fq.gz"
    TRIMMED_READ2="$SAMPLE_OUTPUT_DIR/trimmomatic/${SAMPLE_ID}_2.trimmed.fq.gz"
    SAM_FILE="$SAMPLE_OUTPUT_DIR/bwa/${SAMPLE_ID}.sam"
    SORTED_BAM="$SAMPLE_OUTPUT_DIR/bwa/${SAMPLE_ID}.sorted.bam"
    MARKDUP_BAM="$SAMPLE_OUTPUT_DIR/mark_duplicates/${SAMPLE_ID}.sorted.markdup.bam"
    GVCF_FILE="$SAMPLE_OUTPUT_DIR/vcf/${SAMPLE_ID}.g.vcf"

      # 检查输入文件是否存在
    if [ ! -f "$READ1" ] || [ ! -f "$READ2" ] || [ ! -f "$REF_FA" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Error: Input files not found for sample $SAMPLE_ID" | tee -a "$LOG_FILE"
        exit 1  # 在子shell中退出，不会影响其他样本
    fi
    # 创建样本特定的输出子目录
    mkdir -p "$SAMPLE_OUTPUT_DIR/fastqc" \
             "$SAMPLE_OUTPUT_DIR/trimmomatic" \
             "$SAMPLE_OUTPUT_DIR/bwa" \
             "$SAMPLE_OUTPUT_DIR/mark_duplicates" \
             "$SAMPLE_OUTPUT_DIR/result" \
             "$SAMPLE_OUTPUT_DIR/vcf" || exit 1

    # 1. 数据准备
    echo "### 1. 数据准备 ### - Sample: $SAMPLE_ID" | tee -a "$LOG_FILE"

    # 检查参考基因组索引（若无则构建）
    if [ ! -f "$REF_FA.fai" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Creating SAMtools index for $REF_FA..." | tee -a "$LOG_FILE"
        "$SAMTOOLS_PATH" faidx "$REF_FA" || exit 1
    fi
    if [ ! -f "$REF_FA.bwt" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Creating BWA index for $REF_FA..." | tee -a "$LOG_FILE"
        bwa index "$REF_FA" || exit 1
    fi
     # 创建序列字典
    if [ ! -f "$SAMPLE_INPUT_DIR/${SAMPLE_ID}.dict" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Creating Sequence Dictionary for $REF_FA..." | tee -a "$LOG_FILE"
        gatk CreateSequenceDictionary -R "$REF_FA" -O "$SAMPLE_INPUT_DIR/${SAMPLE_ID}.dict" || exit 1
    fi
    # 2. 质量控制
    echo "### 2. 质量控制 ### - Sample: $SAMPLE_ID" | tee -a "$LOG_FILE"

    # 原始数据质量评估
    if [ ! -f "$SAMPLE_OUTPUT_DIR/fastqc_Before/${SAMPLE_ID}_1_fastqc.html" ] || [ ! -f "$SAMPLE_OUTPUT_DIR/fastqc_Before/${SAMPLE_ID}_2_fastqc.html" ];
     then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running FastQC on raw data for $SAMPLE_ID..." | tee -a "$LOG_FILE"
        mkdir -p "$SAMPLE_OUTPUT_DIR/fastqc_Before"
        fastqc -t 16 -o "$SAMPLE_OUTPUT_DIR/fastqc_Before" "$READ1" "$READ2" # FastQC: 8 线程
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Skipping FastQC on raw data for $SAMPLE_ID (already done)." | tee -a "$LOG_FILE"
    fi

    # 接头去除和低质量过滤
    if [ ! -f "$TRIMMED_READ1" ] || [ ! -f "$TRIMMED_READ2" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running Trimmomatic for trimming $SAMPLE_ID..." | tee -a "$LOG_FILE"
        "$TRIMMOMATIC_PATH" PE -phred33 -threads $THREADS \
            -trimlog "$SAMPLE_OUTPUT_DIR/trimmomatic/${SAMPLE_ID}_trimmomatic.log" \
            "$READ1" "$READ2" \
            "$TRIMMED_READ1" "$SAMPLE_OUTPUT_DIR/trimmomatic/${SAMPLE_ID}_1.unpaired.fq.gz" \
            "$TRIMMED_READ2" "$SAMPLE_OUTPUT_DIR/trimmomatic/${SAMPLE_ID}_2.unpaired.fq.gz" \
            SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:50 || exit 1
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Skipping Trimmomatic for $SAMPLE_ID (already done)." | tee -a "$LOG_FILE"
    fi

    # 过滤后数据质量评估
     if [ ! -f "$SAMPLE_OUTPUT_DIR/fastqc_after_filt/${SAMPLE_ID}_1.trimmed_fastqc.html" ] || [ ! -f "$SAMPLE_OUTPUT_DIR/fastqc_after_filt/${SAMPLE_ID}_2.trimmed_fastqc.html" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running FastQC on trimmed data for $SAMPLE_ID..." | tee -a "$LOG_FILE"
        mkdir -p "$SAMPLE_OUTPUT_DIR/fastqc_after_filt"
        fastqc -t 16 -o "$SAMPLE_OUTPUT_DIR/fastqc_after_filt" \
            "$TRIMMED_READ1" \
            "$TRIMMED_READ2"  # FastQC: 8 线程
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Skipping FastQC on trimmed data for $SAMPLE_ID (already done)." | tee -a "$LOG_FILE"
    fi

    # 3. 序列比对
    echo "### 3. 序列比对 ### - Sample: $SAMPLE_ID" | tee -a "$LOG_FILE"

    # BWA 比对
    if [ ! -f "$SAM_FILE" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running BWA alignment for $SAMPLE_ID..." | tee -a "$LOG_FILE"
        bwa mem -t 32 -R "@RG\tID:${SAMPLE_ID}\tPL:UNKNOWN\tLB:library\tSM:${SAMPLE_ID}" \
            "$REF_FA" "$TRIMMED_READ1" \
            "$TRIMMED_READ2" \
            > "$SAM_FILE" || exit 1  # BWA: 20 线程
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Skipping BWA alignment for $SAMPLE_ID (already done)." | tee -a "$LOG_FILE"
    fi

    # SAM 转 BAM 并排序（管道优化）
    if [ ! -f "$SORTED_BAM" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Converting SAM to sorted BAM for $SAMPLE_ID..." | tee -a "$LOG_FILE"
        "$SAMTOOLS_PATH" view -S -b "$SAM_FILE" | \
            "$SAMTOOLS_PATH" sort -@ 20 -m 2G -O bam -o "$SORTED_BAM" - || exit 1  # Samtools sort: 20 线程
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Skipping SAM to sorted BAM conversion for $SAMPLE_ID (already done)." | tee -a "$LOG_FILE"
    fi

    # 构建BAM 索引
    if [ ! -f "$SORTED_BAM.bai" ];then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Indexing sorted BAM for $SAMPLE_ID..." | tee -a "$LOG_FILE"
        "$SAMTOOLS_PATH" index "$SORTED_BAM" || exit 1
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Skipping sorted BAM indexing for $SAMPLE_ID (already done)." | tee -a "$LOG_FILE"
    fi

    # 4. 统计分析
    echo "### 4. 统计分析 ### - Sample: $SAMPLE_ID" | tee -a "$LOG_FILE"

    # 计算测序深度
    if [ ! -f "$SAMPLE_OUTPUT_DIR/result/${SAMPLE_ID}_all_position_depth.txt" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Calculating depth for $SAMPLE_ID..." | tee -a "$LOG_FILE"
        "$SAMTOOLS_PATH" depth -a "$SORTED_BAM" \
            > "$SAMPLE_OUTPUT_DIR/result/${SAMPLE_ID}_all_position_depth.txt" || exit 1
    else
         echo "[$(date '+%Y-%m-%d %H:%M:%S')] Skipping depth calculation for $SAMPLE_ID (already done)." | tee -a "$LOG_FILE"
    fi

    # 计算覆盖度
     if [ ! -f "$SAMPLE_OUTPUT_DIR/result/${SAMPLE_ID}_coverage.txt" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Calculating coverage for $SAMPLE_ID..." | tee -a "$LOG_FILE"
        "$SAMTOOLS_PATH" coverage "$SORTED_BAM" \
            > "$SAMPLE_OUTPUT_DIR/result/${SAMPLE_ID}_coverage.txt" || exit 1
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Skipping coverage calculation for $SAMPLE_ID (already done)." | tee -a "$LOG_FILE"
    fi
    # 5. 标记重复序列
    echo "### 5. 标记重复序列 ### - Sample: $SAMPLE_ID" | tee -a "$LOG_FILE"

    if [ ! -f "$MARKDUP_BAM" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running Picard MarkDuplicates for $SAMPLE_ID..." | tee -a "$LOG_FILE"
        java -jar "$PICARD_JAR" MarkDuplicates \
            --INPUT "$SORTED_BAM" \
            --OUTPUT "$MARKDUP_BAM" \
            --METRICS_FILE "$SAMPLE_OUTPUT_DIR/mark_duplicates/${SAMPLE_ID}.markdup_metrics.txt" \
             || exit 1
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Skipping Picard MarkDuplicates for $SAMPLE_ID (already done)." | tee -a "$LOG_FILE"
    fi

    # 为标记重复后的 BAM 文件建立索引
     if [ ! -f "$MARKDUP_BAM.bai" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Indexing marked duplicates BAM for $SAMPLE_ID..." | tee -a "$LOG_FILE"
        "$SAMTOOLS_PATH" index "$MARKDUP_BAM" || exit 1
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Skipping marked duplicates BAM indexing for $SAMPLE_ID (already done)." | tee -a "$LOG_FILE"
    fi

    # 6. 变异检测
    echo "### 6. 变异检测 ### - Sample: $SAMPLE_ID" | tee -a "$LOG_FILE"

    if [ ! -f "$GVCF_FILE" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running GATK HaplotypeCaller for $SAMPLE_ID..." | tee -a "$LOG_FILE"
        gatk HaplotypeCaller \
            -R "$REF_FA" \
            -I "$MARKDUP_BAM" \
            -O "$GVCF_FILE" \
            --emit-ref-confidence GVCF \
            --ploidy 1 
            --native-pair-hmm-threads 32|| exit 1
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Skipping GATK HaplotypeCaller for $SAMPLE_ID (already done)." | tee -a "$LOG_FILE"
    fi
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished processing sample: $SAMPLE_ID" | tee -a "$LOG_FILE"
}

# 控制并行任务
for SAMPLE_ID in "$@"; do
    # 启动样本处理任务
    process_sample "$SAMPLE_ID" &

    # 检查当前运行的后台任务数
    while [ $(jobs -r | wc -l) -ge "$MAX_JOBS" ]; do
        sleep 1  # 等待 1 秒后重新检查
    done
done

# 等待所有样本处理完成
wait

echo "[$(date '+%Y-%m-%d %H:%M:%S')] ### 脚本执行完毕 ###" | tee -a "$LOG_FILE"