#!/bin/bash
# ==============================================
# NGS数据处理全流程脚本（精简统计版）
# ==============================================

set -euo pipefail

# --------------------------
# 1. 参数解析与初始化
# --------------------------
usage() {
    echo "Usage: $0 -r <reference.fa> -1 <R1.fastq.gz> -2 <R2.fastq.gz> [options]"
    echo "Options:"
    echo "  -r FILE    Reference genome (required)"
    echo "  -1 FILE    Fastq read 1 (required)"
    echo "  -2 FILE    Fastq read 2 (required)"
    echo "  -o DIR     Output directory (default: ./output)"
    echo "  -t INT     Threads for BWA/samtools (default: 4)"
    echo "  -m INT     Minimum match length (default: 80)"
    echo "  -f         Force overwrite existing files"
    echo "  -h         Show this help message"
    exit 1
}

# 默认参数
OUTPUT_DIR="./output"
THREADS=4
MIN_MATCH=80
FORCE_OVERWRITE=false

# 解析命令行参数
while getopts ":r:1:2:o:t:m:fh" opt; do
    case $opt in
        r) REF="$OPTARG" ;;
        1) INPUT_1_FQ="$OPTARG" ;;
        2) INPUT_2_FQ="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        m) MIN_MATCH="$OPTARG" ;;
        f) FORCE_OVERWRITE=true ;;
        h) usage ;;
        \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
        :) echo "Option -$OPTARG requires an argument" >&2; usage ;;
    esac
done

# 检查必需参数
if [[ -z "${REF:-}" || -z "${INPUT_1_FQ:-}" || -z "${INPUT_2_FQ:-}" ]]; then
    echo "Error: Missing required arguments!" >&2
    usage
fi

# 检查文件是否存在
check_file() {
    if [[ ! -f "$1" ]]; then
        echo "Error: File not found: $1" >&2
        exit 1
    fi
}

check_file "$REF"
check_file "$INPUT_1_FQ"
check_file "$INPUT_2_FQ"

# 创建目录结构
BWA_DIR="${OUTPUT_DIR}/BWA"
RESULTS_DIR="${OUTPUT_DIR}/results"
SPLIT_DIR="${OUTPUT_DIR}/split"
TMP_DIR="${OUTPUT_DIR}/tmp"

mkdir -p "${BWA_DIR}" "${RESULTS_DIR}" "${SPLIT_DIR}" "${TMP_DIR}"

# 记录开始时间
START_TIME=$(date +%s)
echo "========================================"
echo "NGS Pipeline Started at: $(date)"
echo "Reference: $REF"
echo "Input R1: $INPUT_1_FQ"
echo "Input R2: $INPUT_2_FQ"
echo "Output Dir: $OUTPUT_DIR"
echo "Threads: $THREADS"
echo "Min Match Length: $MIN_MATCH"
echo "Force Overwrite: $FORCE_OVERWRITE"
echo "========================================"

# --------------------------
# 辅助函数
# --------------------------

# 运行步骤并检查输出
run_step() {
    local step_name="$1"
    local output_file="$2"
    local command="$3"
    
    echo "[STEP] $step_name"
    if [[ $FORCE_OVERWRITE == true ]] || [[ ! -f "$output_file" ]]; then
        echo "Running step (output file missing or force mode)..."
        eval "$command"
        if [[ ! -f "$output_file" ]]; then
            echo "Error: Step failed to produce expected output: $output_file" >&2
            exit 1
        fi
        echo "Step completed successfully"
    else
        echo "Skipping step (output file exists: $output_file)"
    fi
    echo "----------------------------------------"
}

# --------------------------
# 2. 比对流程
# --------------------------
run_step "1/9 Building BWA index" "${REF}.bwt" \
    "bwa index '$REF'"

run_step "2/9 Running BWA alignment" "${BWA_DIR}/alignment.sam" \
    "bwa mem -t '$THREADS' '$REF' '$INPUT_1_FQ' '$INPUT_2_FQ' > '${BWA_DIR}/alignment.sam'"

run_step "3/9 Converting SAM to sorted BAM" "${BWA_DIR}/alignment.sorted.bam" \
    "samtools view -@ '$THREADS' -bS '${BWA_DIR}/alignment.sam' > '${BWA_DIR}/alignment.bam' && \
     samtools sort -@ '$THREADS' '${BWA_DIR}/alignment.bam' -o '${BWA_DIR}/alignment.sorted.bam' && \
     samtools index '${BWA_DIR}/alignment.sorted.bam' && \
     rm -f '${BWA_DIR}/alignment.sam' '${BWA_DIR}/alignment.bam'"

# --------------------------
# 3. 初始统计
# --------------------------
run_step "4/9 Generating initial statistics" "${RESULTS_DIR}/initial_stats.txt" \
    "{
        echo '===== Alignment Summary ====='
        samtools flagstat '${BWA_DIR}/alignment.sorted.bam'
        echo ''
        echo '===== Coverage Summary ====='
        samtools coverage '${BWA_DIR}/alignment.sorted.bam'
    } > '${RESULTS_DIR}/initial_stats.txt'"

# --------------------------
# 4. 基于CIGAR匹配长度过滤
# --------------------------
FILTERED_BAM="${BWA_DIR}/match${MIN_MATCH}_filtered.bam"
run_step "5/9 Filtering by CIGAR match length (>=${MIN_MATCH}bp)" "$FILTERED_BAM" \
    "samtools view -h '${BWA_DIR}/alignment.sorted.bam' | \
     awk -v min_len='$MIN_MATCH' '
     BEGIN {
         FS=\"\\t\"; OFS=\"\\t\";
         total=0; passed=0;
     }
     /^@/ { print; next }
     {
         total++;
         cigar=\$6;
         if (cigar == \"*\") next;
         
         total_M = 0;
         num = \"\";
         for (i=1; i<=length(cigar); i++) {
             c = substr(cigar, i, 1);
             if (c ~ /[0-9]/) {
                 num = num c;
             } else {
                 if (c == \"M\" && num != \"\") {
                     total_M += num;
                 }
                 num = \"\";
             }
         }
         
         if (total_M >= min_len) {
             print \$0;
             passed++;
         }
     }
     END {
         print \"Total reads:\", total > \"${RESULTS_DIR}/filter_stats.txt\";
         print \"Passed reads:\", passed > \"${RESULTS_DIR}/filter_stats.txt\";
         print \"Filtered reads:\", total-passed > \"${RESULTS_DIR}/filter_stats.txt\";
         print \"Filtering percentage:\", (total-passed)*100/total \"%\" > \"${RESULTS_DIR}/filter_stats.txt\";
     }' | samtools view -@ '$THREADS' -b -o '$FILTERED_BAM' -"

run_step "6/9 Processing filtered BAM" "${FILTERED_BAM%.bam}.sorted.bam" \
    "samtools sort -@ '$THREADS' '$FILTERED_BAM' -o '${FILTERED_BAM%.bam}.sorted.bam' && \
     samtools index '${FILTERED_BAM%.bam}.sorted.bam' && \
     {
         echo '===== Filtered Alignment Summary ====='
         samtools flagstat '${FILTERED_BAM%.bam}.sorted.bam'
         echo ''
         echo '===== Filtered Coverage Summary ====='
         samtools coverage '${FILTERED_BAM%.bam}.sorted.bam'
     } > '${RESULTS_DIR}/filtered_stats.txt'"

# --------------------------
# 5. 生成参考序列比对统计
# --------------------------
run_step "7/9 Generating reference statistics" "${RESULTS_DIR}/reference_stats.txt" \
    "{
        echo '===== Reference Sequence Statistics ====='
        echo 'Note: Only showing mapped reads count per reference'
        samtools idxstats '${FILTERED_BAM%.bam}.sorted.bam' | awk '\$1 != \"*\" {print \$1 \"\t\" \$3}'
    } > '${RESULTS_DIR}/reference_stats.txt'"

# --------------------------
# 6. 按参考序列拆分BAM
# --------------------------
echo "[8/9] Splitting BAM by reference sequence..."
REF_NAMES=$(samtools idxstats "${BWA_DIR}/alignment.sorted.bam" | awk '$1 != "*" && $3 > 0 {print $1}')

for ref in $REF_NAMES; do
    SPLIT_BAM="${BWA_DIR}/${ref}.bam"
    run_step "   Processing $ref" "$SPLIT_BAM" \
        "samtools view -@ '$THREADS' -b -h '${FILTERED_BAM%.bam}.sorted.bam' '$ref' > '$SPLIT_BAM' && \
         samtools index '$SPLIT_BAM'"
done

# --------------------------
# 7. 转换为fastq（用于下游分析）
# --------------------------
echo "[9/9] Converting back to FASTQ..."
for ref in $REF_NAMES; do
    SPLIT_BAM="${BWA_DIR}/${ref}.bam"
    BASE_NAME=$(basename "$SPLIT_BAM" .bam)
    
    R1_FASTQ="${SPLIT_DIR}/${BASE_NAME}_R1.fastq.gz"
    run_step "   Converting $ref to FASTQ" "$R1_FASTQ" \
        "samtools sort -@ '$THREADS' -n -o '${TMP_DIR}/${BASE_NAME}.nsorted.bam' '$SPLIT_BAM' && \
         samtools fastq -@ '$THREADS' -n -F 0x900 \
             -1 '${SPLIT_DIR}/${BASE_NAME}_R1.fastq.gz' \
             -2 '${SPLIT_DIR}/${BASE_NAME}_R2.fastq.gz' \
             -s '${SPLIT_DIR}/${BASE_NAME}_singletons.fastq.gz' \
             '${TMP_DIR}/${BASE_NAME}.nsorted.bam' && \
         rm -f '${TMP_DIR}/${BASE_NAME}.nsorted.bam'"
done

# --------------------------
# 8. 生成最终报告
# --------------------------
run_step "Generating final report" "${RESULTS_DIR}/final_report.txt" \
    "{
        echo '===== NGS Processing Summary ====='
        echo 'Processing completed: $(date)'
        echo ''
        echo '===== Input Parameters ====='
        echo \"Reference: $REF\"
        echo \"Read 1: $INPUT_1_FQ\"
        echo \"Read 2: $INPUT_2_FQ\"
        echo \"Min Match Length: ${MIN_MATCH}bp\"
        echo ''
        echo '===== Alignment Statistics ====='
        echo \"Original reads: \$(samtools view -c ${BWA_DIR}/alignment.sorted.bam)\"
        echo \"Filtered reads: \$(samtools view -c ${FILTERED_BAM%.bam}.sorted.bam)\"
        echo ''
        echo '===== Filtering Statistics ====='
        [ -f \"${RESULTS_DIR}/filter_stats.txt\" ] && cat \"${RESULTS_DIR}/filter_stats.txt\" || echo \"No filtering stats available\"
        echo ''
        echo '===== Output Files ====='
        echo \"BAM files: ${BWA_DIR}/\"
        echo \"FASTQ files: ${SPLIT_DIR}/\"
        echo \"Reports: ${RESULTS_DIR}/\"
    } > \"${RESULTS_DIR}/final_report.txt\""

# --------------------------
# 9. 清理和完成
# --------------------------
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo "========================================"
echo "Pipeline successfully completed!"
echo "Total processing time: $(($ELAPSED/60)) minutes $(($ELAPSED%60)) seconds"
echo ""
echo "===== Key Statistics ====="
echo "Original reads: $(samtools view -c ${BWA_DIR}/alignment.sorted.bam 2>/dev/null || echo 0)"
echo "Filtered reads: $(samtools view -c ${FILTERED_BAM%.bam}.sorted.bam 2>/dev/null || echo 0)"
[ -f "${RESULTS_DIR}/filter_stats.txt" ] && cat "${RESULTS_DIR}/filter_stats.txt"
echo ""
echo "Final report: ${RESULTS_DIR}/final_report.txt"
echo "========================================"

exit 0