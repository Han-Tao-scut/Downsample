#!/bin/bash
# ==============================================
# NGS数据处理脚本（比对与基本统计）
# ==============================================

set -euo pipefail # 命令失败时退出, 使用未定义变量时报错, pipe失败时退出

# --------------------------
# 1. 参数解析与初始化
# --------------------------
usage() {
    echo "Usage: $0 -r <reference.fa> -1 <R1.fastq.gz> -2 <R2.fastq.gz> [options]"
    echo "Options:"
    echo "  -r FILE    参考基因组 FASTA 文件 (必需)"
    echo "  -1 FILE    FASTQ read 1 文件 (必需)"
    echo "  -2 FILE    FASTQ read 2 文件 (必需)"
    echo "  -o DIR     输出目录 (默认: ./output)"
    echo "  -t INT     用于 BWA/samtools 的线程数 (默认: 4)"
    echo "  -f         强制覆盖已存在的输出文件"
    echo "  -h         显示此帮助信息"
    exit 1
}

# 默认参数
OUTPUT_DIR="/data/output"
THREADS=4
FORCE_OVERWRITE=false
REF=""
INPUT_1_FQ=""
INPUT_2_FQ=""
SAMTOOLS_PATH="/home/stereonote/samtools-1.21/samtools"

# 解析命令行参数 (移除了 :m)
while getopts ":r:1:2:o:t:fh" opt; do
    case $opt in
        r) REF="$OPTARG" ;;
        1) INPUT_1_FQ="$OPTARG" ;;
        2) INPUT_2_FQ="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        # m) MIN_MATCH="$OPTARG" ;; # 移除 -m 参数处理
        f) FORCE_OVERWRITE=true ;;
        h) usage ;;
        \?) echo "无效选项: -$OPTARG" >&2; usage ;;
        :) echo "选项 -$OPTARG 需要一个参数" >&2; usage ;;
    esac
done
# 检查必需参数
if [[ -z "$REF" || -z "$INPUT_1_FQ" || -z "$INPUT_2_FQ" ]]; then
    echo "错误: 缺少必需的参数 (-r, -1, -2)!" >&2
    usage
fi

# 检查文件是否存在
check_file() {
    if [[ ! -f "$1" ]]; then
        echo "错误: 文件未找到: $1" >&2
        exit 1
    fi
}

check_file "$REF"
check_file "$INPUT_1_FQ"
check_file "$INPUT_2_FQ"

# 创建目录结构 (移除了 SPLIT_DIR 和 TMP_DIR)
BWA_DIR="${OUTPUT_DIR}/BWA"
RESULTS_DIR="${OUTPUT_DIR}/results"

mkdir -p "${BWA_DIR}" "${RESULTS_DIR}"

# 记录开始时间
START_TIME=$(date +%s)
echo "========================================"
echo "NGS Pipeline Started at: $(date)"
echo "Reference: $REF"
echo "Input R1: $INPUT_1_FQ"
echo "Input R2: $INPUT_2_FQ"
echo "Output Dir: $OUTPUT_DIR"
echo "Threads: $THREADS"
# echo "Min Match Length: $MIN_MATCH" # 移除未使用参数的显示
echo "Force Overwrite: $FORCE_OVERWRITE"
echo "========================================"

# --------------------------
# 辅助函数
# --------------------------

# 运行步骤并检查输出
run_step() {
    local step_name="$1"
    local output_file="$2" # 主要输出文件，用于检查存在性
    local command_str="$3" # 命令字符串

    echo "----------------------------------------"
    echo "[STEP] $step_name"
    # 检查主要输出文件是否存在，以及是否未设置强制覆盖
    if [[ $FORCE_OVERWRITE == false ]] && [[ -f "$output_file" ]]; then
        echo "跳过步骤 (输出文件已存在: $output_file)"
    else
        if [[ $FORCE_OVERWRITE == true ]] && [[ -f "$output_file" ]]; then
             echo "运行步骤 (强制模式，覆盖现有输出)..."
        else
             echo "运行步骤 (输出文件不存在或强制模式)..."
        fi
        echo "命令: $command_str"
        # 使用eval执行命令，以支持管道等复杂语法
        if eval "$command_str"; then
            # 再次检查主要输出文件是否确实已生成
            if [[ ! -f "$output_file" ]]; then
                 echo "错误: 步骤完成但未产生预期的输出文件: $output_file" >&2
                 exit 1
            fi
            echo "步骤成功完成。"
        else
            # eval 返回非0表示命令执行失败
            echo "错误: 步骤执行失败，退出码 $?。" >&2
            exit 1
        fi
    fi
}

# --------------------------
# 2. 比对流程
# --------------------------
# 更新步骤计数为总共5步
run_step "1/5 构建 BWA 索引" "${REF}.bwt" \
    "bwa index '$REF'"

# 中间 SAM 文件，将在下一步被处理和删除
SAM_TEMP_FILE="${BWA_DIR}/alignment.sam"
run_step "2/5 运行 BWA 比对" "$SAM_TEMP_FILE" \
    "bwa mem -t '$THREADS' '$REF' '$INPUT_1_FQ' '$INPUT_2_FQ' > '$SAM_TEMP_FILE'"

# Sorted BAM 文件是此步骤的主要最终产物
SORTED_BAM_FILE="${BWA_DIR}/alignment.sorted.bam"
SORTED_BAM_INDEX="${SORTED_BAM_FILE}.bai" # Index file also created
run_step "3/5 SAM 转为 Sorted BAM 并索引" "$SORTED_BAM_FILE" \
    "$SAMTOOLS_PATH view -@ '$THREADS' -bS '$SAM_TEMP_FILE' | \
     $SAMTOOLS_PATH sort -@ '$THREADS' -o '$SORTED_BAM_FILE' - && \
     $SAMTOOLS_PATH index -@ '$THREADS' '$SORTED_BAM_FILE' && \
     rm -f '$SAM_TEMP_FILE'"
     # Note: Removed the intermediate .bam file step for efficiency

# --------------------------
# 3. 初始统计
# --------------------------
STATS_FILE="${RESULTS_DIR}/initial_stats.txt"
run_step "4/5 生成初始统计 (flagstat, coverage)" "$STATS_FILE" \
    "{ \
        echo '===== 比对摘要 (flagstat) ====='; \
        $SAMTOOLS_PATH flagstat -@ '$THREADS' '$SORTED_BAM_FILE'; \
        echo ''; \
        echo '===== 覆盖度摘要 (coverage) ====='; \
        $SAMTOOLS_PATH coverage '$SORTED_BAM_FILE'; \
    } > '$STATS_FILE'"

DEPTH_FILE="${RESULTS_DIR}/depth_origin.txt"
run_step "5/5 计算测序深度 (depth)" "$DEPTH_FILE" \
    "$SAMTOOLS_PATH depth -@ '$THREADS' -a  '$SORTED_BAM_FILE' > '$DEPTH_FILE'"

# --------------------------
# 4. 完成
# --------------------------
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
echo "========================================"
echo "所有步骤完成。"
echo "Pipeline Finished at: $(date)"
echo "总执行时间: ${DURATION} 秒。"
echo "输出文件位于: $OUTPUT_DIR"
echo "========================================"

exit 0
