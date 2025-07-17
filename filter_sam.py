# filter_sam.py
#
# 该脚本用于过滤SAM文件，不依赖任何外部库。
# 筛选条件如下：
# 1. Read 必须是正确配对 (FLAG 0x2)。
# 2. Read 必须是主要比对 (非 FLAG 0x100)。
# 3. Read 的 CIGAR 字符串中 'M' 操作的总长度必须大于120。
# 4. 一对 reads (pair) 中的两条都必须满足以上所有条件才会被保留。

import sys
import os
import re

# --- 使用说明 ---
#
# 1. 准备输入文件:
#    为了让脚本高效运行且内存占用低，建议先按 read name 对SAM/BAM文件排序。
#    如果输入是BAM: samtools sort -n your_input.bam | samtools view -h > sorted_by_name.sam
#    如果输入是SAM: sort -k1,1 your_input.sam > sorted_by_name.sam
#
# 2. 运行脚本:
#    python filter_sam.py <input.sam> <output.sam>
#
# --- 脚本开始 ---

def filter_sam_file():
    # --- 配置 ---
    if len(sys.argv) != 3:
        print("错误: 参数不足。")
        print("用法: python filter_sam.py <input.sam> <output.sam>")
        sys.exit(1)

    input_sam_path = sys.argv[1]
    output_sam_path = sys.argv[2]
    min_m_length = 100
    # --- 配置结束 ---

    # 检查输入文件是否存在
    if not os.path.exists(input_sam_path):
        print(f"错误: 输入文件未找到: {input_sam_path}")
        sys.exit(1)

    # 该字典用于存储已通过筛选并等待其配对 read 的 read (存储的是完整的行文本)
    pending_reads = {}

    print("开始处理 SAM 文件...")
    line_count = 0
    written_pairs = 0

    # 打开输入和输出文件
    with open(input_sam_path, 'r') as infile, open(output_sam_path, 'w') as outfile:
        # 遍历输入文件中的每一行
        for line in infile:
            line_count += 1
            if line_count % 5000000 == 0:
                print(f"  ...已处理 {line_count:,} 行。")

            # 1. 直接写入头部行
            if line.startswith('@'):
                outfile.write(line)
                continue

            # 2. 解析比对行
            fields = line.strip().split('\t')
            if len(fields) < 11:  # 跳过格式不正确的行
                continue

            try:
                qname = fields[0]
                flag = int(fields[1])
                cigar = fields[5]
            except (ValueError, IndexError):
                continue # 如果FLAG不是数字或列数不够，则跳过

            # 3. 执行筛选
            # 条件1: 检查 FLAG (位运算)
            # 必须是正确配对 (flag & 2) 且不是非主要比对 !(flag & 256)
            is_proper_pair = (flag & 2) != 0
            is_secondary = (flag & 256) != 0
            if not is_proper_pair or is_secondary:
                continue

            # 条件2: 检查 CIGAR 中 'M' 的总长度
            if cigar == '*': # 跳过没有CIGAR字符串的未比对read
                continue
            
            # 使用正则表达式查找所有 '...M' 的模式并加总
            m_lengths = re.findall(r'(\d+)M', cigar)
            total_m = sum(int(length) for length in m_lengths)

            if total_m <= min_m_length:
                continue

            # 4. 如果 read 通过了所有筛选，则处理配对逻辑
            if qname in pending_reads:
                # 这对 read 的另一条已经先被找到并通过了筛选
                # 将两条 reads 的行文本都写入输出文件
                outfile.write(pending_reads[qname])
                outfile.write(line)
                written_pairs += 1
                
                # 从字典中删除该条目以释放内存
                del pending_reads[qname]
            else:
                # 这是这对 read 中第一条通过筛选的
                # 将它的完整行文本存起来，等待它的配对 read
                pending_reads[qname] = line

    # --- 总结 ---
    print("-" * 30)
    print("过滤完成。")
    print(f"总共处理的行数: {line_count:,}")
    print(f"写入输出文件的配对数量: {written_pairs:,} ({written_pairs * 2:,} 条 reads)")
    if len(pending_reads) > 0:
        print(f"警告: 有 {len(pending_reads):,} 条 reads 通过了筛选，但它们的配对 read 未通过。")
    print(f"输出文件已创建: {output_sam_path}")

if __name__ == "__main__":
    filter_sam_file()
