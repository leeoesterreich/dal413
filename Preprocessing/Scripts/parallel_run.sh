#!/bin/bash

# 目标文件夹路径
TARGET_DIR="/bgfs/alee/LO_LAB/Personal/Daisong/Test_sample/2_trimmed_file/paired"

cd "$TARGET_DIR" || exit

# 储存前9个字符相同的文件对
declare -A file_pairs

# 遍历当前文件夹中的所有文件
for file in *; do
    if [[ -f "$file" ]]; then
        prefix=$(echo "$file" | cut -c1-9)
        if [[ -n "${file_pairs[$prefix]}"  && "${file_pairs[$prefix]: -1}" == "z" ]]; then
            # 如果已经找到一个相同前缀的文件，则传递给另一个脚本
            other_file="${file_pairs[$prefix]}"
            echo "run alignment.sh, with parameters: $other_file, $file"
            sbatch /bgfs/alee/LO_LAB/Personal/Daisong/Test_sample/Scripts/alignment_copy.sh "$other_file" "$file" > alignment_output_AAA.log 2>&1
            echo "already run with parameters: $other_file, $file"
            # 清除匹配记录
            unset file_pairs["$prefix"]
        else
            # 记录当前文件
            file_pairs["$prefix"]="$file"
        fi
    fi
done


wait
echo "All finished"