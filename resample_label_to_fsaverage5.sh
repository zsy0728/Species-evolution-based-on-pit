#!/bin/bash

# 设置变量
# INPUT_METRIC="/home/songyao/data/cross-species/BBW_BigData/spaces/tpl-fsaverage/tpl-fsaverage_hemi-L_den-164k_desc-layer1.shape.gii"
# OUTPUT_METRIC="/home/songyao/data/cross-species/BBW_BigData_results/tpl-fsaverage_hemi-L_den-10k_desc-layer1.shape.gii"

# SPHERE_SRC="/home/songyao/freesurfer/subjects/fsaverage/surf/lh.sphere.reg.fsaverage.surf.gii"
# SPHERE_TRG="/home/songyao/freesurfer/subjects/fsaverage5/surf/lh.sphere.reg.fsaverage5.surf.gii"

# AREA_SRC="/home/songyao/freesurfer/subjects/fsaverage/surf/lh.white.fsaverage.surf.gii"
# AREA_TRG="/home/songyao/freesurfer/subjects/fsaverage5/surf/lh.white.fsaverage5.surf.gii"

# 执行 resample
# wb_command -metric-resample "$INPUT_METRIC" "$SPHERE_SRC" "$SPHERE_TRG" ADAP_BARY_AREA "$OUTPUT_METRIC" -area-surfs "$AREA_SRC" "$AREA_TRG"

#!/bin/bash

# 输入与输出目录
INPUT_DIR="/home/songyao/data/cross-species/BBW_BigData/spaces/tpl-fsaverage"
OUTPUT_DIR="/home/songyao/data/cross-species/BBW_BigData_results"

# 重采样所需的表面文件
SPHERE_SRC="/home/songyao/freesurfer/subjects/fsaverage/surf/lh.sphere.reg.fsaverage.surf.gii"
SPHERE_TRG="/home/songyao/freesurfer/subjects/fsaverage5/surf/lh.sphere.reg.fsaverage5.surf.gii"
AREA_SRC="/home/songyao/freesurfer/subjects/fsaverage/surf/lh.white.fsaverage.surf.gii"
AREA_TRG="/home/songyao/freesurfer/subjects/fsaverage5/surf/lh.white.fsaverage5.surf.gii"

# 查找所有匹配的文件（包含164k且后缀为.shape.gii/.curv/.label.gii）
find "$INPUT_DIR" -type f \( -name "*164k*.shape.gii" -o -name "*164k*.curv" -o -name "*164k*.label.gii" \) | while read -r INPUT_METRIC; do
    # 获取文件名
    FILENAME=$(basename "$INPUT_METRIC")
    
    # 替换164k为10k
    OUTPUT_METRIC="${FILENAME/164k/10k}"
    
    # 拼接输出路径
    OUTPUT_PATH="$OUTPUT_DIR/$OUTPUT_METRIC"

    echo "Resampling: $FILENAME → $OUTPUT_METRIC"

    # 执行重采样
    wb_command -metric-resample "$INPUT_METRIC" "$SPHERE_SRC" "$SPHERE_TRG" ADAP_BARY_AREA "$OUTPUT_PATH" -area-surfs "$AREA_SRC" "$AREA_TRG"
done
