#!/bin/bash
cd ../../; ./prepare_cfg
make;
build/release/bin/vdj_labeler configs/vdj_labeler/configs.info;
#build/release/bin/test_accurate_vdj_labeling;
