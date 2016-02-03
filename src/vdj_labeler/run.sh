#!/bin/bash
(cd ../../; ./prepare_cfg
make vdj -j8
build/release/bin/vdj_labeler configs/vdj_labeler/configs.info)
