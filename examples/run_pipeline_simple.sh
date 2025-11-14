#!/bin/bash
#
# Simple wrapper to run Immcantation pipeline with correct environment
#

# Activate venv
source /home/user/immcantation_env/bin/activate

# Run pipeline
python /home/user/sadie/examples/run_immcantation_pipeline.py "$@"
