#!/bin/bash
set -eux

echo "Building package"
R CMD build --no-manual .
