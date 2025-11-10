#!/bin/bash
# Script to compile, upload, and monitor simple.ino

set -e  # Exit on error

SKETCH_DIR="/Users/tatestaples/Code/Automatic-Laser-Alignment"
SKETCH_NAME="simple"
PORT="/dev/cu.usbmodem1101"
BOARD="arduino:avr:mega"
BAUD_RATE="9600"

echo "=== Running ${SKETCH_NAME}.ino ==="
echo "Copying sketch..."
mkdir -p "${SKETCH_DIR}/${SKETCH_NAME}"
cp "${SKETCH_DIR}/${SKETCH_NAME}.ino" "${SKETCH_DIR}/${SKETCH_NAME}/${SKETCH_NAME}.ino"

echo "Compiling..."
arduino-cli compile --fqbn ${BOARD} "${SKETCH_DIR}/${SKETCH_NAME}"

echo "Uploading..."
arduino-cli upload -p ${PORT} --fqbn ${BOARD} "${SKETCH_DIR}/${SKETCH_NAME}"

echo "Starting serial monitor at ${BAUD_RATE} baud..."
echo "Press CTRL-C to exit"
arduino-cli monitor -p ${PORT} -c baudrate=${BAUD_RATE}
