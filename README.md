# Branch Predictor Simulator

## Overview
This project implements a **Branch Predictor Simulator** that evaluates the performance of different branch prediction techniques, including **Gshare, Bimodal, and Hybrid predictors**. The simulator processes branch instruction traces and determines the effectiveness of each prediction strategy in reducing misprediction rates.

## Features
- Implements **Gshare, Bimodal, and Hybrid** branch prediction algorithms.
- Simulates branch prediction behavior on real-world traces.
- Provides **detailed performance metrics**, such as prediction accuracy and misprediction rate.
- Configurable predictor parameters for experimentation.
- Efficient implementation for fast execution and analysis.

## Installation
### Prerequisites
Ensure you have the following dependencies installed:
- C/C++ compiler (e.g., GCC, Clang)
- Python (if used for analysis or trace generation)
- Make (optional, for build automation)

### Build Instructions
```sh
# Clone the repository
git clone https://github.com/Sachin-Saravana-Kumar/branch-predictor
cd branch-predictor-simulator

# Compile the simulator
make  # If using a Makefile
# OR manually compile using g++
g++ -o sim_bp sim_bp.cc -O2
```

## Usage
### Running the Simulator
```sh
./sim_bp <PREDICTOR_TYPE> <PARAMS> <trace_file>
```
Example:
```sh
.sim_bp bimodal 1024 traces/branch_trace.trace
./sim_bp gshare 1024 10 traces/branch_trace.trace
./sim_bp hybrid 10 1024 10 1024 traces/branch_trace.trace
```

### Command-Line Arguments
- `<PREDICTOR_TYPE>`: Type of branch predictor (`bimodal`, `gshare`, or `hybrid`).
- `<PARAMS>`:
  - **Bimodal**: `<M2>` (Table size for Bimodal predictor).
  - **Gshare**: `<M1>` (Table size) and `<N>` (Number of bits used for history).
  - **Hybrid**: `<K>` (Chooser table size), `<M1>` (Gshare table size), `<N>` (History bits), and `<M2>` (Bimodal table size).
- `<trace_file>`: Path to the branch instruction trace file.

### Sample Output
```
COMMAND
./branch_predictor gshare 1024 10 traces/branch_trace.trace
traces.....
Branch Predictions: 500000
Mispredictions: 32000
Prediction Accuracy: 93.6%
```

## Branch Prediction Techniques
### Bimodal Predictor
- Uses a simple direct-mapped table of 2-bit saturating counters.
- Fast and simple but can struggle with certain branch patterns.

### Gshare Predictor
- Uses **global branch history** XORed with the PC index to improve correlation.
- More effective than bimodal for complex branch patterns.

### Hybrid Predictor
- Combines **Gshare and Bimodal** predictors using a selector table.
- Dynamically selects the best predictor for each branch.

## Performance Evaluation
You can analyze the simulator's performance using different branch traces and predictor configurations. The project includes:
- **Trace files** for testing different branch prediction patterns.
- **Python scripts** for performance visualization.

