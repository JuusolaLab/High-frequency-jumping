# Model
# Overview
This repository provides a Musca (fly) visual system model that simulates photoreceptor and Large Monopolar Cell (LMC) responses to visual stimuli.
As an example, this implementation reproduces Figure 4b from the associated study.

# Hardware requirements

The model requires a computer with sufficient RAM to store simulation variables.
(Exact memory requirements depend on stimulus size and simulation parameters.)

# Software requirements

The model has been tested and verified using:
-MATLAB R2023b
-Parallel Computing Toolbox
-Statistics and Machine Learning Toolbox

Earlier versions of MATLAB may also work but have not been officially tested.

# Installation
1. Download or clone this repository
2. Run BurstyCalculationExample.m

# Example

BurstyCalculationExample.m
  -Demonstrates how to simulate both photoreceptor and LMC responses used in Figure 4b.

Bursty_stim_photo_hres_200.mat
  -Contains the example bursty stimulus input used in the model.
