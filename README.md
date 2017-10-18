# Simple MIMO Simulator in Julia

Please visit the [VIP Group](http://vip.ece.cornell.edu) website for recent publications and information about our research.

### Important information 

If you are thinking of contacting us, please *do not* e-mail the author to ask for download instructions, installation guidelines, or the simulator itself. The code itself is well-documented and provides the essential information about the code. Note that we will NOT help to debug user-generated code that was not included in the provided software package. If, however, you notice a bug in our code, please be so kind to contact the [Author](mailto:studer@cornell.edu). 

The software package is supplied "as is," without any accompanying support services, maintenance, or future updates. We make no warranties, explicit or implicit, that the software contained in this package is free of error or that it will meet your requirements for any particular application. It should not be relied on for any purpose where incorrect results could result in loss of property, personal injury, liability or whatsoever. If you do use our software for any such purpose, it is at your own risk. The authors disclaim all liability of any kind, either direct or consequential, resulting from your use of these programs. 

### Overview

Multiple-input multiple-output (MIMO) technology in combination with spatial multiplexing has established itself as a way of significantly increasing the spectral efficiency compared to single-antenna wireless communication systems. One of the main drawbacks of MIMO communication is the computational complexity of optimal data detection. Hence, most practical implementations for MIMO technology rely on low-complexity algorithms that provide sub-optimal performance. During the last decade, a plethora of algorithms and corresponding hardware designs have been proposed in the literature.

This software provides a simple yet modular simulation framework for MIMO wireless communication suitable for less experienced undergraduate and graduate students. The code is written in [Julia](https://julialang.org/) and can be extended to incorporate various data detection algorithms. The simulator generates error-rate performance curves for hard-output data detectors, such as zero forcing (ZF), linear minimum mean-square error (LMMSE), and maximum likelihood (ML) detection via sphere detection. 

### Package details

This repository contains a single-file Julia simulator that performs Monte-Carlo simulations to extract error-rate vs. signal-to-noise ratio (SNR) curves. The simulator only supports hard-output MIMO detectors, but is set up such that you can add your own extensions (e.g., algorithms, channel models, codes, etc.). The code is written by C. Studer with help from E. Gonultas, and is available for free trial, non-commercial research or education purposes, and for non-profit organizations. If you plan on using the code or parts thereof for commercial purposes or if you intend to re-distribute the code or parts thereof, you must contact the [Author](mailto:studer@cornell.edu). If you are using the code or parts thereof for your scientific work, you *must* provide a reference to this repository. 

### How to run a simple simulation

Open the Julia console and simply type

```sh
include("simpleMIMOsim.jl")
```

or execute the code in Atom/Juno by pressing `Ctrl+Shift+Enter` (or `Command+Shift+Enter` in MacOS). This starts a simulation with predefined parameters and data detector, i.e., 4x4 MIMO system with 16-QAM for an i.i.d. Rayleigh fading scenario. You can provide your own system and simulation parameters by either modifying the code or by passing your own `par`-type (see the sourcecode for an example). We highly recommend to inspect the code step-by-step in order to get a detailed understanding of the simulator. 

### Version history
* Version 0.1: studer@cornell.edu - initial version for GitHub release

