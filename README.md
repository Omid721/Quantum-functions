# Quantum-Inspired Algorithms Using Tensor Networks
This repository employs **Tensor Network** approaches to develop quantum-inspired algorithms that emulate the behavior of quantum computing algorithms. Specifically, the examples provided focus on the use of Tensor Networks for calculating Matsubara sums, a key concept in quantum many-body physics. 


**Overview**

Tensor Networks offer a powerful framework to efficiently represent and manipulate high-dimensional tensors, which are essential for simulating quantum systems with low entangelment. In this repository, we use these networks to mimic certain quantum computing algorithms, allowing for classical simulations of quantum-inspired algorithms. These techniques are widely applicable in condensed matter physics and quantum information theory and more.

The examples here demonstrate the use of Tensor Network methods for calculating Matsubara sums, an important task in finite-temperature quantum field theory and many-body physics. (For more details see the [Tensor Network Guide](Tensor_Network_2.pdf))

**Dependencies**

To run the code, the following Julia packages are required for handling tensor contractions and decompositions:

•	ITensor.jl: A flexible library for implementing tensor network algorithms.

•	TensorCrossInterpolation.jl: A package designed to optimize tensor cross interpolation (making MPS for given function).

Make sure you have these dependencies installed in your Julia environment.

Install the required packages in Julia by running:

    using Pkg
    Pkg.add("ITensor")
    Pkg.add("TensorCrossInterpolation")

To run the example codes, you also need to include the required functions by adding the 

    include("functions.jl")

Repository structure

    ├── src/
    │   ├── functions.jl        
    ├── examples/
    │   ├── code_examples         
    │   ├── Notebook_examples        
    ├── test/
    │   ├── test_functions.jl  
    ├── README.md               
    ├── LICENSE                 
    ├── .gitignore              



