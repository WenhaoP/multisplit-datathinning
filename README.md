# multisplit-datathinning

Count splitting, which is proposed by Neufeld et al. (2022), is a framework that provides valid inference by avoiding using the same data for both latent variable estimation and differential expression analysis under a Poisson assumption. Later, Neufeld et al. (2023) generalize this idea to convolution-closed
distributions and propose the data thinning algorithm. In this project, we aim to apply the rank-transformed subsampling algorithm, which is proposed by Guo and Shah (2023), to alleviate thereplication issue of and improve the power and stability of the hypothesis test carried out by the data thinning while maintaining the control of the Type 1 error.

We first directly apply the rank-transformed subsampling algorithm to the data thinning algorithm and use the combined algorithm to redo the simulation study in (Neufeld et al., 2022)’s paper. If thereis any numerical improvement, we will start the theory development of the combined algorithm. We are interested to see if the theoretical results from Guo and Shah (2023)’s paper still apply to the combined algorithm.

## Plan

* Code up the unparalleled version of the algorithm
* Run the algorithm on the small-sized simulation data and check the numerical output to validate the algorithm
* Plot the type-1 error and power 
* Run the algorithm on the original-sized simulation data and plot the plots
* Parallelize the algorithm
    * Might actually need to do this before the (iv) step 

## Repository Structure 

* `algo.R` contains the combined algorithm.
* `multi.R` contains the helper functions directly copied from the public GitHub [repo](https://github.com/richardkwo/MultiSplit) of (Guo and Shah,2023).
* `main.R` contains the code for running the simulation study.

## References
Guo, F. R. and Shah, R. D. (2023). Rank-transformed subsampling: inference for multiple data splitting and exchangeable p-values. arXiv preprint arXiv:2301.02739

Neufeld, A., Dharamshi, A., Gao, L. L., and Witten, D. (2023). Data thinning for convolution-closed distributions. arXiv preprint arXiv:2301.07276

Neufeld, A., Gao, L. L., Popp, J., Battle, A., and Witten, D. (2022). Inference after latent variable estimation for single-cell rna sequencing data. arXiv preprint arXiv:2207.00554

