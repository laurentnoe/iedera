# iedera

(from http://bioinfo.lifl.fr/yass/iedera.php) 

Iedera is a program to select and design subset seed and vectorized subset seed patterns. Spaced seeds and transition constrained spaced seeds can be perfectly represented in the subset seed model. Moreover, BLASTP-like seeds, and more generally Vector seeds (seeds with cumulative score constraint) can also be represented by vectorized subset seeds (which is a subset seed with an additional score constraint).

Iedera is applied to both lossy and lossless seed design. It is already used to design spaced seeds templates for read mappers (see below), to design subset seed templates for protein sequences or nucleic sequences. Recent applications also include sub profile-HMM read mapping problem, alignment-free distances and SVM string kernels with a new coverage criterion, and parameter-free spaced seed design.

This tool is experimental and is provided for research purposes only !! It natively supports two probabilistic models to describe the alignment sequence distribution (Bernoulli and Markov) and other models (as HMM) can be represented by probabilistic automata (see below).
