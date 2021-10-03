# G2 - a package for comparing counts
This package implements the G^2 test for comparing counts as described in [this paper](https://aclanthology.org/J93-1003/).

This is useful because it doesn't depend on an assumption of normal distribution that makes the normal Pearson's $\chi^2$ 
test drastically over-estimate how interesting a situation is. This commonly happens when, for instance, comparing the 
frequencies of words in different contexts or deciding whether a user action is an interesting indicator that of a 
certain purchase behavior.

[![Build Status](https://github.com/tdunning/G2.jl/workflows/CI/badge.svg)](https://github.com/tdunning/G2.jl/actions/workflows/run-tests.yml?query=branch%3Amain)
