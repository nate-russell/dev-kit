# Sequential 



## Off-Policy Evaluation for BO: Suragate Models + Aquisition Functions



#### Resources
 - https://github.com/uber/bayesmark
 - https://botorch.org/
 - https://ax.dev/
 - https://scikit-optimize.github.io/stable/index.html
 - https://gpytorch.ai/
 - https://emukit.github.io/



#### Benchmarking
- [Revisiting Bayesian Optimization in the light of the COCO benchmark](https://arxiv.org/pdf/2103.16649.pdf)
- [Bayesian Optimization: Model Comparison With Different Benchmark Functions](https://ieeexplore.ieee.org/document/9707005)
- [Efficient tuning of online systems using Bayesian optimization](https://research.facebook.com/blog/2018/09/efficient-tuning-of-online-systems-using-bayesian-optimization/)
- [Benchmarking the performance of Bayesian optimization across multiple experimental materials science domains](https://www.nature.com/articles/s41524-021-00656-9)
- [Benchmarking Bayesian optimisation](https://l2s.centralesupelec.fr/wp-content/uploads/uqsay/uqsay25_slides_vpicheny.pdf)
- [Benchmarking surrogate-based optimisation algorithms on expensive black-box functions](https://www.sciencedirect.com/science/article/pii/S1568494623007627)
#### Off-Policy Eval
- [A Review of Off-Policy Evaluation in Reinforcement Learning](https://arxiv.org/abs/2212.06355)
- [ParallelrobustBayesianoptimization withoff-policyevaluations]()
- [BAYESIANCOUNTERFACTUALMEANEMBEDDINGSAND OFF-POLICYEVALUATION]()
- [Offline Policy Selection under Uncertainty]()
- [UniversalOff-PolicyEvaluation]()
- [Active Offline Policy Selection](https://arxiv.org/abs/2106.10251)
- []()




### Problem Settings
1. Start off with K samples of some unknown but consistent policy, we want to estimate wether or not switching policies is a good idea. real-world example: started experimenting with human expert plan, now wondering if it is time to switch to a different BO Policy for another K' samples
2. Same as 1 but we know the functional form of the sampling policy, maybe it was random or GP-UCB
