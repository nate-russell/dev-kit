# Optimal Subset Selection for Multi-Criteria experiments
 
## Intro
Experimentation is a central pillar of science. In the physical and life sciences, their design is often one of the most important decisions made by any lab/organization. They are often expensive, time consuming, complex and even if you do everything right, they often yield uncertain results. Rigourous statistical analysis is often applied to experimental results but too often that same level of rigour is missing from the selection/design process. 

By 1) thinking critically about the goals of an experiment 2) expressing these goals mathematically and 3) treating experimental design as an optimization problem  we can make tractable, defensible and optimal experimental designs quickly and reproducibly. By combining ideas from bayesian optimization, active learning, and optimal experimental design under the umbrella of multi-objective optimization we can understand trade-offs between designs and make informed decisions that satisfy multiple goals and constraints. 

This work focusses on subset selection scenarios that are often encoutered in perturbational assays/screens for hit finding, lead op, target discovery, probing biological systems and machine learning model development. Given a set of items S, find subset E &sub; S that maximizes obejctives F<sub>1,2,...k</sub>(E). We could be picking molecular probes to go into a hit finding assay, designing perturbations for a lead optimiztion campaign, doing a chemogenomic screen looking to learn an MOA, etc.  

#### Project Goals: 
1) demontstrate the value of treating experimental design and screening library selection as an optimization problem and 
2) provide tools to find and understand ~optimal solutions in a computatioanlly tractable manner

#### Notation/terminolgy:
- *Experiments* take many forms. There is the traditional control variables and levels viewpoint but in this work we focus on sub-set selection experiemnts. Examples) picking N molecules for a screening assay OR selecting cell line + CRISPER combinations to understand gene function OR picking ligands, concentrations, temps and times to optimize a reaction. A traditional K factorial experimental scenario can be represented as a set of arms but often doing so isn't natural and there are simpler alternatives to modeling such experimental designs as optimization problems. 

- *Arms* in this context, arm = experiemntal arm and sometimes we use arm as in multi armed bandit problem arm. In the subset selection for experimental design context they are pretty much the same thing.

- We often refer to *maximization* but we mean this in the most generic way. Maximizing F(x) is the same as Mimimizing -1 * F(X) and a < F(X) < b can be phrased as maximization ex) -1 * step_wise_penaltyF(X) or smoothed_penaltyF(X). 

# Questions to ask ourseleves

#### 1. What are the goals of this experiment/library screen?
- Are we trying to find the experiemntal arm that maximizes some observable properties? **Black-Box / Bayesian optimization** 
- Gather data for testing statistical hypothesis of relationships between variables? **Optimal Experimental Design** 
- Building a training data set for a model? **Active Learning**
- collecting data to benchmark a model?
- Drug discovery? Finding hits? targets? optomizing lead series?
- All of the above?

#### 2. What makes one set of N arms better for our experiment goals than another? 
- Say we were interested in hit finding
    - What makes a compound a hit?
    - What makes a set of N compounds more likely to contain hits?
- And if we are interested in learning about a biological system?
    - What does it mean to learn about a biological system?
    - Can we reduce the concept of 'learning' to estimating paramters of a gerealized linear model? or testing statsitical hypothesis?  
    - What can we learn if we test 1 compound? What can we learn if we test 2, 3, …N, N+1
    

#### 3. How to represent experiment goals as [set functions](https://en.wikipedia.org/wiki/Set_function)
- Can we devise a function? or do we only have a black-box oracle? 
- Are functions [monotomic](https://en.wikipedia.org/wiki/Monotonic_function) with respect to [cardinality](https://en.wikipedia.org/wiki/Cardinality)? 
- Are they modular or [submodular](https://en.wikipedia.org/wiki/Submodular_set_function)?
- Can we make proxy functions that are monotomic to our true objectives but easier to optomize?
- What is the wall time of a function? how does it scale with set cardinality?  


## Related Topics
- Optimal Experimental Design
- Bayesian Optimization 
- Active Learning
- Chemogenomics
- Decison Science

## Related Tools / Publications
Treating experimental design and assay data selection as an optimization problem is nothing new. There are already some great tools on the topic. This work focusses on subset selection scenarios that are often encoutered in perturbational assays for hit finding, lead op, target discovery, systems understanding and model development. Below are some great other resources, many of which are dependencies for this repo. 

- [pymoo](https://pymoo.org/index.html)
- [Optimal experimental Design](https://www.nature.com/articles/s41592-018-0083-2)
- [DGEMO](https://github.com/yunshengtian/DGEMO?tab=readme-ov-file)
- [pyDOE2](https://github.com/clicumu/pyDOE2)
- [DoEgen](https://github.com/sebhaan/DoEgen)
- [GPdoemd](https://www.sciencedirect.com/science/article/abs/pii/S0098135419300468)
- [AutoOED](https://autooed.readthedocs.io/en/latest/index.html)
- [Generalized Subset Designs in Analytical Chemistry](https://pubs.acs.org/doi/10.1021/acs.analchem.7b00506)


