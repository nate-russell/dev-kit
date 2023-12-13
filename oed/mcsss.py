from pymoo.algorithms.soo.nonconvex.ga import GA
from pymoo.operators.crossover.pntx import TwoPointCrossover
from pymoo.operators.mutation.bitflip import BitflipMutation
from pymoo.operators.sampling.rnd import BinaryRandomSampling
from pymoo.optimize import minimize
from pymoo.problems.single.knapsack import create_random_knapsack_problem
from pymoo.core.problem import Problem
from pymoo.algorithms.moo.nsga2 import NSGA2
import numpy as np
from pymoo.core.problem import ElementwiseProblem

from pymoo.core.crossover import Crossover
from pymoo.core.mutation import Mutation
from pymoo.core.sampling import Sampling


class MySampling(Sampling):

    def _do(self, problem, n_samples, **kwargs):
        X = np.full((n_samples, problem.n_var), False, dtype=bool)

        for k in range(n_samples):
            I = np.random.permutation(problem.n_var)[:problem.n_max]
            X[k, I] = True

        return X


class BinaryCrossover(Crossover):
    def __init__(self):
        super().__init__(2, 1)

    def _do(self, problem, X, **kwargs):
        n_parents, n_matings, n_var = X.shape

        _X = np.full((self.n_offsprings, n_matings, problem.n_var), False)

        for k in range(n_matings):
            p1, p2 = X[0, k], X[1, k]

            both_are_true = np.logical_and(p1, p2)
            _X[0, k, both_are_true] = True

            n_remaining = problem.n_max - np.sum(both_are_true)

            I = np.where(np.logical_xor(p1, p2))[0]

            S = I[np.random.permutation(len(I))][:n_remaining]
            _X[0, k, S] = True

        return _X


class MyMutation(Mutation):
    def _do(self, problem, X, **kwargs):
        for i in range(X.shape[0]):
            X[i, :] = X[i, :]
            is_false = np.where(np.logical_not(X[i, :]))[0]
            is_true = np.where(X[i, :])[0]
            X[i, np.random.choice(is_false)] = True
            X[i, np.random.choice(is_true)] = False

        return X

class MultiCriteriaSubSetSolver():

    def __init__(self,objective_functions,ids) -> None:
        self.objective_functions = objective_functions
        self.ids = ids
        pass

    def _build_subset_problem(self,n_max):
        return SubsetProblem(self.objective_functions,self.ids,n_max=n_max)
        

    def solve(self,n_max):
        ssp = self._build_subset_problem(n_max)

        algorithm = NSGA2(
            pop_size=200,
            sampling=MySampling(),
            crossover=BinaryCrossover(),
            mutation=MyMutation(),
            eliminate_duplicates=True
            )

        res = minimize(ssp,
                       algorithm,
                       ('n_gen', 400),
                       seed=1,
                       verbose=True)
        
        return res

    def presolve(self):
        pass


class SubsetProblem(ElementwiseProblem):
    def __init__(self,
                 objective_functions,
                 ids,
                 n_max,
                 ):
        super().__init__(n_var=len(ids), n_obj=len(objective_functions), n_ieq_constr=1)
        self.objective_functions = objective_functions
        self.n_max = n_max

    def _evaluate(self, x, out, *args, **kwargs):
        obj_funcs = []
        self.multipliers = []
        for maxmin, label, f in self.objective_functions:
            #print(maxmin,label)
            #print(f)
            #print(type(x))
            #print(sum(x))
            #print(x)
            m = None
            if maxmin.lower() == 'max':
                m = -1
            elif maxmin.lower() == 'min':
                m  = 1
            else:
                raise ValueError(f"{repr(maxmin)} is unexpected. Should be 'max' or 'min'")
            
            self.multipliers.append(m)
            obj_funcs.append(m*f(x))

        out["F"] = obj_funcs
        out["G"] = (self.n_max - np.sum(x)) ** 2





