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
from pymoo.util.ref_dirs import get_reference_directions
from pymoo.visualization.scatter import Scatter
from tqdm.auto import tqdm as base_tqdm
# or: from tqdm.auto import tqdm as base_tqdm

class tqdm(base_tqdm):
    def update(self, n=1):
        super(tqdm, self).update(1)
        

def selected_to_x(selected,nvar):
            x = np.array([False for i in range(nvar)])
            for i in selected:
                x[i] = True
            return x
    
class GreedySelection(Sampling):

    def _do(self, problem, n_samples, **kwargs):
        
        objective_functions = problem.objective_functions
        n_max = problem.n_max


        population = problem.pop_size
        n_objectives = len(objective_functions)
        #print('Riesz energy reference directions')
        ref_dirs = get_reference_directions("energy", n_objectives, population, seed=1)
        #print(ref_dirs)
        #print(ref_dirs.shape)

        F = []
        M = []
        for j in range(n_objectives):
            maxmin, label, f = objective_functions[j]
            m = None
            if maxmin.lower() == 'max':
                m = -1
            elif maxmin.lower() == 'min':
                m  = 1
            else:
                raise ValueError(f"{repr(maxmin)} is unexpected. Should be 'max' or 'min'")
            M.append(m)
            F.append(f)

        obj_funcs = []
        #print('Same X, different F')
        const_x = selected_to_x(set([1,2,3,4,5]),nvar=problem.n_var)
        for i in range(population):
            #print('pop',i)
            W = []
            for j in range(n_objectives):
                w = ref_dirs[i,j]
                #print(f"\tObj {j+1}: {w}")
                W.append(w)
            
            #wmf = lambda x: np.sum(m*f(x)*w for m,w,f in zip(M,W,F))
            def funcC(W,M,F):
                def func(x):
                    return np.sum((m*f(x)*w for m,w,f in zip(M,W,F)))
                return func

            f = funcC(W,M,F)
            #print(W,f(const_x))
            obj_funcs.append(funcC(W,M,F))


        #print('Same X, different F')
        #for f in obj_funcs:
            #print(f(const_x),f)
        

        all_x = []
        for f in tqdm(obj_funcs,desc='Greedy Selection Init'):
            
            selected = set()
            not_selected = set(range(problem.n_var))
            
            for i in range(n_max):
                i = i+1
                #print(f"Selected: ({len(selected)}) {selected}")
                
                fvals = []
                j_list = []
                for j in not_selected:
                    tmp_selected = set([j])
                    tmp_selected = tmp_selected.union(selected)
                    x = selected_to_x(tmp_selected,nvar=problem.n_var)
                    fval = f(x)
                    fvals.append(fval)
                    j_list.append(j)
                
                #print('min/max/mean = ',min(fvals),max(fvals),np.mean(fvals),len(fvals))
                minj_index = np.argmin(fvals)
                minj = j_list[minj_index]
                selected.add(minj)
                not_selected.remove(minj)
            
            x = selected_to_x(selected,nvar=problem.n_var)
            all_x.append(x)
        X = np.vstack(all_x)
        #print('Greedy X')
        #print(X.shape)
        #print(X)
        import matplotlib.pyplot as plt
        plt.imshow(X)
        plt.show()

        '''
        X = np.full((n_samples, problem.n_var), False, dtype=bool)
        for k in range(n_samples):
            I = np.random.permutation(problem.n_var)[:problem.n_max]
            X[k, I] = True
        '''
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

    def _build_subset_problem(self,n_max,pop_size):
        return SubsetProblem(self.objective_functions,self.ids,n_max=n_max,pop_size=pop_size)
        

    def solve(self,n_max,pop_size=100,n_gens=200):
        ssp = self._build_subset_problem(n_max,pop_size)

        algorithm = NSGA2(
            pop_size=pop_size,
            sampling=GreedySelection(),
            crossover=BinaryCrossover(),
            mutation=MyMutation(),
            eliminate_duplicates=True,
            save_history=True,
            )

        with tqdm(desc='NSGA-II',total=n_gens,unit='generation', unit_scale=True) as t:
            res = minimize(
                    ssp,
                    algorithm,
                    ('n_gen', n_gens),
                    seed=1,
                    callback=t.update,
                    verbose=False
                    )
        
        return res

    def presolve(self):
        pass


class SubsetProblem(ElementwiseProblem):
    def __init__(self,
                 objective_functions,
                 ids,
                 n_max,
                 pop_size,
                 ):
        super().__init__(n_var=len(ids), n_obj=len(objective_functions), n_ieq_constr=1)
        self.objective_functions = objective_functions
        self.n_max = n_max
        self.pop_size = pop_size

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





