# MultiGRN
# Input files: 
## expr_matrix:  expression matrix (row: cell, column: gene)
## cell_adj: Cell-Cell interaction adjacency matrix
## gene_adj: Gene regulatory network adjacency matrix


import numpy as np
import networkx as nx
from collections import defaultdict
import random
import pandas as pd

class MultilayerInfluenceSimulator:
    def __init__(self, expr_matrix, gene_adj, spatial_dist, cell_adj, threshold=0.1, cell_names=None, gene_names=None):
        self.expr = expr_matrix
        self.gene_adj = gene_adj
        self.cell_adj = cell_adj
        self.threshold = threshold
        self.n_cells, self.n_genes = expr_matrix.shape
        self.cell_names = cell_names if cell_names else list(range(self.n_cells))
        self.gene_names = gene_names if gene_names else list(range(self.n_genes))
        self.spatial_dist = spatial_dist
        self.G = self._build_multilayer_graph()
        self.valid_nodes = list(self.G.nodes)

    def _build_multilayer_graph(self):
        G = nx.DiGraph()
        for i in range(self.n_cells):
            expressed_genes = np.where(self.expr[i] > 0)[0]
            for m in expressed_genes:
                for n in expressed_genes:
                    w = self.gene_adj.iloc[m, n]
                    if w != 0:
                        G.add_edge((i, m), (i, n), weight=w)

        for i in range(self.n_cells):
            for j in range(self.n_cells):
                if i == j:
                    continue
                cci = self.cell_adj.iloc[i, j]
                if cci < self.threshold:
                    continue
                dist = self.spatial_dist.iloc[i, j]
                if dist == 0:
                    continue
                for m in range(self.n_genes):
                    if self.expr[i, m] == 0:
                        continue
                    for n in range(self.n_genes):
                        w = self.gene_adj.iloc[m, n]
                        if w == 0:
                            continue
                        influence = dist * cci * w
                        G.add_edge((i, m), (j, n), weight=influence)
        return G

    def one_step_propagation(self, seed_set):
        activated = set(seed_set)
        influence_record = defaultdict(float)
        for i in activated:

            for m in range(self.n_genes):
                expr_val = self.expr[i, m]
                if expr_val == 0:
                    continue
                node = (i, m)
                csname = self.cell_names[i]
                gsname = self.gene_names[m]
                if node not in self.G:
                    continue
                for neighbor in self.G.successors(node):
                    w = self.G[node][neighbor]['weight']
                    j, n = neighbor
                    cname = self.cell_names[j]
                    gname = self.gene_names[n]
                    if csname == cname:
                        continue
                    influence_record[(csname, gsname ,cname, gname)] += w
        return influence_record


    def monte_carlo(self, seed_set=None, T=10, delta=1e-4, check_convergence=True, random_seed=True):
        if seed_set is None:
            if not random_seed:
                raise ValueError("seed_set is null。")
            seed_set = [random.choice(range(self.n_cells))]  
        seed_names = [self.cell_names[i] if isinstance(i, int) else i for i in seed_set]

        print(f"use seed nodes: {seed_names}")

        total_influence = defaultdict(float)
        prev_avg = defaultdict(float)

        for t in range(1, T + 1):
            influence = self.one_step_propagation(seed_set)
            
            for node, val in influence.items():
                total_influence[node] += val

            if check_convergence:
                avg_infl = {k: v / t for k, v in total_influence.items()}
                delta_val = sum(abs(avg_infl[k] - prev_avg.get(k, 0.0)) for k in avg_infl)
                print(f"Iterative {t} ，Δ = {delta_val:.6f}")
                if delta_val < delta:
                    print("✅ stop early")
                    break
                prev_avg = avg_infl

        return {k: v / t for k, v in total_influence.items()}

    def display_top_k(self, result, k=10):
        topk = sorted(result.items(), key=lambda x: -x[1])[:k]
        for (ci, gi), val in topk:
            cname = self.cell_names[ci]
            gname = self.gene_names[gi]
            print(f"Cell: {cname}, Gene: {gname}, Influence Score: {val:.4f}")



class BootstrapInfluenceSimulator:
    def __init__(self, sim, expr_matrix, num_seeds=10, T=100, seed_expr_threshold=0.1):

        self.sim = sim
        self.expr_matrix = expr_matrix
        self.num_seeds = num_seeds
        self.T = T
        self.seed_expr_threshold = seed_expr_threshold
        self.valid_nodes = list(range(self.sim.n_cells))
        
    def run_bootstrap(self, n_rounds=100):

        influence_sum = defaultdict(float)
        all_results = []

        for r in range(n_rounds):
            seed_set = random.choices(self.valid_nodes, k=self.num_seeds)
            result = self.sim.monte_carlo(seed_set=seed_set, T=self.T, check_convergence=True)
            all_results.append(result)

            for k, v in result.items():
                influence_sum[k] += v

        averaged_result = {k: influence_sum[k] / n_rounds for k in influence_sum}
        self.display_top_k(averaged_result, k=10)
        return averaged_result, all_results
        
    def display_top_k(self, result, k=10):
        topk = sorted(result.items(), key=lambda x: -x[1])[:k]
        for (sci, sgi,tci,tgi), val in topk:
            scname = self.sim.cell_names[sci] if isinstance(sci, int) else sci
            sgname = self.sim.gene_names[sgi] if isinstance(sgi, int) else sgi
            tcname = self.sim.cell_names[tci] if isinstance(tci, int) else tci
            tgname = self.sim.gene_names[tgi] if isinstance(tgi, int) else tgi
            print(f"Source Cell: {scname}, Source Gene: {sgname} - Target Cell: {tcname}, Target Gene: {tgname}, Influence Score: {val:.6f}")
            

sim = MultilayerInfluenceSimulator(
  cell_adj=cell_adj,
  gene_adj=gene_adj,
  spatial_dist = spatial_dist,
  expr_matrix=expr_matrix.values,  
  cell_names=expr_matrix.index.tolist(),
  gene_names=expr_matrix.columns.tolist(),
  threshold=0.1
)

bootstrap_sim = BootstrapInfluenceSimulator(
    sim=sim,                         
    expr_matrix=expr_matrix.values,       
    num_seeds=5,                    
    T=50,                          
    seed_expr_threshold=0.1        
)

averaged_result, all_runs = bootstrap_sim.run_bootstrap(n_rounds=100)


