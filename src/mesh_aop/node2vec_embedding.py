# -*- coding: utf-8 -*-
"""
node2vec_embedding.py - self-contained Node2Vec graph embedding (vendored).

A dependency-light reimplementation of the walk-generation and Word2Vec embedding
from the `node2vec` package (Elior Cohen, MIT-licensed), used by the Node2Vec
dendrogram figure. It reproduces the exact algorithm - biased second-order random
walks (return parameter `p`, in-out parameter `q`) fed to gensim's skip-gram
Word2Vec - but drops the `node2vec` package so the pipeline is free of its
transitive import constraints (its `edges.py` imports `pkg_resources`, removed in
setuptools 81). This module needs only gensim, networkx, numpy and the standard
library, and therefore works with modern numpy / networkx / gensim / Python.

Calculations are identical to the upstream package; only the packaging changes.
Walk generation is run serially here (upstream splits it across joblib workers
purely for speed - the algorithm and its statistics are unchanged), and the
gensim version check now reads `gensim.__version__` instead of `pkg_resources`.
"""

import random
from collections import defaultdict

import gensim
import networkx as nx
import numpy as np

try:
    from tqdm.auto import tqdm
except Exception:  # tqdm is optional; degrade to a no-op iterator wrapper
    def tqdm(iterable=None, **kwargs):
        return iterable if iterable is not None else []

_FIRST_TRAVEL_KEY = 'first_travel_key'
_PROBABILITIES_KEY = 'probabilities'
_NEIGHBORS_KEY = 'neighbors'


def _edge_weight(graph, u, v, weight_key):
    """Return the edge weight exactly as upstream node2vec does (1 if absent)."""
    try:
        if graph[u][v].get(weight_key):
            return graph[u][v].get(weight_key, 1)
        # MultiGraph fallback: AtlasView keyed by edge id, e.g. {0: {'weight': 0.1}}
        edge = list(graph[u][v])[-1]
        return graph[u][v][edge].get(weight_key, 1)
    except Exception:
        return 1


def _precompute_probabilities(graph, p, q, weight_key, quiet):
    """Precompute biased transition probabilities for each (previous, current) edge.

    Mirrors node2vec's `_precompute_probabilities`: for each source node and each
    of its neighbours `current_node`, store the normalized probability of stepping
    to each of `current_node`'s neighbours given that we arrived from `source`
    (weighting return moves by 1/p and out moves by 1/q), plus the first-step
    weights and neighbour list for `source`.
    """
    d_graph = defaultdict(dict)
    nodes_iter = graph.nodes() if quiet else tqdm(graph.nodes(), desc='Computing transition probabilities')

    for source in nodes_iter:
        if _PROBABILITIES_KEY not in d_graph[source]:
            d_graph[source][_PROBABILITIES_KEY] = dict()

        for current_node in graph.neighbors(source):
            if _PROBABILITIES_KEY not in d_graph[current_node]:
                d_graph[current_node][_PROBABILITIES_KEY] = dict()

            unnormalized_weights = list()

            for destination in graph.neighbors(current_node):
                weight = _edge_weight(graph, current_node, destination, weight_key)

                if destination == source:               # backwards (return) move
                    ss_weight = weight * 1 / p
                elif destination in graph[source]:       # destination also adjacent to source
                    ss_weight = weight
                else:                                    # out move
                    ss_weight = weight * 1 / q

                unnormalized_weights.append(ss_weight)

            unnormalized_weights = np.array(unnormalized_weights)
            d_graph[current_node][_PROBABILITIES_KEY][source] = unnormalized_weights / unnormalized_weights.sum()

        first_travel_weights = []
        for destination in graph.neighbors(source):
            first_travel_weights.append(graph[source][destination].get(weight_key, 1))
        first_travel_weights = np.array(first_travel_weights)
        d_graph[source][_FIRST_TRAVEL_KEY] = first_travel_weights / first_travel_weights.sum()

        d_graph[source][_NEIGHBORS_KEY] = list(graph.neighbors(source))

    return d_graph


def _generate_walks(d_graph, walk_length, num_walks, quiet):
    """Generate `num_walks` biased random walks of length `walk_length` from every node.

    Identical to node2vec's `parallel_generate_walks` (run serially): for each walk
    pass, shuffle the nodes and walk from each, choosing the next node with the
    precomputed first-travel or (previous, current) transition probabilities.
    """
    walks = list()
    iterator = range(num_walks) if quiet else tqdm(range(num_walks), desc='Generating walks')

    for _ in iterator:
        shuffled_nodes = list(d_graph.keys())
        random.shuffle(shuffled_nodes)

        for source in shuffled_nodes:
            walk = [source]

            while len(walk) < walk_length:
                walk_options = d_graph[walk[-1]].get(_NEIGHBORS_KEY, None)
                if not walk_options:  # dead-end node
                    break

                if len(walk) == 1:
                    probabilities = d_graph[walk[-1]][_FIRST_TRAVEL_KEY]
                else:
                    probabilities = d_graph[walk[-1]][_PROBABILITIES_KEY][walk[-2]]

                walk_to = random.choices(walk_options, weights=probabilities)[0]
                walk.append(walk_to)

            walks.append(list(map(str, walk)))  # tokens must be strings for Word2Vec

    return walks


class Node2Vec:
    """Drop-in replacement for `node2vec.Node2Vec` covering the features the pipeline uses.

    Precomputes biased transition probabilities and generates the random walks on
    construction; `fit(**word2vec_params)` returns a trained gensim Word2Vec model
    whose `.wv.index_to_key` / `.wv.vectors` feed the dendrogram. `workers` is
    forwarded to gensim for training parallelism (walk generation is serial).
    """
    def __init__(self, graph: nx.Graph, dimensions: int = 128, walk_length: int = 80, num_walks: int = 10,
                 p: float = 1, q: float = 1, weight_key: str = 'weight', workers: int = 1,
                 quiet: bool = False, seed: int = None):
        self.dimensions = dimensions
        self.workers = workers

        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)

        d_graph = _precompute_probabilities(graph, p, q, weight_key, quiet)
        self.walks = _generate_walks(d_graph, walk_length, num_walks, quiet)

    def fit(self, **skip_gram_params) -> gensim.models.Word2Vec:
        """Train skip-gram Word2Vec on the generated walks and return the gensim model."""
        if 'workers' not in skip_gram_params:
            skip_gram_params['workers'] = self.workers

        # gensim renamed 'size' -> 'vector_size' in v4.0.0.
        try:
            gensim_major = int(gensim.__version__.split('.')[0])
        except Exception:
            gensim_major = 4
        size_key = 'vector_size' if gensim_major >= 4 else 'size'
        if 'vector_size' not in skip_gram_params and 'size' not in skip_gram_params:
            skip_gram_params[size_key] = self.dimensions

        if 'sg' not in skip_gram_params:
            skip_gram_params['sg'] = 1

        return gensim.models.Word2Vec(self.walks, **skip_gram_params)
