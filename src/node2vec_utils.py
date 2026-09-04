"""Shared node2vec runner for the H1/H2/H3 similarity graphs.

Uses `pecanpy` (the `SparseOTF` "on the fly" transition-probability variant)
to simulate the biased random walks, then gensim's `Word2Vec` (via pecanpy's
`embed()`) to train the skip-gram embeddings. pecanpy is actively maintained
and targets modern gensim (>= 4.0) natively -- unlike the older PyPI
`node2vec` package, which needed a compatibility shim to work with gensim's
renamed `size`/`iter` kwargs. See `notebooks/NODE2VEC_ENV_NOTES.md` for that
superseded approach and the rest of the environment history.
"""

from __future__ import annotations

from pathlib import Path
from typing import Union

import networkx as nx
import pandas as pd
from pecanpy.pecanpy import SparseOTF


def run_node2vec(
    G: nx.Graph,
    dimensions: int = 64,
    walk_length: int = 30,
    num_walks: int = 10,
    p: float = 1.0,
    q: float = 1.0,
    weight_key: str = "weight",
    workers: int = 1,
    window: int = 10,
    epochs: int = 5,
    seed: int = 42,
) -> pd.DataFrame:
    """Run node2vec (via pecanpy) on `G`, returning one embedding row per node, indexed by node id."""
    node_ids = list(G.nodes())
    adj_mat = nx.to_numpy_array(G, nodelist=node_ids, weight=weight_key)

    graph = SparseOTF.from_mat(adj_mat, node_ids, p=p, q=q, workers=workers, random_state=seed)
    vectors = graph.embed(
        dim=dimensions,
        num_walks=num_walks,
        walk_length=walk_length,
        window_size=window,
        epochs=epochs,
    )

    embeddings = pd.DataFrame(vectors, index=node_ids, columns=[f"dim_{i}" for i in range(dimensions)])
    embeddings.index.name = "drug_id"
    return embeddings


def embed_graphml(path: Union[str, Path], **node2vec_kwargs) -> pd.DataFrame:
    """Load a GraphML file and run node2vec on it."""
    G = nx.read_graphml(path)
    return run_node2vec(G, **node2vec_kwargs)


def save_embeddings(embeddings: pd.DataFrame, out_path: Union[str, Path]) -> None:
    """Save an embedding DataFrame (as returned by `run_node2vec`) to parquet."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    embeddings.to_parquet(out_path)
