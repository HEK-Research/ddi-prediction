# node2vec environment notes (`ddi_stable`)

Status as of 2026-09-02 (updated -- previously `pecanpy` was missing, see history below):

- **`pecanpy` (2.0.9):** installed and working. This is now the node2vec backend used by [`src/node2vec_utils.py`](../src/node2vec_utils.py) -- it simulates the random walks (`SparseOTF`, on-the-fly transition probabilities, no precompute step) and trains embeddings via `pecanpy`'s built-in `embed()`, which calls gensim's `Word2Vec` internally using gensim's current (>= 4.0) kwarg names natively.
- **`gensim` (4.4.0):** installed and working -- no compatibility shim needed anymore, since `pecanpy` already targets modern gensim.
- **`node2vec` (PyPI package, 0.3.0):** still installed but **no longer used**. It's unmaintained and calls `gensim.models.Word2Vec(..., size=..., iter=...)`, using kwarg names gensim renamed to `vector_size` / `epochs` in its 4.0 release -- this raised `TypeError: Word2Vec.__init__() got an unexpected keyword argument 'size'` when called directly against gensim 4.4.0. `src/node2vec_utils.py` previously worked around this with a monkeypatch shim; that approach has been replaced by the `pecanpy`-based implementation now that `pecanpy` installs cleanly in this environment.

**Resolved:** the three `node2vec_embedding.ipynb` notebooks (in `h1_biological_overlap/`, `h2_pharmacological_similarity/`, `h3_structural_similarity/`) now run node2vec via `pecanpy` + gensim through `src/node2vec_utils.py`, with no known compatibility issues.
