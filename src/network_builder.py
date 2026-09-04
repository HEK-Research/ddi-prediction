"""Network construction and diagnostics for DDI similarity graphs.

Turns a dense pairwise similarity matrix (identity/GO/pathway Tversky or Jaccard
scores, ATC similarity, Morgan/Tanimoto similarity, etc.) into a sparsified
``networkx`` graph, and validates that the resulting topology is a sane input
for node2vec: edge weights should reflect biology rather than profile-size
artifacts, degree distribution should be reasonable, and graph topology
should actually separate known DDI pairs from random pairs.

Typical usage::

    sim_df = pd.DataFrame(sim_matrix, index=drug_ids, columns=drug_ids)
    G = build_sparsified_graph(sim_df, method="hybrid", k=15, min_sim=0.05)
    diag = NetworkDiagnostics(G, size_vector=size_vector, ddi_edges=ddi_edges, sim_df=sim_df)
    report = diag.run_all(contained_pair=("DB001", "DB002"), unrelated_pair=("DB001", "DB999"))
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional, Sequence, Set, Tuple, Union

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.metrics import roc_auc_score

DrugId = str
Edge = Tuple[DrugId, DrugId]

# ============================================================================
# Shared helpers
# ============================================================================


def derive_size_vector(profile_sets: Mapping[DrugId, Set[Any]]) -> Dict[DrugId, int]:
    """Derive a size vector (profile cardinality per drug) from ``profile_sets``."""
    return {drug_id: len(features) for drug_id, features in profile_sets.items()}


def _as_labeled_frame(
    sim_matrix: Union[np.ndarray, pd.DataFrame],
    node_ids: Optional[Sequence[DrugId]] = None,
) -> pd.DataFrame:
    """Normalize a raw ndarray or DataFrame into a square, drug-id-labeled DataFrame."""
    if isinstance(sim_matrix, pd.DataFrame):
        sim_df = sim_matrix
    else:
        if node_ids is None:
            raise ValueError("`node_ids` must be provided when `sim_matrix` is a plain numpy array.")
        sim_df = pd.DataFrame(sim_matrix, index=node_ids, columns=node_ids)

    if list(sim_df.index) != list(sim_df.columns):
        raise ValueError("`sim_matrix` index and columns must match (same drug ID order).")
    if sim_df.shape[0] != sim_df.shape[1]:
        raise ValueError("`sim_matrix` must be square.")
    return sim_df


def _upper_triangle_values(sim_df: pd.DataFrame) -> np.ndarray:
    """Flat array of the strictly-upper-triangular (i.e. one-per-pair) similarity values."""
    arr = sim_df.to_numpy()
    iu = np.triu_indices_from(arr, k=1)
    return arr[iu]


def _safe_spearmanr(x: Sequence[float], y: Sequence[float]) -> float:
    """Spearman rho, returning NaN (instead of scipy's ConstantInputWarning) when either side has no variance."""
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    if x_arr.size < 2 or np.all(x_arr == x_arr[0]) or np.all(y_arr == y_arr[0]):
        return float("nan")
    rho, _ = spearmanr(x_arr, y_arr)
    return float(rho)


def edge_list_to_matrix(
    df: pd.DataFrame,
    weight_col: str,
    drug1_col: str = "drug1_id",
    drug2_col: str = "drug2_id",
) -> pd.DataFrame:
    """Pivot a long (drug1_id, drug2_id, weight) pair table into a dense, symmetric matrix.

    Intended for datasets that only contain a *sample* of pairs (e.g. curated adverse or
    non-interacting DDI pairs) rather than the full N-choose-2 universe of drug pairs -- any
    pair absent from `df` is filled with 0.0, so "not sampled" and "observed as dissimilar"
    are both encoded the same way. Duplicate (a, b) rows are collapsed by averaging `weight_col`.
    """
    pairs = df[[drug1_col, drug2_col, weight_col]].copy()
    pairs = pairs.groupby([drug1_col, drug2_col], as_index=False)[weight_col].mean()

    drug_ids = pd.Index(sorted(set(pairs[drug1_col]) | set(pairs[drug2_col])), name="drug_id")
    # Build on a plain numpy array (not `DataFrame.values`) -- pandas' copy-on-write mode can
    # hand back a read-only array from `.values`, which would make in-place assignment fail.
    matrix = np.zeros((len(drug_ids), len(drug_ids)), dtype=float)

    i = drug_ids.get_indexer(pairs[drug1_col])
    j = drug_ids.get_indexer(pairs[drug2_col])
    w = pairs[weight_col].to_numpy(dtype=float)
    matrix[i, j] = w
    matrix[j, i] = w
    np.fill_diagonal(matrix, 0.0)
    return pd.DataFrame(matrix, index=drug_ids, columns=drug_ids)


# ============================================================================
# Module 1: Network sparsification
# ============================================================================


def find_knee_threshold(sim_df: pd.DataFrame, n_steps: int = 200) -> float:
    """Find the knee of the edges-kept-vs-threshold curve (max distance to the endpoint chord)."""
    values = _upper_triangle_values(sim_df)
    values = values[~np.isnan(values)]
    if values.size == 0:
        return 0.0

    thresholds = np.linspace(values.min(), values.max(), n_steps)
    edge_counts = np.array([(values > t).sum() for t in thresholds], dtype=float)

    x = (thresholds - thresholds.min()) / (thresholds.max() - thresholds.min() + 1e-12)
    y = (edge_counts - edge_counts.min()) / (edge_counts.max() - edge_counts.min() + 1e-12)

    x1, y1, x2, y2 = x[0], y[0], x[-1], y[-1]
    numerator = np.abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1)
    denominator = np.sqrt((y2 - y1) ** 2 + (x2 - x1) ** 2) + 1e-12
    distances = numerator / denominator

    knee_idx = int(np.argmax(distances))
    return float(thresholds[knee_idx])


def sparsify_global(
    sim_df: pd.DataFrame,
    threshold: Optional[float] = None,
    percentile: Optional[float] = None,
    use_knee: bool = False,
) -> nx.Graph:
    """Keep edges with weight > threshold. Exactly one of the three threshold sources must be given."""
    n_sources = sum(x is not None for x in (threshold, percentile)) + int(use_knee)
    if n_sources != 1:
        raise ValueError("Specify exactly one of `threshold`, `percentile`, or `use_knee=True`.")

    if use_knee:
        resolved_threshold = find_knee_threshold(sim_df)
    elif percentile is not None:
        resolved_threshold = float(np.percentile(_upper_triangle_values(sim_df), percentile))
    else:
        resolved_threshold = float(threshold)  # type: ignore[arg-type]

    G = nx.Graph()
    G.add_nodes_from(sim_df.index)

    arr = sim_df.to_numpy()
    ids = sim_df.index.to_numpy()
    iu = np.triu_indices_from(arr, k=1)
    row_idx, col_idx = iu
    mask = arr[row_idx, col_idx] > resolved_threshold
    for a, b, w in zip(ids[row_idx[mask]], ids[col_idx[mask]], arr[row_idx[mask], col_idx[mask]]):
        G.add_edge(a, b, weight=float(w))

    G.graph["threshold_used"] = resolved_threshold
    return G


def sparsify_knn(sim_df: pd.DataFrame, k: int) -> nx.Graph:
    """Keep each node's k highest-similarity neighbors, symmetrized (union) into an undirected graph."""
    G = nx.Graph()
    ids = sim_df.index.to_numpy()
    G.add_nodes_from(ids)

    arr = sim_df.to_numpy()
    n = len(ids)
    k_eff = min(k, n - 1)
    if k_eff <= 0:
        return G

    for i in range(n):
        row = arr[i].copy()
        row[i] = -np.inf  # a node is never its own neighbor
        neighbor_idx = np.argpartition(row, -k_eff)[-k_eff:]
        for j in neighbor_idx:
            # symmetric sim_matrix => weight is identical regardless of which side proposed the edge
            G.add_edge(ids[i], ids[j], weight=float(arr[i, j]))

    return G


def sparsify_hybrid(sim_df: pd.DataFrame, k: int, min_sim: float) -> nx.Graph:
    """Top-k neighbors per node, but drop any candidate edge below `min_sim` (up to k, never below the bar)."""
    G = nx.Graph()
    ids = sim_df.index.to_numpy()
    G.add_nodes_from(ids)

    arr = sim_df.to_numpy()
    n = len(ids)
    k_eff = min(k, n - 1)
    if k_eff <= 0:
        return G

    for i in range(n):
        row = arr[i].copy()
        row[i] = -np.inf
        neighbor_idx = np.argpartition(row, -k_eff)[-k_eff:]
        for j in neighbor_idx:
            weight = arr[i, j]
            if weight >= min_sim:
                G.add_edge(ids[i], ids[j], weight=float(weight))

    return G


def build_sparsified_graph(
    sim_matrix: Union[np.ndarray, pd.DataFrame],
    method: str = "global",
    node_ids: Optional[Sequence[DrugId]] = None,
    threshold: Optional[float] = None,
    percentile: Optional[float] = None,
    use_knee: bool = False,
    k: Optional[int] = None,
    min_sim: Optional[float] = None,
) -> nx.Graph:
    """Dispatch to `sparsify_global` / `sparsify_knn` / `sparsify_hybrid` based on `method`."""
    sim_df = _as_labeled_frame(sim_matrix, node_ids)

    if method == "global":
        return sparsify_global(sim_df, threshold=threshold, percentile=percentile, use_knee=use_knee)
    if method == "knn":
        if k is None:
            raise ValueError("`k` is required for method='knn'.")
        return sparsify_knn(sim_df, k=k)
    if method == "hybrid":
        if k is None or min_sim is None:
            raise ValueError("`k` and `min_sim` are both required for method='hybrid'.")
        return sparsify_hybrid(sim_df, k=k, min_sim=min_sim)
    raise ValueError(f"Unknown method '{method}'. Expected one of: 'global', 'knn', 'hybrid'.")


# ============================================================================
# Module 2: Network diagnostics
# ============================================================================


@dataclass
class NetworkDiagnostics:
    """Diagnostic suite validating that a sparsified similarity graph is fit for node2vec."""

    graph: nx.Graph
    size_vector: Mapping[DrugId, float]
    ddi_edges: Sequence[Edge] = field(default_factory=list)
    sim_df: Optional[pd.DataFrame] = None  # required for the sanity probe, enrichment, and AUROC

    # --- Diagnostic A: meaningful edge weights -----------------------------

    def size_similarity_correlation(self) -> Dict[str, float]:
        """Spearman rho between edge weight and max size / size diff / size ratio. Pass signal: rho ~ 0."""
        edges = list(self.graph.edges(data="weight"))
        if not edges:
            return {"rho_max_size": float("nan"), "rho_size_diff": float("nan"), "rho_size_ratio": float("nan")}

        weights: List[float] = []
        max_sizes: List[float] = []
        size_diffs: List[float] = []
        size_ratios: List[float] = []
        for a, b, w in edges:
            size_a, size_b = self.size_vector[a], self.size_vector[b]
            smaller, larger = min(size_a, size_b), max(size_a, size_b)
            weights.append(w)
            max_sizes.append(larger)
            size_diffs.append(larger - smaller)
            size_ratios.append(smaller / larger if larger > 0 else np.nan)

        rho_max_size = _safe_spearmanr(weights, max_sizes)
        rho_size_diff = _safe_spearmanr(weights, size_diffs)
        rho_size_ratio = _safe_spearmanr(weights, size_ratios)
        return {"rho_max_size": rho_max_size, "rho_size_diff": rho_size_diff, "rho_size_ratio": rho_size_ratio}

    def known_pair_sanity_probe(self, contained_pair: Edge, unrelated_pair: Edge) -> Dict[str, float]:
        """Similarity scores for a manually chosen contained pair and unrelated pair."""
        if self.sim_df is None:
            raise ValueError("`sim_df` is required for the known-pair sanity probe.")
        return {
            "contained_pair_similarity": float(self.sim_df.loc[contained_pair]),
            "unrelated_pair_similarity": float(self.sim_df.loc[unrelated_pair]),
        }

    def plot_weight_distribution(self, ax: Optional[plt.Axes] = None, bins: int = 30) -> plt.Axes:
        """Histogram of the graph's non-zero edge weights."""
        weights = [w for _, _, w in self.graph.edges(data="weight")]
        if ax is None:
            _, ax = plt.subplots(figsize=(6, 4))
        sns.histplot(weights, bins=bins, ax=ax)
        ax.set_xlabel("Edge weight (similarity)")
        ax.set_ylabel("Count")
        ax.set_title("Distribution of non-zero edge weights")
        return ax

    # --- Diagnostic B: reasonable degree distribution ----------------------

    def degree_diagnostics(self) -> Dict[str, float]:
        """Isolate count and 5th-percentile degree. Pass signal: 0 isolates, p5 degree >= 5."""
        degrees = np.array([d for _, d in self.graph.degree()], dtype=float)
        if degrees.size == 0:
            return {"isolate_count": 0, "p5_degree": float("nan"), "median_degree": float("nan")}
        return {
            "isolate_count": int((degrees == 0).sum()),
            "p5_degree": float(np.percentile(degrees, 5)),
            "median_degree": float(np.median(degrees)),
        }

    def degree_vs_size(self, plot: bool = True, ax: Optional[plt.Axes] = None) -> Dict[str, float]:
        """Correlation + linear slope between node degree and profile size. Pass signal: flat slope."""
        nodes = list(self.graph.nodes())
        degrees = np.array([self.graph.degree(n) for n in nodes], dtype=float)
        sizes = np.array([self.size_vector[n] for n in nodes], dtype=float)

        rho = _safe_spearmanr(degrees, sizes) if len(nodes) > 1 else float("nan")
        slope = float(np.polyfit(sizes, degrees, 1)[0]) if len(set(sizes)) > 1 else float("nan")

        if plot:
            if ax is None:
                _, ax = plt.subplots(figsize=(6, 4))
            sns.regplot(x=sizes, y=degrees, ax=ax, scatter_kws={"alpha": 0.4})
            ax.set_xlabel("Profile size")
            ax.set_ylabel("Node degree")
            ax.set_title("Degree vs. profile size")

        return {"spearman_rho": rho, "slope": slope}

    # --- Diagnostic C: topology reflects DDI --------------------------------

    def _sample_negative_pairs(self, n: int, random_state: int) -> List[Edge]:
        """Sample `n` random drug pairs that are not in `ddi_edges`."""
        if self.sim_df is None:
            raise ValueError("`sim_df` is required to sample negative pairs.")
        rng = np.random.default_rng(random_state)
        ids = self.sim_df.index.to_numpy()
        known = {frozenset(e) for e in self.ddi_edges}

        negatives: List[Edge] = []
        seen: Set[frozenset] = set()
        max_attempts = max(n * 50, 1000)
        attempts = 0
        while len(negatives) < n and attempts < max_attempts:
            a, b = rng.choice(ids, size=2, replace=False)
            pair_key = frozenset((a, b))
            if pair_key not in known and pair_key not in seen:
                negatives.append((a, b))
                seen.add(pair_key)
            attempts += 1
        return negatives

    def positive_pair_enrichment(self, n_random: Optional[int] = None, random_state: int = 42) -> Dict[str, float]:
        """Compare similarity scores of known DDI pairs vs. a random sample of non-interacting pairs."""
        if self.sim_df is None:
            raise ValueError("`sim_df` is required for enrichment analysis.")
        positive_pairs = [e for e in self.ddi_edges if e[0] in self.sim_df.index and e[1] in self.sim_df.index]
        positive_scores = np.array([self.sim_df.loc[a, b] for a, b in positive_pairs])

        n_neg = n_random or len(positive_scores)
        negative_pairs = self._sample_negative_pairs(n_neg, random_state)
        negative_scores = np.array([self.sim_df.loc[a, b] for a, b in negative_pairs])

        stat, p_value = mannwhitneyu(positive_scores, negative_scores, alternative="greater")
        return {
            "positive_mean": float(positive_scores.mean()),
            "negative_mean": float(negative_scores.mean()),
            "positive_median": float(np.median(positive_scores)),
            "negative_median": float(np.median(negative_scores)),
            "mannwhitney_u": float(stat),
            "p_value": float(p_value),
        }

    def auroc_ddi_prediction(self, n_random: Optional[int] = None, random_state: int = 42) -> float:
        """AUROC using similarity score as the predictor for known-DDI (1) vs. random-pair (0). Pass: > 0.5."""
        if self.sim_df is None:
            raise ValueError("`sim_df` is required for AUROC calculation.")
        positive_pairs = [e for e in self.ddi_edges if e[0] in self.sim_df.index and e[1] in self.sim_df.index]

        n_neg = n_random or len(positive_pairs)
        negative_pairs = self._sample_negative_pairs(n_neg, random_state)

        scores = [self.sim_df.loc[a, b] for a, b in positive_pairs] + [self.sim_df.loc[a, b] for a, b in negative_pairs]
        labels = [1] * len(positive_pairs) + [0] * len(negative_pairs)
        return float(roc_auc_score(labels, scores))

    # --- Convenience wrapper -------------------------------------------------

    def run_all(
        self,
        contained_pair: Optional[Edge] = None,
        unrelated_pair: Optional[Edge] = None,
        n_random: Optional[int] = None,
        random_state: int = 42,
    ) -> Dict[str, Any]:
        """Run every diagnostic that has the data it needs and return a single flat results dict."""
        results: Dict[str, Any] = {}
        results.update(self.degree_diagnostics())

        if self.size_vector:
            results.update({f"size_{k}": v for k, v in self.size_similarity_correlation().items()})
            results.update({f"degree_vs_size_{k}": v for k, v in self.degree_vs_size(plot=False).items()})

        if contained_pair is not None and unrelated_pair is not None:
            results.update(self.known_pair_sanity_probe(contained_pair, unrelated_pair))

        if self.sim_df is not None and len(self.ddi_edges) > 0:
            results.update(self.positive_pair_enrichment(n_random=n_random, random_state=random_state))
            results["auroc_ddi_prediction"] = self.auroc_ddi_prediction(n_random=n_random, random_state=random_state)

        return results


# ============================================================================
# Module 3: Sweep & report
# ============================================================================


@dataclass
class LevelData:
    """Everything `run_diagnostic_sweep` needs for one hierarchical level (e.g. `go_bp`, `pathway_native`)."""

    sim_matrix: Union[np.ndarray, pd.DataFrame]
    size_vector: Mapping[DrugId, float]
    ddi_edges: Sequence[Edge] = field(default_factory=list)
    node_ids: Optional[Sequence[DrugId]] = None


def run_diagnostic_sweep(
    levels: Mapping[str, LevelData],
    param_grid: Sequence[Mapping[str, Any]],
    contained_pair: Optional[Edge] = None,
    unrelated_pair: Optional[Edge] = None,
    n_random: Optional[int] = None,
    random_state: int = 42,
) -> pd.DataFrame:
    """Sparsify + diagnose every (level, param combo) pair and summarize results in one DataFrame.

    `param_grid` entries are kwargs for `build_sparsified_graph`, e.g.
    `{"method": "global", "percentile": 90}`, `{"method": "knn", "k": 15}`, or
    `{"method": "hybrid", "k": 15, "min_sim": 0.05}`.
    """
    rows: List[Dict[str, Any]] = []

    for level_name, level in levels.items():
        sim_df = _as_labeled_frame(level.sim_matrix, level.node_ids)

        for params in param_grid:
            method = params["method"]
            graph_kwargs = {k: v for k, v in params.items() if k != "method"}
            G = build_sparsified_graph(sim_df, method=method, **graph_kwargs)

            diagnostics = NetworkDiagnostics(
                graph=G, size_vector=level.size_vector, ddi_edges=level.ddi_edges, sim_df=sim_df
            )
            result = diagnostics.run_all(
                contained_pair=contained_pair,
                unrelated_pair=unrelated_pair,
                n_random=n_random,
                random_state=random_state,
            )

            threshold_or_k = graph_kwargs.get("k", graph_kwargs.get("threshold", graph_kwargs.get("percentile")))
            rows.append(
                {
                    "level": level_name,
                    "metric": method,
                    "threshold_or_k": threshold_or_k,
                    "spearman_rho_size_disparity": result.get("size_rho_size_diff"),
                    "isolate_count": result.get("isolate_count"),
                    "degree_vs_size_slope": result.get("degree_vs_size_slope"),
                    "auroc_ddi_prediction": result.get("auroc_ddi_prediction"),
                }
            )

    return pd.DataFrame(rows)


if __name__ == "__main__":
    # Small synthetic smoke test -- not real DDI data, just confirms the pipeline runs end to end.
    rng = np.random.default_rng(0)
    n_drugs = 60
    drug_ids = [f"DB{i:04d}" for i in range(n_drugs)]

    profile_sets = {d: set(rng.integers(0, 200, size=rng.integers(3, 30))) for d in drug_ids}
    size_vector = derive_size_vector(profile_sets)

    raw = rng.random((n_drugs, n_drugs))
    sym = (raw + raw.T) / 2
    np.fill_diagonal(sym, 0.0)
    sim_df = pd.DataFrame(sym, index=drug_ids, columns=drug_ids)

    ddi_edges = [(drug_ids[i], drug_ids[i + 1]) for i in range(0, 20, 2)]

    levels = {"demo_level": LevelData(sim_matrix=sim_df, size_vector=size_vector, ddi_edges=ddi_edges)}
    param_grid = [
        {"method": "global", "percentile": 90},
        {"method": "knn", "k": 5},
        {"method": "hybrid", "k": 5, "min_sim": 0.6},
    ]

    report = run_diagnostic_sweep(levels, param_grid, n_random=len(ddi_edges))
    print(report.to_string(index=False))
