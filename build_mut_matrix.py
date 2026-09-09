"""Build the high-dimensional mutation matrix from parse_mutations.py output.
Isolates x specific amino-acid alleles (binary, alleles seen in >=3 isolates).
Outputs cryptic/mut_matrix.npz + mut_rows.txt + mut_cols.txt."""
import collections, numpy as np, scipy.sparse as sp

rows, sets = [], []
for line in open("cryptic/mut_tokens.jsonl"):
    parts = line.rstrip("\n").split("\t")
    rows.append(parts[0])
    sets.append(set(t for t in (parts[1].split(",") if len(parts) > 1 and parts[1] else []) if t))

cnt = collections.Counter(t for s in sets for t in s)
vocab = sorted([t for t, c in cnt.items() if c >= 3])          # drop singletons/doubletons
vidx = {t: i for i, t in enumerate(vocab)}

ri, ci = [], []
for i, s in enumerate(sets):
    for t in s:
        j = vidx.get(t)
        if j is not None:
            ri.append(i); ci.append(j)
M = sp.csr_matrix((np.ones(len(ri), np.float32), (ri, ci)),
                  shape=(len(rows), len(vocab)), dtype=np.float32)

sp.save_npz("cryptic/mut_matrix.npz", M)
open("cryptic/mut_rows.txt", "w").write("\n".join(rows))
open("cryptic/mut_cols.txt", "w").write("\n".join(vocab))
print("matrix:", M.shape, "| nnz:", M.nnz,
      "| density %.4f" % (M.nnz / (M.shape[0] * M.shape[1])))
