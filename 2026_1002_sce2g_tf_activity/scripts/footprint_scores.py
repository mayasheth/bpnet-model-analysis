"""Site footprint scores for 2.3 / 2.4: TOBIAS's own (`K562_Xu_score` in the BINDetect overview) or
PRINT's at the same sites (3.1 output), so the two methods are tested on identical sites."""
import glob
import os

import pandas as pd

PRINT_DIR = os.path.join(os.environ["SCRATCH"], "scprinter_k562_xu")
SITE = ["TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_strand"]
_print = None


def _load_print():
    global _print
    if _print is None:
        files = sorted(glob.glob(f"{PRINT_DIR}/print_tfbs_sites.shard*.tsv.gz"))
        names = [os.path.basename(f) for f in files]
        expected = ["print_tfbs_sites.shard1of4.tsv.gz"] + [
            f"print_tfbs_sites.shard{s}of4_sub{j}of4.tsv.gz" for s in (0, 2, 3) for j in range(4)]
        missing = sorted(set(expected) - set(names))
        assert not missing, f"PRINT shards missing: {missing}"
        P = pd.concat([pd.read_csv(f"{PRINT_DIR}/{f}", sep="\t", usecols=SITE + ["motif", "print_score"])
                       for f in expected if f in names], ignore_index=True)
        assert not P.duplicated(SITE + ["motif"]).any(), "duplicate PRINT sites across shards"
        _print = P
    return _print


def site_scores(ov, source, motif):
    """Return the BINDetect overview `ov` with `K562_Xu_score` set from `source` ('tobias' | 'print')."""
    if source == "tobias":
        return ov
    assert source == "print", source
    P = _load_print()
    P = P[P["motif"] == motif].drop(columns="motif")
    out = ov.merge(P, on=SITE, how="left")
    assert len(out) == len(ov)
    frac = out["print_score"].notna().mean()
    assert frac > 0.999, f"{motif}: only {frac:.4f} of sites have a PRINT score"
    out = out[out["print_score"].notna()].copy()
    out["K562_Xu_score"] = out.pop("print_score")
    return out
