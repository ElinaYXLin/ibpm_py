"""Parse UIUC LSAT-format .DRG (alpha/Cl/Cd) and .LFT (alpha/Cl/Cm) files
into a list of {Re, alphas, cls, cds_or_cms} blocks.
"""
import re


def parse_blocks(path, kind):
    """kind: 'drg' (alpha, Cl, Cd[, spanwise Cds]) or 'lft' (alpha, Cl, Cm)."""
    text = open(path).read()
    lines = text.splitlines()
    blocks = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("Average Reynolds"):
            re_val = float(lines[i + 1].strip())
            # find "Number of angles of attack:" then count, then header, then data
            j = i + 2
            while not lines[j].strip().startswith("Number of angles"):
                j += 1
            n = int(lines[j + 1].strip())
            j += 2
            while "alpha" not in lines[j]:
                j += 1
            j += 1
            alphas, cls_, seconds = [], [], []
            for k in range(n):
                parts = lines[j + k].split()
                alphas.append(float(parts[0]))
                cls_.append(float(parts[1]))
                seconds.append(float(parts[2]))
            blocks.append(dict(Re=re_val, alpha=alphas, Cl=cls_,
                                **({"Cd": seconds} if kind == "drg" else {"Cm": seconds})))
            i = j + n
        else:
            i += 1
    return blocks


def nearest_block(blocks, target_re):
    return min(blocks, key=lambda b: abs(b["Re"] - target_re))


if __name__ == "__main__":
    import sys
    b = parse_blocks(sys.argv[1], sys.argv[2])
    for blk in b:
        print(blk["Re"], len(blk["alpha"]))
