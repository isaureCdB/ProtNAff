# Author: Sjoerd de Vries, CNRS. Copyright 2023.

import sys
import numpy as np

bead3_G = [4.00146615, 2.03821757, 0.03968358]
bead3_G[2] = 0  # impose base planarity
bead3_G = np.array(bead3_G)

bead3_U = [1.6724041757167758, 1.810610794424034, 0.00135061477068403]
bead3_U[2] = 0  # impose base planarity
bead3_U = np.array(bead3_U)


def mutate_AtoG(coor):
    """mutates A to G in-place"""
    assert coor.shape == (7, 3), coor.shape
    v1 = coor[4] - coor[3]
    v2 = coor[6] - coor[3]
    v3 = np.cross(v1, v2)
    v2a = np.cross(v3, v1)
    x = v1 / np.linalg.norm(v1)
    y = v2a / np.linalg.norm(v2a)
    z = v3 / np.linalg.norm(v3)
    mat = np.stack((x, y, z))

    bead3 = bead3_G.dot(mat) + coor[3]
    coor[5] = bead3


def mutate_CtoU(coor):
    """mutates C to U in-place"""
    assert coor.shape == (6, 3), coor.shape
    v1 = coor[4] - coor[3]
    v2 = coor[5] - coor[3]
    v3 = np.cross(v1, v2)
    v2a = np.cross(v3, v1)
    x = v1 / np.linalg.norm(v1)
    y = v2a / np.linalg.norm(v2a)
    z = v3 / np.linalg.norm(v3)
    mat = np.stack((x, y, z))

    bead3 = bead3_U.dot(mat) + coor[3]
    coor[5] = bead3


if __name__ == "__main__":

    def err(*args, **kwargs):
        print(*args, **kwargs, file=sys.stderr)
        exit(1)

    import argparse

    parser = argparse.ArgumentParser(
        description="Mutate coarse-grain fragment or fragment library"
    )
    parser.add_argument(
        "input_sequence",
        help="Sequence of the input fragment(s). Must consist of A and/or C",
    )
    parser.add_argument(
        "output_sequence",
        help="Sequence of the output fragment(s). Must consist of A, C, G and/or U",
    )
    parser.add_argument(
        "input_data",
        type=argparse.FileType("rb"),
        help="""Input fragment coordinates in .npy format.
The dtype must be float, and the shape must be (X, 3) or (Y, X, 3).
X is the number of coarse-grain beads in the fragment, which must match with the sequence:
7 beads for each A, 6 beads for each C.""",
    )
    parser.add_argument(
        "output_data",
        type=argparse.FileType("wb"),
        help="""Output fragment coordinates in .npy format.
The shape will be the same as the input data.""",
    )

    args = parser.parse_args()

    sizes = {"A": 7, "C": 6}

    totsize = 0
    for char in args.input_sequence:
        if char not in ("A", "C"):
            err("Input sequence must contain only A and C")
        totsize += sizes[char]

    for char in args.output_sequence:
        if char not in ("A", "C", "G", "U"):
            err("Input sequence must contain only A, C, G, U")

    if len(args.input_sequence) != len(args.output_sequence):
        err("Input and output sequence are not of the same length")

    arr = np.load(args.input_data)
    if arr.ndim not in (2, 3) or arr.shape[-1] != 3:
        err(f"Input data has the wrong shape: {arr.shape}")
    
    if arr.shape[-2] != totsize:
        err(f"Input data has the wrong number of atoms for the sequence: {arr.shape[-2]} instead of {totsize}")
    if arr.ndim == 2:
        arr2 = arr[None, :, :]
    else:
        arr2 = arr
    mutate = arr2.copy()

    offset = 0

    for nuc1, nuc2 in zip(args.input_sequence, args.output_sequence):
        size = sizes[nuc1]
        if nuc1 != nuc2:
            if (nuc1, nuc2) not in (("A", "G"), ("C", "U")):
                err(f"Cannot mutate {nuc1} to {nuc2}")
            for frag in range(len(mutate)):
                to_mutate = mutate[frag, offset : offset + size]
                if (nuc1, nuc2) == ("A", "G"):
                    mutate_AtoG(to_mutate)
                elif (nuc1, nuc2) == ("C", "U"):
                    mutate_CtoU(to_mutate)
        offset += size

    if arr.ndim == 2:
        mutate = mutate[0]

    np.save(args.output_data, mutate)
