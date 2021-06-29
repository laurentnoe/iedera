#!/usr/bin/env python3

import math

flatten = lambda l: [item for sublist in l for item in sublist]

def tamura1992substitutionMatrix(gcPercent, kappa, pam):
    s = 0.005 * gcPercent  # prob(c) = prob(g)
    w = 0.5 - s            # prob(a) = prob(t)
    t = pam * 0.01 / (4 * w * s * kappa + 0.5)
    i = math.exp(-t)
    j = 2 * math.exp(-(kappa + 1) * t / 2)
    wwr = 1 + i + j*s/w
    ssr = 1 + i + j*w/s
    tir = 1 + i - j
    tvr = 1 - i

    wwp = w * w * wwr
    ssp = s * s * ssr
    percentIdentity = 200 * (wwp + ssp)
    #print("percent identity: " + str(percentIdentity))

    #                       A         C         G         T
    probFromRowToColumn = [[w * wwr,  s * tvr,  s * tir,  w * tvr],  # A
                           [w * tvr,  s * ssr,  s * tvr,  w * tir],  # C
                           [w * tir,  s * tvr,  s * ssr,  w * tvr],  # G
                           [w * tvr,  s * tir,  s * tvr,  w * wwr]]  # T
    return probFromRowToColumn




def acgt_cross_tamuratamura1992substitutionMatrix(gcPercent, kappa, pam):
    s = 0.005 * gcPercent  # prob(c) = prob(g)
    w = 0.5 - s            # prob(a) = prob(t)
    original = [w, s, s, w]
    mutated = tamura1992substitutionMatrix(gcPercent, kappa, pam)    
    result = []
    for i,o in enumerate(original):
        result.append([])
        for m in mutated[i]:
            result[i].append(o * m)
    return result




import argparse

if __name__ == "__main__":
    usage = "%(prog)s [options] length"
    descr = "provides the TAM92 associated probabilities for the \"./iedera -iupac \'-f\'\" command-line"
    ap = argparse.ArgumentParser(usage=usage, description=descr,
                                 formatter_class=
                                 argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("-p", "--pam", type=float, default=30.0,
                    help="PAM distance between related sequences")
    ap.add_argument("-k", "--kappa", type=float, default=1.0, metavar="K",
                    help="transition/transversion rate ratio")
    ap.add_argument("--gc", type=float, default=50.0, metavar="PERCENT",
                    help="percent G+C")
    args = ap.parse_args()


    t = acgt_cross_tamuratamura1992substitutionMatrix(args.gc, args.kappa, args.pam)
    t = flatten(t)
    print(','.join('{}'.format(e) for e in t))




