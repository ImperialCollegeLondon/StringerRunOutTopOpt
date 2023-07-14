import math
import numpy as np


def log(msg, level=0, title="", padBefore=0, padAfter=0):
    levelStr = "--"*(level)
    if levelStr:
        levelStr += " "
    padBeforeStr = "\n"*padBefore
    padAfterStr = "\n"*padAfter
    if title != "":
        stringVal = "{0}({1}){2}{3}{4}".format(
            padBeforeStr, title, levelStr, msg, padAfterStr)
    else:
        stringVal = "{0}{1}{2}{3}".format(
            padBeforeStr, levelStr, msg, padAfterStr)
    print(stringVal)


def log_list(name, list):
    levelStr = " "*6
    print("")
    print(levelStr + ">> " + name + ':')
    for val in list:
        print(levelStr + val)
    print("")


def str_from_list(alist, index=False, ids=[]):
    out = []
    if type(alist[0]) is list:
        if (not ids):
            ids = list(range(1, len(alist)+1))
        for isl, sublist in enumerate(alist):
            if len(sublist) > 16:
                lines = math.ceil(len(sublist) / 16)
                current = 0
                l = []
                for iln in range(lines):
                    if iln < lines - 1:
                        l += [sublist[current:(current+16)]]
                        current += 16
                    else:
                        l += [sublist[current:]]
                for il, li in enumerate(l):
                    if il == 0:
                        out += [str(ids[isl]) + ', ' + str(li)[1:-1] +
                                ',\n'] if index else [str(li)[1:-1] + ',\n']
                    else:
                        out += ['    ' + str(li)[1:-1] + '\n']
            else:
                out += [str(ids[isl]) + ', ' + str(sublist)[1:-1] +
                        '\n'] if index else [str(sublist)[1:-1] + '\n']
    elif len(alist) > 16:
        l = []
        lines = math.ceil(len(alist) / 16)
        current = 0
        for iln in range(lines):
            if iln < lines - 1:
                l += [alist[current:(current+16)]]
                current += 16
            else:
                l += [alist[current:]]
        for il, li in enumerate(l):
            if il == 0:
                out += [str(li)[1:-1] + ',\n']
            else:
                out += ['    ' + str(li)[1:-1] + '\n']
    else:
        if (not ids):
            ids = [1]
        out += [str(ids[0]) + ', ' + str(alist)[1:-1] +
                '\n'] if index else [str(alist)[1:-1] + '\n']

    return out


def list_from_str(astr, atype, index=False):
    out = []
    if type(astr) is list:
        alist = astr
    else:
        alist = astr.split('\n')
    if len(alist) > 1:
        for entry in alist:
            out += [[atype(x.strip()) for x in entry.strip().split(",") if x != ''][1:]
                    ] if index else [atype(x.strip()) for x in entry.strip().split(",") if x != '']
    else:
        out += [[atype(x.strip()) for x in astr[0].strip().split(",") if x != ''][1:]
                ] if index else [atype(x.strip()) for x in astr[0].strip().split(",") if x != '']
    return out


def dotp(lista, listb):
    prod = [x * y for (x, y) in zip(lista, listb)]
    return sum(prod)


def cross(lista, listb):
    retlist = [
        (lista[1]*listb[2]-lista[2]*listb[1]),
        -(lista[0]*listb[2]-lista[2]*listb[0]),
        (lista[0]*listb[1]-lista[1]*listb[0])
    ]
    return [x if abs(x) > 1e-9 else 0 for x in retlist]


def matvecmul(mat, vec):
    return [dotp(mat[0], vec), dotp(mat[1], vec), dotp(mat[2], vec)]


def norm(list):
    # print(list)
    # print(math.sqrt(sum([x*x for x in list])))
    return math.sqrt(sum([x*x for x in list]))

def sumv(lista, listb):
    return [lista[x] + listb[x] for x in range(len(lista))]

def diffv(lista, listb):
    return [lista[x] - listb[x] for x in range(len(lista))]

def distv(lista, listb):
    return norm(diffv(lista, listb))

def scalev(lista, scalar):
    return [scalar*x for x in lista]

def normalise(list):
    return [x/norm(list) for x in list]


def transpose(list):
    return np.array(list).T.tolist()
