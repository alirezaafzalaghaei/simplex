import numpy as np
from tabulate import tabulate
from sympy import Matrix
import sympy as sp

def solve(s, method='bigm'):
    A, b, c, B, two_phase = preprocess(s)

    c1 = Matrix([1 if i == None else 0 for i in c])
    c2 = Matrix([i for i in c if i != None])
    if two_phase and method=='twophase':
        print("Phase 1:")
        z, A, b, c, x, B, BT = tableau_simplex(c1, A, b, B)
        A = A[:,:c2.shape[0]]
        print("Phase 2:")
        res = tableau_simplex(c2, A, b, B)    
    elif two_phase and method=='bigm':
        M = sp.Symbol('M', positive=True, zero=False)
        c = Matrix([i if i != None else M for i in c])
        res = tableau_simplex(c, A, b, B)        
    else:
        res = tableau_simplex(c2, A, b, B)
    return res

def print_tableau(T, B):

    m, n = T.shape
    R = np.empty((m + 1, n + 1), dtype=object)
    R[1:, 1:] = [list(map(str, T.row(i))) for i in range(m)]
    R[0, 0] = ''
    R[0, 1] = R[1, 0] = 'z'
    R[0, 2:n] = list(map(lambda i: 'x_%d' % i, range(1, n - 1)))
    R[0, -1] = 'RHS'

    R[2:, 0] = list(map(lambda i: 'x_%d' % i, B))

    return tabulate(R,
                    tablefmt='fancy_grid',
                    headers="firstrow",
                    numalign='center',
                    stralign='center')


def nCm(n, m):
    import math
    return math.factorial(n) // math.factorial(m) // math.factorial(n - m)

def showLP(A, b, c):
    m, n = A.shape
    row = ''
    row += 'min ' + ' + '.join(
        map(lambda i: '%.2g x_%d' % (c[i - 1], i), range(1, n + 1))) + '\n'
    A = [list(map(str, A.row(i))) for i in range(m)]
    for r, _ in zip(A, b):
        row += '    ' + ' + '.join(
            map(lambda i: '%s x_%d' %
                (r[i - 1], i), range(1, n + 1))) + ' <= ' + str(_) + '\n'

    row = row.replace('+ -', '-')
    row += ''
    return row


def preprocess(s):
    s = list(filter(None, s.strip().split('\n')))
    ch = s.pop(0).split(' ')
    c = Matrix(ch[1:])
    
    if ch[0] == 'max':
        c = -c
    A = []
    b = []
    typ = []
    for r in s:
        t = r.split(' ')
        A.append(t[:-2])
        b.append(t[-1])

        ty = t[-2].strip()

        if ty == '<=':
            typ.append(0)
        elif ty == '>=':
            typ.append(1)
        else:
            raise ValueError(t[-2])

    A = Matrix(A)
    b = Matrix(b)

    m, n = A.shape
    Ahat = Matrix.hstack(A, sp.eye(m))
    c = Matrix.vstack(c, sp.zeros(m, 1))
    eye = sp.eye(m)
    two_phase = False
    B = []
    for i in range(m):
        if typ[i] == 1:
            Ahat[i, n + i] = -1
            #             b[i] = -b[i] # if B-1b <0?
            Ahat = Ahat.col_insert(Ahat.shape[1], eye.col(i))
            c = c.row_insert(c.shape[0], Matrix([None]))
            two_phase = True
            B.append(Ahat.shape[1])
        else:
            B.append(n + i + 1)

    if not two_phase:
        B = (n + np.arange(m, dtype=int) + 1).tolist()
    return Ahat, b, c, B, two_phase

def argmax(z):       
    try:        
        tmp = np.array(z).flatten()        
        i = tmp.argmax()
        return i, tmp[i]
    except:
        from collections import defaultdict
        
        M = sp.Symbol('M', positive=True, zero=False)
        
        d1 = defaultdict(int)
        d2 = {}
        for i in range(2,8):
            pp = z.subs({M:10**i})
            
            j,zj = argmax(pp)
            d1[j] += 1
            d2[j] = zj
        mx = max(d1.items(),key=lambda x:x[1])[0]
        return mx, d2[mx]
        
def tableau_simplex(c, A, b, Basis):
    BasisTrace = []
    m, n = A.shape
    B = Basis

    max_iterations = nCm(n, m)
    print('max iterations:', max_iterations)

    T = sp.zeros(m + 1, n + 2)

    T[0, 0] = 1
    T[0, 1:-1] = -c.transpose()
    T[1:, 1:-1] = A
    T[1:, -1] = b

    print('%dth iteration:' % 1)
    print(print_tableau(T, B))
    BasisTrace.append(set(B))

    for base in B:
        if T[0, base] != 0:
            print('Making z_%d zero.' % (base))
            Ma = T[1:, base]            
            for i in range(Ma.shape[0]):
                if Ma[i] == 1:
                    break
            else:
                print('There is no 1 in this column.')
                break
            T[0, :] += -T[1 + i, :]*T[0, base]
            print(print_tableau(T, B))
    print('Beginning Problem solving:')
    for iteration in range(2, max_iterations + 1):
        print('%dth iteration:' % iteration)
        zj_cj = T[0, 1:-1]
        index_max,max_val = argmax(zj_cj)
        
        if max_val <= 0:
            print('Negative Zj-Cj. Optimal Solution found. Z=%.2g' % T[0, -1])
            break

        entering_index = index_max

        i = 0

        min_ratio = np.array([
            T[1 + i, -1] /
            T[1 + i, index_max + 1] if T[1 + i, index_max + 1] > 0 else np.inf
            for i in range(m)
        ])

        if np.all(np.isinf(min_ratio.astype(float))):
            print('Unbounded Problem.')
            break
        exiting_index = min_ratio.argmin()

        print('Entering variable: x_%d' % (entering_index + 1))

        print('Exiting  variable: x_%d' % (B[exiting_index]))

        print('Pivot: %s' % T[exiting_index + 1, entering_index + 1])

        B.remove(B[exiting_index])
        B.insert(exiting_index, entering_index + 1)
        if set(B) not in BasisTrace:
            BasisTrace.append(set(B))
        else:
            BasisTrace.append(set(B))
            print('A Basis Cycle has been found! Trace:')
            print(*BasisTrace, sep='\n')
            break
        for i in range(m + 1):  # handle first row
            if i == exiting_index + 1:
                continue
            T[i, :] += -T[i, entering_index +
                          1] / T[exiting_index + 1,
                                 entering_index + 1] * T[exiting_index + 1, :]
        T[exiting_index + 1, :] /= T[exiting_index + 1, entering_index + 1]

        
        print(print_tableau(T, B))

    c = T[0, 1:-1]
    z = T[0, -1]
    A = T[1:, 1:-1]
    b = T[1:, -1]
    x = [b[B.index(i)] if i in B else 0 for i in range(1, 1 + n)]
    return z, A, b, c, x, B, BasisTrace



