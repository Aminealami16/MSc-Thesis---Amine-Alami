import numpy as np

def bisection(f, a, b, tol=1e-6, max_iter=100):
    fa, fb = f(a), f(b)
    if fa * fb > 0:
        return None  # no sign change, no guaranteed root
    for _ in range(max_iter):
        c = 0.5*(a+b)
        fc = f(c)
        if abs(fc) < tol or (b - a)/2 < tol:
            return c
        if fa * fc < 0:
            b, fb = c, fc
        else:
            a, fa = c, fc
    return 0.5*(a+b)

def find_all_roots(f, a, b, step=0.1, tol=1e-6):
    xs = np.arange(a, b+step, step)
    roots = []

    for i in range(len(xs)-1):
        x1, x2 = xs[i], xs[i+1]
        y1, y2 = f(x1), f(x2)

        # detect sign change or hitting exactly zero
        if y1 == 0:
            roots.append(x1)
        elif y1 * y2 < 0:  # sign change
            root = bisection(f, x1, x2, tol=tol)
            if root is not None:
                # avoid duplicates due to resolution
                if not roots or abs(root - roots[-1]) > tol:
                    roots.append(root)


                    
    
    return roots

