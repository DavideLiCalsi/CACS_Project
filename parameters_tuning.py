"""
python script to search optimal parameters: alpha, v, b
"""
from msilib.schema import Error
from sympy import symbols, log

def H_2(x):
    return (-x*log(x,2)-(1-x)*log(1-x,2))

def complexity(R, delta_0, alpha, b, v):
    expr = (1-R)*(1-H_2((delta_0-alpha)/(1-R))) - max(0.5*R*H_2(alpha/R)-0.5*v, (b-H_2((delta_0-alpha)/(1-R-(b-1)*v)))*v)
    return expr


if __name__ == '__main__':

    #x, y = symbols('x y')

    n=20
    k=10
    R=k/n

    if R == 0.5:
        delta_0 = 0.11002786443835955
    else:
        raise Error

    min_alpha = max(0, delta_0+R-1)
    max_alpha = min(delta_0, R)

    #print(min_alpha, max_alpha)

    best_alpha = 0
    best_v = 0
    best_b = 0
    best_complexity = [999999999]
    best_params = [[-1,-1,-1]]

    alpha = min_alpha + 0.001
    count = 0
    
    while alpha < max_alpha:
        v=0
        b=2
        while b <= (n-k):
            while v <= 1-R:    # 0 <= y <= n-k --> 0 <= v*n <= n-k --> 0 <= v <= 1-R
                if (delta_0-alpha) >= 1-R-(b-1)*v:
                    v += 0.01
                    continue
                try:
                    curr_complexity = complexity(R, delta_0, alpha, b, v)
                except Exception as ex:
                    print(ex)
                    continue
                if (curr_complexity > 0 and curr_complexity <= best_complexity[-1]):
                    best_complexity.append(curr_complexity)
                    best_params.append([alpha,b,v])
                    best_alpha = alpha
                    best_b = b
                    best_v = v
                v += 0.01
            b += 1
        alpha += 0.001

for index, val in enumerate(best_complexity):
    print("%d --- %f" %(index, val))

print("---------------------------------")

for index, val in enumerate(best_params):
    print("%d ---" %index, val)


