from sympy import *
from math import ceil, floor, comb
from scipy.special import binom
import sys

def H_2(x):
    return (0-x*log(x,2)-(1-x)*log(1-x,2))

def epsilon1(R,alpha,delta_0):
    expr = H_2(delta_0)-R*H_2(alpha/R)-(1-R)*H_2( (delta_0-alpha)/(1-R) )
    return expr

def epsilon2(R,v,alpha,b,delta_0):
    expr = R*H_2(alpha/R)+v*H_2( (delta_0-alpha)/(1-R - (b-1)*v ) )
    return expr

def complexity(R,delta_0,alpha,b,v):
    expr = epsilon1(R,alpha,delta_0) + max( 0.5 * epsilon2(R,v,alpha,b,delta_0)+v, epsilon2(R,v,alpha,b,delta_0) -b*v )
    return expr

def compute_e2(delta_0,alpha,b,v,y):
    expr = ((delta_0-alpha)/(1-R-(b-1)*v))*y
    return expr

def gilbertVashamovDistance(n, k, q):

    y = pow(q, n-k);
    distance = 1;
    partial_sum = 0;

    while(partial_sum <= y):
        distance += 1
        partial_sum += comb(n,distance-1)*pow(q-1,distance-1);
    
    return distance


def computeIterations(n,k,e,GV_dist):
    expr= n*log(n,2)*( binom(n,GV_dist)/( binom(k,e)*binom(n-k,GV_dist-e) ) )
    return expr


#def complexity(R, delta_0, alpha, b, v):
#     if alpha > delta_0:
#        print("Error")
#    expr = (1-R)*(1-H_2((delta_0-alpha)/(1-R))) - max(0.5*R*H_2(alpha/R)-0.5*v, (b-H_2((delta_0-alpha)/(1-R-(b-1)*v)))*v)
#    return expr

def search_parameters(n,k,d):

    R=k/n
    if R == 0.5:
        delta_0 = 0.11002786443835955

    y_low=1
    y_high=n-k
    b_low=1
    b_high=None # it depends on y
    e1_low=1
    e1_high= floor(n* min(R,delta_0) )

    print(e1_high)

    best_complexity=-1
    best_param={
        'y':-1,
        'e1':-1,
        'b':-1
    }

    for e1 in range(e1_low,e1_high+1):
        for y in range(y_low,y_high+1):

            b_high = ceil( (n-k)/y )
            for b in range(b_low,b_high+1):

                alpha=e1/n

                v=y/n

                if delta_0-alpha<(1-R -(b-1)*v) and alpha < min(delta_0,R) and alpha > max(0,delta_0+R-1):
                    print("Trying y:",y,"e1:",e1,"b:",b, end=" Complexity: ")
                    new_complexity=complexity(R,delta_0,alpha,b,v)
                    print(N(new_complexity), end=" Iterations: ")
                    print( ceil( N(computeIterations(n,k,e1,d))) )

                    if best_complexity==-1 or new_complexity < best_complexity:
                        best_complexity=N(new_complexity)
                        best_param['b']=b
                        best_param['e1']=e1
                        best_param['y']=y
    
    print("FOUND PARAMETERS")
    print("\ty:",best_param['y'])
    print("\te1:",best_param['e1'])
    print("\tb:",best_param['b'])
    best_alpha=best_param['e1']/n
    best_v=best_param['y']/n
    print("\te2:",compute_e2(delta_0,best_alpha,b,best_v,best_param['y']))
    print("Expected complexity:", N(complexity(0.5,delta_0,best_alpha,best_param['b'],best_v)))
    print("Expected number of iterations:", ceil( N(computeIterations(n,k,best_param['e1'],d))) )



if __name__ == "__main__":

    USAGE = "usage: python find_param <n>"
    EXAMPLE =  "example: python find_param 30"

    args = sys.argv[1:]

    if len(args) < 1:
        print(USAGE)
        print(EXAMPLE)
        sys.exit(-1)

    n = int(args[0])
    k=n//2
    R=k/n
    print("CODE RATE",R)
    if R == 0.5:
        delta_0 =n* 0.11002786443835955
    else:
        raise Exception("Assumption: n == 2*k")

    search_parameters(n, k, gilbertVashamovDistance(n, k, 2) )