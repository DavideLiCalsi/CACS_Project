from sympy import symbols, log
from math import ceil, floor

def H_2(x):
    #print("Log of",x)
    return (-x*log(x,2)-(1-x)*log(1-x,2))

def epsilon1(R,alpha,delta_0):
    expr=H_2(delta_0)-R*H_2(alpha/R)-(1-R)*H_2( (delta_0-alpha)/(1-R) )
    return expr

def epsilon2(R,v,alpha,b,delta_0):
    expr=R*H_2(alpha/R)+v*H_2( (delta_0-alpha)/(1-R - (b-1)*v ) )
    return expr

def complexity(R,delta_0,alpha,b,v):
    expr = epsilon1(R,alpha,delta_0) + max( 0.5 * epsilon2(R,v,alpha,b,delta_0)+v, epsilon2(R,v,alpha,b,delta_0) -b*v )
    return expr

def compute_e2(delta_0,alpha,b,v,y):
    expr= ((delta_0-alpha)/(1-R-(b-1)*v))*y
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
    e1_high=d #floor(n* min(R,delta_0) )

    print(e1_high)

    best_complexity=-1
    best_param={
        'y':-1,
        'e1':-1,
        'b':-1
    }

    #else:
     #   raise Error

    for e1 in range(e1_low,e1_high+1):
        for y in range(y_low,y_high+1):

            b_high = ceil( (n-k)/y )
            for b in range(b_low,b_high+1):

                alpha=e1/n

                v=y/n

                if delta_0-alpha<(1-R -(b-1)*v) and alpha < min(delta_0,R) and alpha > max(0,delta_0+R-1):
                    print("Trying y:",y,"e1:",e1,"b:",b, end=" Complexity: ")
                    new_complexity=complexity(R,delta_0,alpha,b,v)
                    print(new_complexity)

                    if best_complexity==-1 or new_complexity < best_complexity:
                        best_complexity=new_complexity
                        best_param['b']=b
                        best_param['e1']=e1
                        best_param['y']=y
    
    print(best_param,best_complexity)
    best_alpha=best_param['e1']/n
    best_v=best_param['y']/n
    print("e2",compute_e2(delta_0,best_alpha,b,best_v,best_param['y']))

n=30
k=n//2
R=k/n

if R == 0.5:
    delta_0 = 0.11002786443835955
else:
    raise Error


search_parameters(n,k,12)