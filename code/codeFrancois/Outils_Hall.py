@cached_function
def Delta_Mac(self,g):
        selfH = H(self)
        if self == 0:
            return 0
        return s(sum(cmu*g(Un*B_mu(mu))*H(mu) for mu,cmu in selfH))

def e_mn(m,n=None):
    if n==None:
        return E_mn(m,m)
    if m<0: return 0
    else: return E_mn(m,n)
    
@func_persist
def E_mn(m,n=None):
    if n==None:
        return E_mn(m,m)
    else:
        d=gcd(m,n)
        a=m//d
        b=n//d
        if m<0: 
            return 0
        else: 
            return Theta(e[d],a,b)

@cached_function
def Eval_Qmn_q(m,n):
    return q**((m-1)*(n-1)/2)*Eval1(Q_mn(m,n),q+1/q)

def Formel_P_mn(m,n):
    return Formal_P(Partlist(m,n).conjugate(),decalage(m,n))

@cached_function
def Formule_Qmn_q(m,n):
    return q_int(gcd(m,n))/q_int(m)*s(e[n](X*q_int(m)))

def H_mn(m,n):
    d=gcd(m,n)
    a=m//d
    b=n//d
    return Theta(h_hat(d),a,b)

def M_mu(mu,a=1,b=1):
    d=mu.size()
    return Theta((-1)**(d-mu.length())*m(mu),a,b)

def P_mn(k,n):
    _d=gcd(k,n)
    _a=k//_d
    _b=n//_d
    return Theta(p_hat(_d),_a,_b)

@func_persist
def Q(g,mu,c,d):
    mu=Partition(mu)
    if mu.is_empty(): return g
    else: 
        nu=Partition([mu[k] for k in range(1,mu.length())])
        return Q(Qk(g,mu[0]*c,mu[0]*d),nu,c,d)

@cached_function
def Qk(g,m,n):
    if m == 0: return g*e[1]
    elif n == 0: return g-(1-t)*(1-q)*Delta_Mac(g,e[1])
    else : 
        [u,v]=split(m,n)
        [u1,u2]=u
        [v1,v2]=v
        return s(1/((1-t)*(1-q))*(Qk(Qk(g,u1,u2),v1,v2)-Qk(Qk(g,v1,v2),u1,u2)))
    
def Q_mn(m,n):
    d=gcd(m,n)
    a=m//d
    b=n//d
    return Theta(add((-1/(q*t))**iota(mu)*s(mu) for mu in hooks(d)),a,b)

def split(m,n): 
    if m==1: return [vector([m,n-1]),vector([0,1])]
    elif n==1: return [vector([1,0]),vector([m-1,n])]
    elif gcd(m,n)<>1: 
        d=gcd(m,n) 
        v=split(m/d,n/d)[0]
        u=vector([m,n])-v
        return [v,u]
    L=[[i-floor(i*m/n)*n/m,floor(i*m/n),i] for i in range(1,n)]
    s=min([L[i][0]] for i in range(n-1))[0]
    v=vector([[t[1],t[2]] for t in L if t[0]==s][0])
    u=vector([m,n])-v
    d=Matrix([u,v]).det()
    if d<0: return [v,u]
    else: return [u,v]
    
@cached_function    
def Theta(g,c=1,d=1):
    n=g.degree()
    return Formal_F(sum(g.scalar(q_mu(mu))*Q(1,mu,c,d) for mu in Par(n)))
