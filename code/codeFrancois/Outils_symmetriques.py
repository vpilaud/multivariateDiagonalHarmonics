def addcol(k,f):
    if f==0:
        return 0
    else:
        return add(c*s(plus(mu,conjugate([k]))) for mu,c in f)

def addmu(nu,f):
    if f==0: 
        return 0
    else: 
        return add(c*s(plus(mu,nu)) for mu,c in f)
    
def add_row(k,f):
    if f==0:
        return 0
    else:
        return add(c*s(mu+[k]) for mu,c in f)

def add_rows(rho,f):
    return add(c*s(sorted(list(mu)+list(rho),reverse=true)) for mu,c in f)

def alt(f):
    if f==0: return 0
    else: 
        n=f.degree()
        return f.scalar(e[n])

def dim(f): 
    if f==0: return 0
    else:
        n=f.degree()
        return f.scalar(h[1]**n)

@func_persist
def Delta(g,f):
    if f==0:
        return 0
    else:
        n=f.degree()
        return Formal_F(s(add(c*g(B_mu(mu)*Un)*H(mu) for mu,c in H(f))))

def Deltaq(g,f):
    n=f.degree()
    return s(add(c.substitute({t:1/q})*g(B_muq(mu)*Un)*Hq(mu) for mu,c in H(f)))

@cached_function
def delta_ek_emn(k,m,n):
    return Delta(e[k],Eval1(E_mn(m,n),q+t))

@cached_function
def Delta_em_pn_q(m,n):
    return Deltaq(e[m],p[n])

@cached_function
def delta_f_emn(f,m,n):
    return Delta(f,Eval1(E_mn(m,n),q+t))

@func_persist
def DeltaPrime(g,f):
    if f==0:
        return 0
    else:
        n=f.degree()
        return Formal_F(s(add(c*g((B_mu(mu)-1)*Un)*H(mu) for mu,c in H(f))))

@cached_function
def deltaprime_ek_emn(k,m,n):
    return DeltaPrime(e[k],Eval1(E_mn(m,n),q+t))

def Degree_bound(k,F):
    if F==0: return 0
    else: return add(c*s(mu) for mu,c in F if mu.size()<k)
    
def D_k(k,f):
    if k<0:
        return (D_k(0,f.skew_by(p[-k]))-D_k(0,f).skew_by(p[-k]))
    elif k>0:
        return (D_k(0,f*p[k])-D_k(0,f)*p[k])/(1-q^k)/(1-t^k)
    else:
        return f-M*Delta_Mac(f,e[1])
    
def e_Prob(n):
    return 1/mul(add(Kostka(nu,mu) for nu in Partitions(n)) for mu in Partitions(n))

@cached_function
def e_h(k,l):
    return s(e([k])(h[l])).restrict_partition_lengths(2,exact=false) 

def Ell(F):
    if F==0: return 0
    else: return max([ell(Partition(mu)) for mu,c in F])
    
@cached_function
def e_sum_h(j,n):
    return add(s(mul(e_h(c[k],k+1) for k in range(n))).restrict_partition_lengths(2,exact=false) 
               for c in IntegerVectors(j, min_part=0,length=n))

@cached_function
def E_sum_h(j,n):
    return s(colone(j))(add(s[i] for i in range(1,n+1))).restrict_partition_lengths(j,exact=false) 

def extend_first_row(k,f):
    if f==Un: return s[1]
    else: return add(c*s([mu[0]+k]+[mu[j] for j in range(1,mu.length())]) for mu,c in f)

def Formal_coeff(G,mu):
    if G.scalar(s(mu))==0: return 0
    else: return  s(Sym(S.from_polynomial(G.scalar(s(mu)).numerator())))

def Formal_q(pol):
    return add(pol.coefficient({q:k})*s[k] for k in range(pol.degree()+1))

def Foulk(mu,nu):
    return s(mu)(s(nu))-s(nu)(s((mu)))

def Hankel(n,k):
    return matrix([[h[i+j+k] for i in range(n)] for j in range(n)])

def Hankele(n,k):
    return matrix([[e[i+j+k] for i in range(n)] for j in range(n)])

@cached_function
def hh(n):
    return add(s[j](add(h[k] for k in range(1,(n+1)-j+1))) 
               for j in range(1,n+1)).restrict_degree(n,exact=false)        

def h_hat(d):
    return (-1/(q*t))**(d-1)*h[d]

def Hq(mu):
    return s(H(mu)).map_coefficients(lambda c:c.substitute({t:1/q}))

@func_persist
def H_delta(n): 
    mu=Partition([n-k for k in range(n)])
    return Formal_F(s(H(mu)))

def Id(F):
    return F

@cached_function
def Inv_hh(n):
    return Un+add((-1)**j*(hh(n+1-j))**j 
               for j in range(1,n+1)).restrict_degree(n,exact=false)

def Inv_Delta(g,F):
    n=F.degree()
    return s(add((1/(g(B_mu(mu)*Un).scalar(Un))*c*H(mu)
               for mu,c in H(F))))

def isotypique(F,mu,k=None): 
    mu=Partition(mu)
    if k==None: return isotypique(F,mu,2)
    else: return s(Sym(S.from_polynomial(numerator(F[mu])))).restrict_partition_lengths(k,exact=false)

def longueur(f):
    if f==0: return 0
    return max({mu.__len__() for mu,c in s(f)})

@func_persist
def Nabla(g,k=None): 
    if k is None: k=1
    return Formal_F(s(g.nabla(power=k)))

def newFoulkes(n):
    return add(Foulk(hook(n,i),hook(2,i%2)) for i in range(n))

def Omega(k,z=None):
    if z==None: return Omega(k,1)
    else: return add(s[i]*z**i for i in range(1,k+1))

def Omega_even(k,z=None):
    if z==None: return Omega_even(k,1)
    else: return add(s[2*i]*z**i for i in range(1,k+1))

def Omega_inv(k,z=None):
    if z==None: return Omega_inv(k,1)
    else: return add(e[i]*(-z)**i for i in range(k+1))

def Omega_odd(k,z=None):
    if z==None: return Omega_odd(k,1)
    else: return add(s[2*i+1]*z**i for i in range(k))
    
def pi_n(n):
    return add((-1/(q*t))**(n-k-1)*s(hook(n,k)) for k in range(n))

def pol_to_s(pol):
    return add(c*s[d.degree()] for c,d in pol)

def p_hat(d):
    return (-1)**(d-1)*p[d]

@cached_function
def p_k_h(k,n):
    return p[k](add(p(s[j]) for j in range((n+1)//k+1))).restrict_degree(n,exact=false)

def p_nu_h(nu,n):
    return mul(p_k_h(k,n) for k in nu).restrict_degree(n,exact=false)

@cached_function
def q_eval(f,n):
    return f(q_int(n)*Un).scalar(Un)

def q_mu(mu): 
    mu=Partition(mu)
    return (((q*t-1)/(q*t))**(mu.length()))*f(mu)(X*q*t/(q*t-1))

def restrict_to_hooks(F):
    if F==Un or F==0: 
        return F
    else: 
        return add(c*s(mu) for mu,c in F 
               if mu==[] or mu.length()==1 or mu[1]==1)

def restreint(f,k):
    if f==0: return 0*Un
    else: return f.restrict_partition_lengths(k)
    
def restreint_plus_trois(F):
    return add(c*s(mu) for mu,c in F if mu.length()>2)

def sans2(f): return f-f.restrict_partition_lengths(2,exact=false)

def Scalar(f,g):
    if f==0: return 0
    else: return f.scalar(g)

def Schur_hat(mu):
    return (-1/(q*t))**iota(mu)*s(mu)

def SchurProb(n):
    return 1/mul(add(Kostka(mu,nu) for nu in Partitions(n)) for mu in Partitions(n))

def Skew(f,g):
    if f==0: return 0
    else: return f.skew_by(g)
    
def Skew_by_F(f,g):
    if f==q*Un or f==0:
        return f
    else:
        return f.skew_by(g)

def Strip(F):
    if F==0 or F==1:
        return F
    else:
        return  add(c*t**(strip(mu).size())*q**(mu.size()) for mu,c in F)
    
def Stripx(F):
    if F==0 or F==1:
        return F
    else:
        return  add(c*t**(strip(mu).size())*x**(mu.size()) for mu,c in F)

def s_hook(n,k):
    if k>n-1 or k<0:
        return 0
    else:
        return s(hook(n,k))

def s_mu_h(mu):
    mu=Partition(mu)
    n=mu.size()
    return s(add(c*p_nu_h(nu,n) for nu,c in p(s(mu))).restrict_degree(n,exact=false))

def to_eigenH(F):
    G=H(eval_tensor1(F,q+t))
    n=G.degree()
    return add((G.coefficient(mu))/(H(e[n]).coefficient(mu))*H(mu) for mu in Partitions(n))

def ToPol(F):
    return  add(c*t**(strip(mu).size())*q**(mu.size()) for mu,c in F)

def UnStrip(F):
    return  add(add(c*s([qd.degree()-mu.size()]+list(mu)) for c,qd in cq.numerator()) for mu,cq in F)

def updown(f,a,b):
    return Skew(addmu(colone(a),f),e[b])

def ValPropreH(F):
    n=F.degree()
    G=H(F)
    E=H(e[n])
    return add(G.coefficient(mu)/E.coefficient(mu)*H(mu) for mu in Partitions(n))