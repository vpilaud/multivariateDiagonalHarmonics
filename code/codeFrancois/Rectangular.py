from sage.combinat.tamari_lattices import GeneralizedTamariLattice

def Aire(mu,m,n):
    mu=Partition(mu)
    return staircase(m,n).size()-mu.size()

def B_mu(mu):
    mu=Partition(mu)
    return sum((t**i*q**j) for i,j in mu.cells())

def compact(mu):
    if mu==Partition([]): return 0
    else: return (''.join(mystr(i) for i in mu))


Partition._latex_= compact
Composition._latex_= compact


def a(mu,m,n=None):
    if n==None:
        return a(mu,m,m)
    else:
        nu=list(mu)+[0 for i in range(n-mu.length())]
        dmn=staircase(m,n)
        rho=list(dmn)+[0 for i in range(n-dmn.length())]
        return [rho[i]-nu[i] for i in range(n)]

def area(mu,m,n):
    return staircase(m,n).size()-mu.size()

@cached_function
def bigraded(m,n):
    return s(add((-1)**(n-ell(mu))*(m*n+1)**(ell(mu)-1)*p(mu)/zee(mu) for mu in Partitions(n)))

def binqn(n):
    return add(stirling_number2(n,k)*factorial(k)*q**k for k in range(n+1))

def Bizley(m,n):
    d=gcd(m,n)
    a=m/d
    b=n/d
    return s(add(1/zee(mu)*mul(1/a*e[k*b](k*a*X) for k in mu) for mu in Partitions(d))) 

def bizley(m,n):
    d=gcd(m,n)
    a=m/d
    b=n/d
    return add(1/zee(mu)*mul(1/(a+b)*binomial(k*(a+b),k*a) for k in mu) for mu in Partitions(d))

def bizley_e(m,n):
    d=gcd(m,n)
    a=m/d
    b=n/d
    return add(1/zee(mu)*mul(1/a*e[k*b](k*a*X) for k in mu) for mu in Partitions(d))

def bizley_h(m,n):
    d=gcd(m,n)
    a=m/d
    b=n/d
    return (-1)**(n-1)*add((-1)**(d-mu.length())/zee(mu)*mul(1/a*e[k*b](k*a*X) for k in mu) 
                          for mu in Partitions(d))

def bottom(n):
    if n==0: return []
    else: return  bottom(n-1)+[1,0]
    
def Catq(m,n):
    return add(q**area(mu,m,n) for mu in Dyck(m,n))

def catalan(n): return (1/(n+1)*binomial(2*n,n))

def CC(k,j): return SR(k).binomial(j, hold=True)

@cached_function
def chaines(n,k):
    return [[nu for nu in Intervale(mu,[0],n) if longueur_mu(nu,n)==k] for mu in debuts(n,k)]

def chemin_christofell(m,n):
    r=staircase(m,n)
    return 1

def Christofell_h(m,n): 
    _C=[floor(i*(m/n-1/(n+1)**2)) for i in range(n)]
    A=Set(_C)
    return Partition(sorted([_C.count(i) for i in A],reverse=1))

def Christofell_e(_a,_b): 
    _C=[floor(i*(_a/_b)) for i in range(_b)]
    _A=Set(_C)
    return Partition(sorted([_C.count(i) for i in _A],reverse=1))

def Christo(a,b):
    return Partition([a for k in range(b//a)]+[b-a*(b//a)])

def coeff_nu(nu):
    nu=Partition(nu)
    n=nu.size()
    return add((sigma(mu,n)*(e_rho_mu(mu,n).scalar(nu))) for mu in Dyck(n,n))

@cached_function
def composante(mu,n,k):
    return restreint(sigma_prime(mu,n),k).skew_by(e[k])

def debuts(n,k):
    return [mu for mu in Dyck(n) if longueur_mu(mu,n)==k and minimal(mu,n,k)]


@cached_function
def defaut_precede(n):
    return [mu for mu in Dyck(n) 
            if sigma_prime(mu,n).skew_by(e[1])<>add(sigma_prime(nu,n) 
                                                    for nu in precede(mu,n))]

def descentes(mu,n):
    mu=Partition(mu)
    L=list(mu)+[0 for k in range(n-mu.length())]
    return Set([i for i in range(n-1) if L[i]>L[i+1]])

def DimTam(m,n): return (m+1)**n*(m*n+1)**(n-2)

def DimAltTam(m,n): return (m+1)/(n*(m*n+1))*binomial((m+1)**2*n+m,n-1)


def Dyck(m,n=None):
    if n==None: return Dyck(m,m)
    else:
        N=staircase(m,n).size()
        return [mu for k in range(N+1) for mu in Partitions(k, outer=staircase(m,n))]

def Dyck_forme(mu):
    mu=Partition(mu)
    n=mu.size()
    return [g for g in Dyck(n,n) if e_rho_mu(g,n)==e(mu)]


@cached_function
def exces(nu,n):
    return add(sigma_prime(mu,n) for mu in to_ideal(nu,n))-sigma_prime(nu,n)(X+Un)


def e_rho_mu(mu,n):
    mu=Partition(mu)
    return mul(e[k] for k in risers(mu,n))

def E_coeff_mu(mu):
    mu=Partition(mu)
    n=mu.size()
    return add(Kostka(mu.conjugate(),riserpart(nu,n))*sigma_prime(nu,n)
               for nu in Dyck(n))(X-Un)

def fins(n,k):
    return [mu for mu in Dyck(n) if longueur_mu(mu,n)==k and maximal(mu,n,k)]


def F_coeff_mu(rho):
    rho=Partition(rho)
    n=rho.size()
    return add(sigma_prime(mu,n) for mu in Dyck(n) 
               if riserpart(nu,n)== rho)

@cached_function
def Formule_F_mn(n):
    return add(Tenseur(sigma(mu,n),e_rho_mu(mu,n)) 
                                  for mu in Dyck(n,n))

def Formule_prime_F_mn(m,n=None):
    if n==None:
        return Formule_prime_F_mn(m,m)
    elif m==n-1:
        return add(Tenseur(sigma_prime(mu,n-1),e_rho_mu(mu,n)) 
                                  for mu in Dyck(m,n))
    else:
        return add(Tenseur(sigma_prime(mu,m,n),e_rho_mu(mu,n)) 
                                  for mu in Dyck(m,n))

@cached_function
def Haglund(m,n,k):
    return add(Jim(mu,m,n,k)*e_rho_mu(mu,n) for mu in Dyck(m,n))

@cached_function
def Haglund_F(n,k):
    return Formal_F(add(Jim(Partition(mu),n,n,k)*llt for (mu,llt) in ListeLLT[n-1]))

def hook_rel(k,n):
    return add(Kostka(Partition([n-k]+list(colone(k))),nu)*F_coeff_mu(nu)
               for nu in Partitions(n))

def Intervale(mu,nu,n):
    nu=Partition(nu)
    mu=Partition(mu)
    v=to_dyckword(mu,n)
    w=to_dyckword(nu,n)
    return list(reversed([to_partition(mu) 
                          for mu in DyckWord(v).tamari_interval(w).dyck_words()]))

@cached_function
def Jim(mu,m,n,k):
    b=a(mu,m,n)
    return add(mul(q**(b[i]) for i in J) for J in Set(range(n-1)).subsets() 
               if J.issuperset(descentes(mu,n)) and J.cardinality()==k)


@cached_function
def longueur_mu(mu,n):
    return n-risers(mu,n).__len__()

def maximal(mu,n,k):
    return all([longueur_mu(nu,n)>k for nu in succede(mu,n)])


def minimal(mu,n,k):
    return all([longueur_mu(nu,n)<k for nu in precede(mu,n)])


def Park(mu,n):
    return StandardSkewTableaux([plus(mu,colone(n)),mu])

def ParkFrob_mu(mu,n):
    return s(plus(mu,colone(n))).skew_by(s(mu))

def ParkFrob(m,n):
    return add(ParkFrob_mu(mu,n) for mu in Dyck(m,n))

def ParkFrobq(m,n):
    N=staircase(m,n).size()
    return add(q**area(mu,m,n)*ParkFrob_mu(mu,n) for mu in Dyck(m,n))

def Part_un_Fmn(m,n):
    return add(Tenseur(s[Aire(mu,m,n)],mul(e[k] for k in risers(mu,n))) for mu in Dyck(m,n))

def plus(mu,nu):
    a=add(1 for j in mu)
    b=add(1 for j in nu)
    if a<=b: 
        return Partition([mu[j]+nu[j]
                          for j in range(a)]+[nu[j] 
                                              for j in range(a,b)])
    else:
        return Partition([mu[j]+nu[j] 
                          for j in range(b)]+[mu[j]
                                              for j in range(b,a)])
@cached_function
def precede(mu,n):
    rep=[]
    for nu in to_ideal(mu,n):
        if Intervale(nu,mu,n).__len__()==2:
            rep=rep+[nu]
    return rep

def q_cat(n):
    return add(q^k*Partitions(k, outer=[n-1-k for k in range(n-1)]).cardinality() 
               for k in range(binomial(n,2)+1))

def q_contain(mu):
    mu=Partition(mu)
    n=mu.size()
    return add(q^k*Partitions(k, outer=mu).cardinality() for k in range(n+1))

def risers(mu,n):
    nu=list(mu)+[0 for i in range(n-mu.length())]
    return [nu.count(i) for i in Set(nu)]

def riserpart(mu,n):
    return Partition(list(reversed(sorted(risers(mu,n)))))

def rho(a,b):
    a=Integer(a)
    b=Integer(b)
    return add(floor(i*(a/b)) for i in range(b))

def sigma_k(k,mu,m,n=None):
    if n==None:
        return sigma_k(k,mu,m,m)
    else:
        return pol_to_s((q**Aire(mu,m,n)*e[k](Un*add(q**(-a(mu,m,n)[i]) 
                                          for i in range(n-1) 
                                           if a(mu,m,n)[i]>a(mu,m,n)[i+1])).scalar(Un)).numerator())
def sigma(mu,m,n=None):
    if n==None:
        return  sigma(mu,m,m)
    else:
        mu=Partition(mu)
        if mu==Partition([0]): return CalA_mn(m,n)
        elif m==n and mu==colone(n-1) or mu==Partition([n-1]): return sigma([0],m,n-1)
        elif m==n and mu.length()==1:
            k=mu[0]
            return addmu([n-1-k],sigma([n-1],m,n))
        elif m==n and (mu.conjugate()).length()==1: return sigma(mu.conjugate(),m,n)
        elif mu==staircase(m,n): return Un
        else:
            rep=0
            for j in range(n-1):
                k=n-1-j
                temp=sigma_k(k,mu,m,n)-restreint(Skew(rep,e[k]),1)
                rep=rep+addmu(colone(k),temp)
        return rep

def sigma_prime(mu,m,n=None):
    if n==None:
        return sigma_prime(mu,m,m)
    elif m==5 and n==5:
        if mu==Partition([1]): return sigma(mu,n)+s[3,3]+s[5,2]
        elif mu==Partition([1,1]): return sigma(mu,n)+s[4,2]
        elif mu==Partition([2]): return sigma(mu,n)+s[4,2]
        elif mu==Partition([1,1,1]): return sigma(mu,n)+s[3,2]
        elif mu==Partition([3]): return sigma(mu,n)+s[3,2]
        elif mu==Partition([2,1]): return sigma(mu,n)+s[3,2]
        elif mu==Partition([2,2]): return sigma(mu,n)+s[2,2]
        elif mu==Partition([3,1,1]): return sigma(mu,n)+s[2,2]
        else: return sigma(mu,m,n)
    else:
        return sigma(mu,m,n)
    

def staircase(m,n):
    return Partition([(m*(n-k))//n for k in range(1,n+1)])

@cached_function
def succede(mu,n):
    rep=[]
    for nu in Intervale(mu,[0],n):
        if Intervale(mu,nu,n).__len__()==2:
            rep=rep+[nu]
    return rep

def teste_ordre(mu,nu,n):
    try:
        Intervale(mu,nu,n)
    except ValueError:
        return(false)
    else:
        return(true)

def test_sigma(n):
    return Set([sigma(mu,n)==sigma(mu.conjugate(),n) for mu in Dyck(n,n)])

def test_sigma_prime(n):
    for mu in Dyck(n,n):
        if sigma_prime(mu,n)<>sigma_prime(mu.conjugate(),n):
            show((mu,e_rho_mu(mu,n),mu.conjugate(),e_rho_mu(mu.conjugate(),n)))

def to_area(mu,n):
    v=[0 for i in range(n-mu.length())]+list(reversed(list(mu))) 
    return [i-v[i] for i in range(n)]

def to_ens(mu):
    nu=list(mu)+[0]
    return Set([k for k in range(1,nu.__len__()) if nu[k-1]>nu[k]])

def to_ideal(nu,n):
    nu=Partition(nu)
    w=to_dyckword(nu,n)
    return [to_partition(mu) for mu in DyckWord(bottom(n)).tamari_interval(w).dyck_words()]
   
def to_dyckword(mu,n):
    return DyckWord(list(DyckWords().from_area_sequence(to_area(mu,n))))

def tobinom(pol):
    pol=numerator(pol)
    return add(pol.coefficient({q:k})*binqn(k) for k in range(pol.degree()+1))

def ToBin(pol):
    pol=numerator(pol)
    return add(a*CC('k',b.degree()) for a,b in pol)

def to_partition(w):
    c=w.to_area_sequence()
    n=c.__len__()
    return Partition([n-1-i-c[n-1-i] for i in range(n)])

@cached_function
def trigraded(r,n): 
    return s(add((-1)**(n-ell(mu))*(r*n+1)**(ell(mu)-2)
                 *mul(binomial((r+1)*k,k) for k in mu)*p(mu)/zee(mu) for mu in Partitions(n)))

@cached_function
def trigraded_prime(r,n): 
    return s(add((-1)**(n-ell(mu))*((r+1)*n-1)**(ell(mu)-2)
                 *mul(binomial((r+1)*k-1,k) for k in mu)*p(mu)/zee(mu) for mu in Partitions(n)))


########  LLT a revoir

def CLLT(mu,m,n,u=None):
    if u==None: return CLLT(mu,m,n,t) 
    else: return add(t^dinv(pf,m,n)*SchurComp(comp(pf,m,n)) for pf in Park(mu,n))

def dinv(pf,m,n):
    mu=forme(pf)
    return tdinv(pf,m,n)+pdinv(mu,m,n)-maxtdinv(mu,m)
    
def pdinv(mu,m,n):
    return add(mu.arm_length(i,j)/(mu.leg_length(i,j)+1)<=m/n and
               mu.leg_length(i,j)*m/n<(mu.arm_length(i,j)+1) for i,j in mu.cells())

def rang(d,m,n):
    return m*d[1]-n*d[0]+floor(d[0]*gcd(m,n)/m)

def cellule(pf,i):
    c=(pf.cells_containing(i))[0]
    return (c[1],c[0])

def cellules(pf):
    return tuple((c[1],c[0]) for c in pf.cells())

def tdinv(pf,m,n):
    return add(rang(fr(cellule(pf,i),n),m,n)<rang(fr(cellule(pf,j),n),m,n)
               and rang(fr(cellule(pf,j),n),m,n)<rang(fr(cellule(pf,i),n),m,n)+m 
               for i in range(1,n) for j in range(i,n+1))

def perm(pf,m,n):
    r=[rang(fr(c,n),m,n) for c in cellules(pf)]
    return [position(rang(fr(cellule(pf,i),n),m,n),r) for i in range(1,n+1)]

def inverse(perm):
    n=perm.__len__()
    return [position(i+1,perm) for i in range(n)]

def desc(L):
    n=L.__len__()
    return sorted([i+1 for i in range(n-1) if L[i]>L[i+1]])

def comp(pf,m,n):
    S=desc(inverse(perm(pf,m,n)))+[n]
    return Composition([S[0]]+[S[i+1]-S[i] for i in range(S.__len__()-1)])


def position(a,L):
    for i in range(L.__len__()):
        if L[i]==a: return i+1

def fr(c,n):
    return tuple((c[0],n-1-c[1]))

def maxtdinv(mu,m):
    return max([tdinv(pf,m,n) for pf in Park(mu,n)])

def forme(pf):
    return Partition([c[0] for c in cellules(pf)])

def SchurComp(L):
    n=L.__len__()
    def hh(k):
        if k<0: return 0
        else: return h[k]
    return s(matrix([[hh(L[j]+i-j) for i in range(n)] for j in range(n)]).determinant())

