########################  FONCTIONS SUR LES PARTAGES
def compact(mu):
    if mu==Partition([]): return 0
    else: return (''.join(mystr(i) for i in mu))

CL=['blue', 'gold', 'green','red', 'orange', 'purple', 'yellow','darkblue', 'darkcyan', 'darkorange', 'darkred', 'darkturquoise', 'darkviolet', 
           'deepskyblue', 'firebrick','indigo',
           'magenta',   'violet',  'yellowgreen']

Partitions.options.convention="french"
Partition._latex_= compact
Composition._latex_= compact

def AntiBar2(F,n): 
    if F==0: 
        return 0
    else: 
        return add(c*Tenseur(s(mu),s(plus(nu,colone(n-nu.size())))) for (mu,nu),c in F)
    
def Bar(mu):
    if mu.size()==0: 
        return mu
    else: 
        return Partition([c-1 for c in mu])

def buthooks(n):
    return sorted(list(Set(Partitions(n)).difference(Set(hooks(n)))))

def B_mu(mu):
    mu=Partition(mu)
    return sum((t**i*q**j) for i,j in mu.cells())

def B_muq(mu):
    mu=Partition(mu)
    return sum(q**(j-i) for i,j in mu.cells())

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

def CoBar(mu):
    if mu.size()==0: 
        return mu
    else:
        return Bar(mu.conjugate()).conjugate()

def concat(mu,nu):
    return Partition(sorted(list(mu)+list(nu),reverse=true))

def colone(n):
    return conjugate(Partition([n]))

def conjugate(mu):
    mu=Partition(mu)
    return mu.conjugate()

def conj_skew(mu):
    return Partition([k-1 for k in mu if k>1]).conjugate()

def dessine_chaine(c,m=None):
    if m==None:
        n=c.__len__()
        return dessine_chaine(c,n+add(c[i][0] for i in range(n)))
    else:
        n=c.__len__()
        def undiag(k):
            return Ferrers(c[k],CL[k],k+add(c[i][0] for i in range(k)))
        show(Ferrers(Partition([]),'white',m)
                                   +add(undiag(n-1-k) for k in range(n)))
        
def ell(mu): 
    mu=Partition(mu)
    return mu.length()

def fathooks(n):
    return Set(mu for mu in buthooks(n) if mu.length()==2 or mu[2]==1)


def Ferrers(mu,couleur='blue',decalage=0):
    mu=Partition(mu)
    m=partition_to_matrix(mu,decalage)
    return matrix_plot(m, cmap=['white',couleur],
            subdivisions=True, 
            subdivision_style=dict(color='white',thickness=3),
            axes=False,
            frame=False,
            origin='lower')

def hook(n,k):
    return [n-k]+[1 for i in range(k)]

def hooks(n):
    return sorted([Partition(hook(n,k)) for k in range(n)])

def iota(mu):
    mu=Partition(mu)
    return add(mu[k]-k-1 for k in range(mu.length()) if mu[k]>k)

def is_hook(mu):
    mu=Partition(mu)
    if mu.length()-1==mu.size()-mu[0]:
        return true
    else:
        return false
    
def K(mu,nu):
    mu=Partition(mu)
    nu=Partition(nu)
    return SemistandardTableaux(mu,nu).cardinality()

def Kostka(mu,nu):
    return s(mu).scalar(h(nu))

def LLT(mu,n,u=None):
    if u==None:
        return LLT(mu,n,t)
    else:
        if n>=7: correction=1
        else: correction=0
        mu=Partition(mu)
        if mu==Partition([0]) and n==1: return e[1]
        else:
            for (nu,llt) in ListeLLT[n-1]:
                if Partition(nu)==mu:
                    return e(llt.map_coefficients(lambda c: c.substitute({t:u-correction})))
                
def mystr(i): 
    if i<10: 
        return str(i) 
    else: 
        return ''.join([str(i),"."])
    
def Mult_part(mu): return prod(k for k in mu)

def n_mu(mu):
    mu=Partition(mu)
    return add(vector(c) for c in mu.cells())[1]

def No_One(mu):
    return Partition([k for k in mu if k<>1])

def Par(n): 
    return Partitions(n).list()

def Part_Euclid(a,b):
    return Partition([a for k in range(b//a)]+[b-a*(b//a)])

def part_premier(mu): 
    nu=Partition(mu).conjugate()
    return Partition([nu[k] for k in range(1,ell(nu))])

def part_second(mu):
    k=add(1 for j in mu)
    if k==2 and mu[1]>1: 
        return Partition([mu[0]-1,mu[1]-1,1])
    elif k>=3 and mu[1]>1 and mu[2]==1: 
        return Partition([mu[0]-1,mu[1]-1]+[1 for j in range(2,k+1)])
    else: return Partition([])

def partition_to_matrix(mu,decalage=0):
    mu=Partition(mu)
    if mu==[]:
        return matrix([])
    n=mu[0]
    m=matrix([[0 for i in range(decalage)]
              +[1 for i in range(mu[j])]
              +[0 for i in range(n-mu[j])] 
              for j in range(mu.length())])
    m.subdivide(range(1,mu.length()),range(1,n+decalage))
    return m
    
def Parts_un(mu): return add(k for k in mu if k==1)

def plus(mu,nu):
    mu=Partition(mu)
    nu=Partition(nu)
    a=mu.length()
    b=nu.length()
    if a<=b: 
        return Partition([mu[j]+nu[j]
                          for j in range(a)]+[nu[j] 
                                              for j in range(a,b)])
    else:
        return Partition([mu[j]+nu[j] 
                          for j in range(b)]+[mu[j]
                                              for j in range(b,a)])

def minus(mu,nu):
    mu=Partition(mu)
    nu=Partition(nu)
    a=mu.length()
    b=nu.length()
    if mu.contains(nu): 
        return [mu[j]-nu[j] 
                          for j in range(b)]+[mu[j]
                                              for j in range(b,a)]
    else:
        return []


def q_kostka(mu,nu):
    return add(q^(t.charge()) for t in SemistandardTableaux(mu,nu))

def rotate(mu,m):
    mu=Partition(mu)
    return Partition([m-k for k in range(1,m) for i in range(list(mu).count(k))])

def rest(mu):
    if mu.size()==0:
        return mu
    elif mu.length()==1:
        return Partition([])
    else:
        return Partition([mu[k] for k in range(1,mu.length())])  
    
def staircase(m,n):
    return Partition([(m*(n-k))//n for k in range(1,n+1)])

def strip(mu):
    if mu==[]: return Partition([0])
    else: return Partition([mu[j] for j in range(1,Partition(mu).length())])
    
def tonum(mu):
    Partition(mu).length()
    return add(10 (k-j-1)*mu[j] for j in range(k))

def zee(mu): 
    mu=Partition(mu)
    return mu.aut()