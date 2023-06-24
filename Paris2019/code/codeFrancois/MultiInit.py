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

def catalan(n): return (1/(n+1)*binomial(2*n,n))

def CC(k,j): return SR(k).binomial(j, hold=True)

def Cellules(c):
    m=c.__len__()
    return [(a,b) for a in range(m) for b in range(c[a])]

def Coefficient(pol,t,k):
    return pol.numerator().coefficient({t:k})

def CoInversions(p):
    n=p.__len__()
    return add(1 for i in range(n-1) for j in range(i+1,n) if p[i]<p[j])

def C_r_mu(r,mu):
    return mul(1/(r*k+1)*binomial((r+1)*k,k) for k in mu)

def decalage(m,n):
    return (add(vector(c) for c in Cellules(FloorList(m,n)))
            -add(vector(c) for c in Cellules(Partlist(m,n))))[0]

def Defaut(m,n):
    return CoInversions(FloorList(m,n))

def DegAlt(a,b): 
    return (a*b-a-b+gcd(a,b))/2

@cached_function
def degmn(a,b):
    a=Integer(a)
    b=Integer(b)
    return add(floor(i*(a/b)) for i in range(b))

def DimTam(m,n): return (m+1)**n*(m*n+1)**(n-2)

def DimAltTam(m,n): return (m+1)/(n*(m*n+1))*binomial((m+1)**2*n+m,n-1)

def Factorise(F):
    if F==0: return 0
    elif F==1: return 1
    else: 
        return factor(F)
    
def flip(mu,k,n):
    mu=Partition(mu)
    l=mu.length()
    return [k for j in range(n-1-mu.length())]+[k-mu[l-j-1] for j in range(l)]

def FloorList(m,n):
    return Composition(list(floor(i*m/n) for i in range(n)).count(k) for k in range(m))

def tobinom(pol):
    pol=numerator(pol)
    return add(pol.coefficient({q:k})*binqn(k) for k in range(pol.degree()+1))

def limite(m,n): return min(n-1-m,m)

def monfactor(pol):
    pol=numerator(pol)
    if pol.degree()==0: return pol
    else: return factor(pol)
    
def N_mu(c):
    return add(a for (a,b) in Cellules(c))

def Parkj(k,j,n):
    return [mu for mu in Partitions(j,parts_in={i for i in range(1,k+1)}) 
                if mu.length()<n]

def Partlist(m,n):
    return Partition(sorted(FloorList(m,n),reverse=true))

def polynome(g):
    return simplify(g).numerator()

def Positif(f):
    if f==0: return true
    elif min(f.coefficients()) > 0: return true
    else: return false

def rho(a,b):
    a=Integer(a)
    b=Integer(b)
    return add(floor(i*(a/b)) for i in range(b))

def ToBin(pol):
    pol=numerator(pol)
    return add(a*CC('k',b.degree()) for a,b in pol)

######################## Frobenius divers

@cached_function
def bigraded(m,n):
    return s(add((-1)**(n-ell(mu))*(m*n+1)**(ell(mu)-1)*p(mu)/zee(mu) for mu in Partitions(n)))

@cached_function
def Coinvariant(n,k=None): 
    r=2
    if k is None: k=1
    return add(Tenseur(pol_to_s(c),s(mu)) 
               for mu,c in s(mul((1-q**i) 
                                 for i in range(1,n+1))*h[n](X/(1-q))).map_coefficients(polynome))  
@cached_function
def trigraded(r,n): 
    return s(add((-1)**(n-ell(mu))*(r*n+1)**(ell(mu)-2)
                 *mul(binomial((r+1)*k,k) for k in mu)*p(mu)/zee(mu) for mu in Partitions(n)))


######################## Composantes et Multi-Frobenius

@cached_function
def a_mnz(m,n):
    return add(RH(a_mn(m,n,k))*z**k for k in range(min(m,n)))

@cached_function
def b_mnz(m,n):
    return add(RH(b_mn(m,n,k))*z**k for k in range(min(m,n)))

def b_mn(m,n,k):
    if n>m:
        return b_mn(n,m,k)
    else:
        return (Scalar2(e_mn(m,n),s(hook(n,n-1-k)))
            -Skew(add(addcol(j,a_mn(m,n,j)) 
                      for j in range(k+1,n)),e[k])).restrict_partition_lengths(k,exact=false)

def Alt_mn(m,n):
    if min(m,n)<=3:
        return Scalar2(E_mn(m,n),e[n])
    elif m==n+1:
        return Alt_mn(n,n)
    elif m==n-1:
        return Alt_mn(n-1,n-1)
    elif n==4:
        return (addmu(colone(n-1),Alt_mn(m-n,n))
                +Scalar2(E_mn(m,n),e[n]).restrict_partition_lengths(n-2,exact=false))
    elif m<=n:
        return add(Dk(k,m,n) for k in range(n))
    elif m>n:
        return Alt_mn(n,m)
    else:
        return 0*Un
    
def Approx(mu,m=None): 
    k=add(1 for j in mu)
    n=add(j for j in mu)
    if m==None: 
        Approx(mu,n)
    elif m==n: 
        return Cmu(colone(n),n).skew_by(s(Partition(mu)).skew_by(e[k]).omega())-second(mu)
    else: 
        return Cmu(colone(n),m).skew_by(s(Partition(mu)).skew_by(e[k]).omega())

def C(mu):
    n=add(k for k in mu)
    return Scalar2(Phi[n],s(mu))

@cached_function
def Cj(mu,m,j):
    return restreint(C_m(m,mu),j)

def Cmu(mu,m=None):
    n=add(j for j in mu)
    if m==None: 
        return Cmu(mu,n)
    elif mu==colone(n) and m==n-1: 
        return Cmu(colone(n-1),n-1)
    else: 
        return Scalar2(Fqt(m,n),s(mu))

def Cn_211(n): 
    if n<3: return 0
    else: return alpha(n)+e[1]*add(s[k] for k in range(2,n))+s[1,1]

def Cnab(mu):
    n=add(k for k in mu)
    return Scalar2(Nabla(n),s(mu))

@cached_function
def Coeff_mu(mu):
    n=add(j for j in mu)
    if max(mu)==1: return Phi_Alt(n)
    else:
        nu=No_One(mu)
        k=n-add(j for j in nu)
        return add(Phi_Alt(n-c).skew_by(g) for rho,[c,g] in Regles(n) if rho==nu)

def Coeff_mup(mu):
    n=add(j for j in mu)
    if max(mu)==1: return Formule_Altp(n)
    else:
        nu=No_One(mu)
        k=n-add(j for j in nu)
        return add(Formule_Altp(n-c).skew_by(g) for rho,[c,g] in Regles(n) 
                   if rho==nu)

def Correction(m,n):
    if m<n or n<4:
        return 0
    elif n==4 and m==4:
        return Tenseur(Alt_mn(3,4),s[2,2])
    elif n==5 and m==5:
        return (Tenseur(Alt_mn(4,5),s[2,2,1])
                +Tenseur(Alt_mn(4,5).skew_by(e[1]),s[3,2]))
    elif n==6 and m==6:
        return (Tenseur(Alt_mn(5,6),s[2,2,1,1])
                +Tenseur(Alt_mn(5,6).skew_by(e[2]),s[4,2])
                +Tenseur(Alt_mn(5,6).skew_by(e[1]),s[3,2,1])
                +Tenseur(Alt_mn(5,6).skew_by(s[2])-Alt_mn(4,5),s[3,3])
                +Tenseur(Alt_mn(6,6).skew_by(s[4]),s[2,2,2]))
    elif n==7 and m==7:
        return (Tenseur(Alt_mn(6,7),s[2,2,1,1,1])
                +Tenseur(Alt_mn(6,7).skew_by(e[3]),s[5,2])
                +Tenseur(Alt_mn(6,7).skew_by(e[2]),s[4,2,1])
                +Tenseur(Alt_mn(6,7).skew_by(e[1]),s[3,2,1,1])
                +Tenseur(Alt_mn(6,7).skew_by(s[2])-Alt_mn(5,6),s[3,3,1])
                +Tenseur(Alt_mn(6,7).skew_by(s[2,1])-Alt_mn(5,6).skew_by(s[1]),s[4,3]))
    elif m==n+1:
        return Correction(n,n)
    else: 
        return 0

def Cp(mu):
    if mu==[4,1,1,1]: return C4111
    else:
        n=add(k for k in mu)
        return Scalar2(Nabla(n),s(mu))
    
def C_From_Alt(mu,m):
    return Premier(mu,m)-Second(mu,m)

def C_m(m,mu):
    n=add(j for j in mu)
    return Scalar2(D_m(m,n),s(mu))  

def F1(n): 
    return tensor([Un,s(e[n])])

def F2(n): 
    return add(tensor([s[n//2].skew_by(h[i]),s(conjugate(Partition([n-i,i])))]) for i in range(n//2+1))

@cached_function
def Formule_Alt(n):
    if n==1: return(Un)
    else:
        return add(addmu(conjugate([k]),delta(k,n)) for k in range(1,n))

@cached_function
def Formule_Altp(n):
    if n==1: return(Un)
    else:
        return add(addcol(k,deltap(k,n)) for k in range(1,n))

@cached_function
def Formule_Phip(n):
    if n<=6: 
        return Phi[n]
    else: 
        return add(tensor([Coeff_mup(mu),s(mu)]) for mu in Partitions(n))

@cached_function
def Formule_Phi(n):
    return add(tensor([Coeff_mu(mu),s(mu)]) for mu in Partitions(n))

@cached_function
def Fqt(m,n):
    d=gcd(m,n)
    a=m//d
    b=n//d
    if n<=3 or m<=3: 
        return Theta(e[d],a,b)
    elif m==n and n<=6: 
        return Formule_Phip(n)
    elif m==4 and n==5: 
        return (add(tensor([C(mu),s(mu+[1])]) for mu in Partitions(4))
                +tensor([C([2,1]),s[3,2]])
                +tensor([C([1,1,1]),s[2,2,1]]))
    else: 
        return 0

def gen_dim_alt(n):
    w=numerator(tobinom(Phi_Alt(n)(q*Un,exclude=[q]).scalar(Un)))
    return add(a*CC('k',b.degree()) for a,b in w)

@cached_function
def gen_dim_altp(n):
    w=numerator(tobinom(Formule_Altp(n)(q*Un,exclude=[q]).scalar(Un)))
    return add(a*CC('k',b.degree()) for a,b in w)

def Gqt(m,n): 
    return add(tensor([Approx(mu,m),s(mu)]) for mu in Partitions(n))

@func_persist
def O_n(n):
    return add(Tenseur(((s_mu_h(mu)
                         *Inv_hh(n)).restrict_degree(n,exact=false)),s(mu))
               for mu in Partitions(n))

@cached_function
def OO_n(n):
    return restrict_degree1(
        restrict_degree1(s([n])(Tenseur(add(h[k] for k in range(n+1)),X)),n)
        *Tenseur(Inv_hh(n),Un),n)

@cached_function
def Phi4(m):
    if m<=4: 
        return D_m(m,4)
    else: 
        return UpArrow1(3,D_m(m-4,4))+E_mn(m,4)

def Phi_Alt(n):
    return C([1 for j in range(n)])    

@cached_function
def PhiMat(m,n):
    return Matrix([[Cj(mu,m,j) for j in range(0,n)] for mu in Partitions(n)])

@cached_function
def PhiMatSkew(m,n):
    return Matrix([[Skew(Cj(mu,m,j),e[j]) for j in range(0,n)] for mu in Partitions(n)])

def premier(mu): 
    k=add(1 for j in mu)
    n=add(j for j in mu)
    return Coeff_mup([1 for j in range(n)]).skew_by(s(Partition(mu)).skew_by(e[k]).omega())

def Premier(mu,m): 
    k=add(1 for j in mu)
    n=add(j for j in mu)
    return Alt_mn(m,n).skew_by(s(Partition(mu)).skew_by(e[k]).omega())

def second(mu):
    k=add(1 for j in mu)
    if k==2 and mu[1]>1: return Coeff_mup([mu[0]-1,mu[1]-1,1])
    elif k>=3 and mu[1]>1 and mu[2]==1: return Coeff_mup([mu[0]-1,mu[1]-1]+[1 for j in range(2,k+1)])
    else: return 0

def Second(mu,m): 
    n=add(j for j in mu)
    if m==n and not mu in hooks(n):
        return Cm(part_second(mu),n)
    else:
        return 0*Un

###############


def Psi(m,n): 
    return (eval_tensor_k(D_m(m,n),q).map_coefficients(tobinom)).map_coefficients(ToBin)

def Psie(m,n): 
    return (e(eval_tensor_k(Plus_un(D_m(m,n)),q).map_coefficients(tobinom))).map_coefficients(ToBin)

def DimPsi(m,n):
    return ToBin(tobinom(dim(eval_tensor_k(D_m(m,n),q))))

def eval_D_m(m,n,q): 
    return eval_tensor_k(D_m(m,n),q)

def Formule_cat(n):     
    return (-1)**(n-1)*add(Mult_part(mu)*catalan(ell(mu)-1)*f(mu) for mu in Partitions(n))

###############

def deltanab(m,n):
    return deltap(m,n).restrict_partition_lengths(2,exact=false)


@cached_function
def delta(m,n):
    if m==n-1: 
        return Un
    else:
        j=n-m-1
        return C(hook(n,j))-add(addcol(k,delta(k,n)).skew_by(e[m]) for k in range(m+1,n))

def Deltak(k,n): 
    if k>n-1: return 0
    else: return addcol(k,deltap(k,n))

@cached_function
def deltap(m,n):
    if m==n-1: 
        return Un
    else:
        j=n-m-1
        return (Cp(hook(n,j))
                -add(addcol(k,deltap(k,n)).skew_by(e[m]) 
                     for k in range(m+1,n))).restrict_partition_lengths(limite(m,n),exact=false)
def D(k,m,n):
    return restreint(Alt_mn(m,n),k)


def Dk(k,m,n):
    return addmu(colone(k),dk(k,m,n))

@cached_function
def dk(k,m,n):
    if m==n:
        return deltap(k,n)
    elif m==n-1:
        return deltap(k,n-1)
    elif k==0:
        return restreint(Scalar2(E_mn(m,n),e[n]),0)
    elif k<=2:
        return Skew(restreint(Scalar2(E_mn(m,n),e[n]),k),e[k])
    elif n-k==1 and m<n:
        return Scalar2(E_mn(m,n),s[n])
    elif n-k==2 and m<n:
        return (Scalar2(E_mn(m,n),s[n-1,1])
                -Skew(addmu(colone(n-1),dk(n-1,m,n)),e[n-2]))
    elif n-k==3 and m<n:
        return (Scalar2(E_mn(m,n),s[n-2,1,1])
                -Skew(addmu(colone(n-2),dk(n-2,m,n)),e[n-3])
                -Skew(addmu(colone(n-1),dk(n-1,m,n)),e[n-3]))
    else:
        return 0*Un

def Deltakmn(k,m,n):
    return restreint(Cmu(conjugate(Partition([n])),m),k)
    
def ksi(n): 
    if n<7: return 0 
    elif n % 3==0: return s[2*(n//3)-1,2*(n//3)-1]
    elif n % 3==2: return s[2*((n+1)//3)-1,2*((n+1)//3)-2]
    else: return 0                       

@cached_function
def alpha(n): 
    if n<7: return Cp(hook(n,2))-(e[1]*add(s[k] for k in range(2,n))+s[1,1])
    else: return (ksi(n)+alpha(n-1)+extend_first_row(1,alpha(n-1))+extend_first_row(2,alpha(n-1))
                 -extend_first_row(1,alpha(n-2))-extend_first_row(2,alpha(n-2))
                 -extend_first_row(3,alpha(n-2))+extend_first_row(3,alpha(n-3)))
def dkj(k,m,n,j):
    return restreint(dk(k,m,n),j)

def d(k,m,n):
    return C_m(hook(n,n-1-k),m)-add(updown(dk(a,m,n),a,k) for a in range(k+1,n))

@cached_function
def deltaMat(m,n):
    return Matrix([[dkj(k,m,n,j) for j in range(0,(n+1)//2)] for k in range(1,m)])

@cached_function
def deltaMatSkew(m,n):
    return Matrix([[Skew(dkj(k,m,n,j),e[j]) for j in range(0,(n+1)//2)] for k in range(1,m)])

@cached_function
def MatAlt(d,N):
    return Matrix([[Scalar(eval_tensor1(Delta_F(k,n),d),e[n]) for k in range(N+1)] 
                   for n in range(1,N+1)])
@cached_function
def MatDim(d,N):
    return Matrix([[Scalar(eval_tensor1(Delta_F(k,n),d),p[1]**n) for k in range(N+1)] 
                   for n in range(1,N+1)])


def E_equerres(m,n):
    return CoBar2(E_hooks(m,n))

def E_hooks(m,n):
    return add(Tenseur(Scalar2(E_mn(m,n),s(mu)),s(mu)) for mu in hooks(n))


