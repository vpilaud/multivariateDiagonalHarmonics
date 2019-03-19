def alpha_mn(k,m,n):
    return Scalar2(epsilon(k,m,n),e[n])

def Addmu1(rho,F):
    return add(c*Tenseur(addmu(rho,s(mu)),s(nu)) for (mu,nu),c in F)

def Addmu2(rho,F):
    return add(c*Tenseur(s(mu),addmu(rho,s(nu))) for (mu,nu),c in F)

def AddParts1(rho,F): 
    return add(c*Tenseur(s(list(mu)+list(rho)),s(nu)) for (mu,nu),c in F)

def AddParts2(rho,F): 
    return add(c*Tenseur(s(mu),s(list(nu)+list(rho))) for (mu,nu),c in F)

def Add_rows1(nu,F):
    return add(c*Tenseur(add_rows(nu,s(mu)),s(nu)) for (mu,nu),c in tensor([s,s])(F))
    
def Add_rows2(nu,F):
    n=Deg2(F)
    return add(Tenseur(Scalar2(F,s(mu)),add_rows(nu,s(mu))) for mu in Partitions(n))

def Alt(F): 
    def tensalt(a,b): return(a*alt(b))
    return F.apply_multilinear_morphism(tensalt)

def Bar1(F):
    if F==0: return 0
    else: return add(Tenseur(c*s(Bar(mu)),s(nu)) for (mu,nu),c in F)

def Bar2(F): 
    if F==0: return 0
    else: return add(c*Tenseur(s(mu),s(Bar(nu))) for (mu,nu),c in F)
    
def CatPart(rho,F):
    if F==0: return 0
    else: return add(c*Tenseur(s(mu),s(nu+rho)) for (mu,nu),c in F)
    
def CoBar1(F):
    if F==0: return 0
    else: return add(Tenseur(c*s(CoBar(mu)),s(nu)) for (mu,nu),c in F)

def CoBar2(F):
    if F==0: return 0
    else: return add(Tenseur(c*s(mu),s(CoBar(nu))) for (mu,nu),c in F)

def Conj1(F):
    return add(Tenseur(c*s(mu.conjugate()),s(nu)) for (mu,nu),c in tensor([s,s])(F))

def Conj2(F):
    return add(Tenseur(c*s(mu),s(nu.conjugate())) for (mu,nu),c in tensor([s,s])(F))
    
def Deg1(F):
    return max([nu.size() for (nu,mu),c in F])

def Deg2(F):
    return max([mu.size() for (nu,mu),c in F])

def Deg_Tens(F):
    return max(s(nu).degree() for (mu,nu),c in F)

def Dim(F): 
    def tensdim(a,b): return(a*dim(b))
    return F.apply_multilinear_morphism(tensdim)

def Eval1(F,k,E=None): 
    if F==0: return 0
    elif E==None:
        return Eval1(F,k,{})
    else:
        return add(c*(s(mu)(k*Un,exclude=E))*s(nu) for (mu,nu),c in tensor([s,s])(F))

def Eval2(F,k,E=None): 
    if F==0: return 0
    elif E==None:
        return Eval2(F,k,{})
    else:
        return add(c*s(mu)*(s(nu)(k*Un,exclude=E)) for (mu,nu),c in tensor([s,s])(F))
    
def Fleche1(rho,F): 
    if F==0:
        return 0
    else:
        return add(c*Tenseur(addmu(rho,s(mu)),s(nu)) for (mu,nu),c in F)

def Formal_F(F):
    n=F.degree()
    return add(tensor([isotypique(F,mu),s(mu)]) for mu in Partitions(n))

def length_component(F,k): 
    def comp_tens(a,b): return tensor([a.restrict_partition_lengths(k),b])
    return F.apply_multilinear_morphism(comp_tens)

def Longueur(F):
    if F==0: 
        return 0    
    else:
        return max(mu.length() for (mu,nu),c in F)
    
def Moins_un(F): 
    def tensor_moins_un(a,b): return tensor([a(X-Un),b])
    return tensor([s,e])(F.apply_multilinear_morphism(tensor_moins_un))

def Mult1(g,F):
    return add(Tenseur(c*s(mu)*g,s(nu)) for (mu,nu),c in tensor([s,s])(F))

def Mult2(g,F):
    return add(Tenseur(c*s(mu),s(nu)*g) for (mu,nu),c in tensor([s,s])(F))
    
def NoHooks(F):
    n=Deg2(F)
    return F-add(Tenseur(restrict_to_hooks(Scalar2(F,s(nu))),s(nu)) for nu in hooks(n))

def odot(F,G):
    return add(add(c*d*Tenseur(s(plus(mu,alpha)),s(concat(nu,beta))) for (mu,nu),c in tensor([s,s])(F)) 
               for (alpha,beta),d in tensor([s,s])(G))

def Omega1(F):
    n=Deg2(F)
    return add(c*Tenseur(s(mu.conjugate()),s(nu)) for (mu,nu),c in F)

def Omega2(F):
    n=Deg2(F)
    return add(Tenseur(Scalar2(F,s(mu)),s(mu.conjugate())) for mu in Partitions(n))

def Only_hooks(F):
    n=Deg2(F)
    return add(Tenseur(Scalar2(F,s(nu)),s(nu)) for nu in hooks(n))

@cached_function
def Plus_un(F): 
    def tensor_plus_un(a,b): return tensor([a(X+Un),b])
    return tensor([s,e])(F.apply_multilinear_morphism(tensor_plus_un))

def pleth_tensor1(F,k): 
    def tens1(a,b): return Tenseur(a(k*Un),b)
    return add(c*tens1(s(mu),s(nu)) for (mu,nu),c in tensor([s,s])(F))

def pleth_tensor2(F,k): 
    def tens2(a,b): return Tenseur(a,b(k*Un))
    return add(c*tens2(s(mu),s(nu)) for (mu,nu),c in tensor([s,s])(F))
    
def Pleth1(F,k,E=None): 
    if F==0:
        return 0
    elif E==None: 
        return Pleth1(F,k,{})
    else:
        return add(c*Tenseur(s(mu)(k*Un,exclude=E),s(nu)) for (mu,nu),c in tensor([s,s])(F))
    
def Pleth2(F,k,E=None): 
    if F==0:
        return 0
    elif E==None: 
        return Pleth2(F,k,{})
    else:
        return add(c*Tenseur(s(mu),s(nu)(k*Un,exclude=E)) for (mu,nu),c in tensor([s,s])(F))

def restreint_largeur1(F,k):
    if F==0:
        return F
    else:
        return add(c*Tenseur(s(mu),s(nu)) for (mu,nu),c in Toss(F) if mu==[] or mu[0]<=k)

def restreint_largeur2(F,k):
    if F==0:
        return F
    else:
        return add(c*Tenseur(s(mu),s(nu)) for (mu,nu),c in Toss(F) if nu==[] or nu[0]<=k)

def restrict_degree1(F,k): 
    return add(c*Tenseur(s(mu),s(nu)) for (mu,nu),c in tensor([s,s])(F) if mu.size()<=k)

def restrict_length1(F,k): 
    if F==0: 
        return 0
    else:
        def restr_tens1(a,b): return tensor([a.restrict_partition_lengths(k,exact=false),b])
        return F.apply_multilinear_morphism(restr_tens1)

def restrict_length2(F,k): 
    if F==0: 
        return 0
    else:
        def restr_tens2(a,b): return tensor([a,b.restrict_partition_lengths(k,exact=false)])
        return F.apply_multilinear_morphism(restr_tens2)
    
def RestrictHooks(F):
    n=Deg2(F)
    return add(Tenseur(restrict_to_hooks(Scalar2(F,s(nu))),s(nu)) for nu in hooks(n))

def Scalar1(U,f): 
    res=add(c*(s(a).scalar(f))*s(b) for (a,b),c in tensor([s,s])(U))
    if res==None: 
        return 0*Un
    else:
        return res

@cached_function
def Scalar2(U,f): 
    res=add(_c*s(_a)*(s(_b).scalar(f)) for (_a,_b),_c in tensor([s,s])(U))
    if res==None: 
        return 0*Un
    else:
        return res

def Skew1(g,F):
    return add(Tenseur(c*s(mu).skew_by(g),s(nu)) for (mu,nu),c in tensor([s,s])(F))

def Skew2(g,F):
    return add(Tenseur(c*s(mu),s(nu).skew_by(g)) for (mu,nu),c in tensor([s,s])(F))
    
def Strip2(F): 
    if F==0: return 0
    else: return add(c*Tenseur(s(mu),s(strip(nu))) for (mu,nu),c in F)

def Tenseur(f,g):
    if f==0 or g==0: return 0
    else: return tensor([f,g])
        
def Tensor_Add_mu(rho,F):
    return add(c*Tenseur(addmu(rho,s(nu)),s(mu)) for (nu,mu),c in F)

def Toss(F):
    return tensor([s,s])(F)

def Trivial(n):
    return Tenseur(Un,s[n])
    
def UpArrow1(k,F): 
    if F==0: return 0
    else: return add(c*Tenseur(addmu(colone(k),s(mu)),s(nu)) for (mu,nu),c in F)
        
def UpArrow2(k,F): 
    if F==0: return 0
    else: return add(c*Tenseur(s(mu),addmu(colone(k),s(nu))) for (mu,nu),c in F)

def UpPart1(rho,F): 
    if F==0: return 0
    else: return add(c*Tenseur(addmu(rho,s(mu)),s(nu)) for (mu,nu),c in F)
    
def width_component(F,k):
    if k==0:
        return add(c*Tenseur(s(mu),s(nu)) for (mu,nu),c in Toss(F) if mu==[])
    else:
        return add(c*Tenseur(s(mu),s(nu)) for (mu,nu),c in Toss(F) if mu<>[] and mu[0]==k)

def width_restrict1(F,k):
    if k==0:
        return add(c*Tenseur(s(mu),s(nu)) for (mu,nu),c in Toss(F) if mu==[])
    else:
        return add(c*Tenseur(s(mu),s(nu)) for (mu,nu),c in Toss(F) if mu==[] or mu[0]<=k)