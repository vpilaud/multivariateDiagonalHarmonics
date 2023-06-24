from functools import wraps # This convenience func preserves name and docstring
Zero=tensor([0*s[0],0*s[0]])
ClasseTenseurs=Zero.__class__

def Print_tensor(F,tau): 
    Liste=Set([nu for (mu,nu),c in F.To_b(b=s)])
    print("Partition(",tau,"):(")
    for mu in sorted(Liste):
        val=F.scalar(s(mu))
        if val==Un:
            print("+tensor([s([0]),s(",mu,")])")
        elif not val==Integer(Integer(Integer(0))):
            print("+tensor([",val,",s(",mu,")])")
    print("),")

def T_Schur_hat(nu):
    return tensor([(-Un/(q*t))**iota(nu),s(nu)])

def Partie_Positive(g):
    if g==0:
        return 0*Un
    else:
        return add(c*s(mu) for mu,c in s(g) if c>0)
    
def Support(g):
    return [nu for nu,c in g]


def hook(n,k):
    return Partition([n-k]+[1 for i in range(k)])

def To_s_tensor(F):
    return tensor([s for i in range(F.length())])(F)

def Is_Sym(F):
    return Sym.is_parent_of(F)




def add_method(cls):
    def decorator(func):
        @wraps(func) 
        def wrapper(self, *args, **kwargs): 
            return func(self,*args, **kwargs)
        setattr(cls, func.__name__, wrapper)
        # return func # returning func means func can still be used normally
    return decorator

def Dishout(F,g=None,Famille=None):
    if F==0: 
        show(0)
    elif g==None:
        Dishout(F,g=s,Famille=Famille)
    else:
        n=F.degree()
        if n==0:
            show(F)
        else:
            if Famille==None:
                Liste=sorted(Partitions(n))
            else:
                Liste=sorted(Famille)
            for mu in Liste:
                First=Liste[0]
                gdual=g.dual_basis()
                val=F.scalar(gdual(mu))
                if mu==First and val.__len__()==1 and not val==0:
                    show(val,LatexExpr("\\otimes"),g(mu))
                elif mu==First and not val==0:
                    show(LatexExpr("\\big("),val,LatexExpr("\\big)"),LatexExpr("\\otimes"),g(mu))
                elif val.__len__()==1 and not val==0:
                    show(LatexExpr("\\qquad +\\ "),val,LatexExpr("\\otimes"),g(mu))
                elif not val==0:
                    show(LatexExpr("\\qquad +\\ "),LatexExpr("\\big("),val,LatexExpr("\\big)")
                         ,LatexExpr("\\otimes"),g(mu))

def Bar(nu):
    if nu.length()==1:
        return zero
    else:
        return Partition(nu[1:])
    

    
@add_method(ClasseTenseurs)
def bar(self,composante=-1,base=e):
    self=tensor([s,base])(self)
    return add(c*tensor([s(mu) for mu in v[:-1]]+[base(Bar(v[-1]))]) for v,c in self)

@add_method(ClasseTenseurs)
def AddCol(self,k,composante=1):
    self=tensor([s,s])(self)
    return self.modify(composante=composante,base=s,transformation=lambda nu:nu.AjouteColonne(k-nu.size()))

@add_method(ClasseTenseurs)
def AddRow(self,k):
    self=tensor([s,e])(self)
    return self.modify(composante=1,
                       base=e,
                       transformation=(lambda nu:nu.AjouteLigne(k-nu.size())))


@add_method(ClasseTenseurs)
def degree(self,composante=1):
    return max([0]+[v[composante].size() for v,c in self])

@add_method(ClasseTenseurs)
def length(self):
    return max([0]+[len(v) for v,c in self])


@add_method(ClasseTenseurs)
def Down(self,composante=1):
    def Sub1(nu):
        return Partition([nu[0]-1]+nu[1:])
    return self.modify(composante=composante,
                    condition=lambda nu:nu.length()==1 or nu[0]>nu[1],
                    transformation=Sub1)




@add_method(ClasseTenseurs)
def etend(self,k):
    if k<0:
        return tensor([0*Un,0*Un])
    elif k==0:
        return self
    else:
        return self.modify(transformation=(lambda nu:list(nu)+[1 for i in range(k)]))

@add_method(ClasseTenseurs)
def Formal_Coefficients(self):
    return add(tensor([InSchur(c)]+[s(nu) for nu in v]) for v,c in To_s_tensor(self))



@add_method(ClasseTenseurs)
def Hooks_Hooks_to_Pol(self):
    G=self.modify(composante=1,condition=Is_Hook)
    return G.modify(composante=0,transformation=lambda nu:Hooks_to_Pol(s(nu))).Eval(1)

@add_method(ClasseTenseurs)
def Hooks_Hooks(self):
    G=self.modify(composante=1,condition=Is_Hook)
    return G.modify(composante=0,condition=Is_Hook)


@add_method(ClasseTenseurs)
def is_zero(self):
    return self==tensor([0*Un,0*Un])

@add_method(ClasseTenseurs)
def is_positive(self,base0=s,base1=s):
    if self==Zero:
        return True
    else:
        return min((tensor([base0,base1])(self)).coefficients())>=0

@add_method(ClasseTenseurs)
def nabla(self,puissance=1):
    n=self.degree()
    return add(add(tensor([Produit_Formel(s(mu),c*s(nu).nabla(power=puissance).scalar(s(rho))),s(rho)]) 
                   for (mu,nu),c in self) for rho in Partitions(n))        

@add_method(ClasseTenseurs)
def map_fonction(self,fonction,composante=0):
    if composante==1:
        return add(tensor([s(mu),fonction(s(nu))]) for (mu,nu),c in tensor([s,s])(self))
    elif composante==0:
        n=self.degree()
        self=tensor([s,s])(self)
        return add(tensor([fonction(self.scalar(s(nu))),s(nu)]) for nu in Partitions(n))
    else:
        self


@add_method(ClasseTenseurs)
def modify(self,composante=1,transformation=(lambda nu:nu),condition=(lambda nu:True),base=s):
    if composante==1:
        return add(c*tensor([s(mu),base(transformation(nu))]) 
                   for (mu,nu),c in tensor([s,base])(self) if condition(nu))
    elif composante==0:
        return add(c*tensor([base(transformation(mu)),s(nu)]) 
                   for (mu,nu),c in tensor([base,s])(self) if condition(mu))
    else:
        self

@add_method(ClasseTenseurs)
def omega(self,k=0):
    return add(c*tensor([s(mu) for mu in v[:k]] + [s(v[k]).omega()] + [s(mu) for mu in v[k+1:]]) for v,c in To_s_tensor(self))

        
@add_method(ClasseTenseurs)
def omega_all(self):
    return add(c*tensor([s(nu).omega() for nu in v]) for v,c in To_s_tensor(self))

        
@add_method(ClasseTenseurs)
def PP(self,base=s,Table=False):
    if self==0:
        return self
    elif self.length()<=2 and not(Table):
        Liste=Set([v[0] for v,c in To_s_tensor(self)])
        return Dishout(self,g=base,Famille=Liste)
    else:
        return TableCoeff(self)


@add_method(ClasseTenseurs)
def PPse(self):
    n=self.degree()
    Famille=Tous_Partages(n)
    First=true
    for mu in sorted(Famille):
        val=self.scalar(f(mu))
        if First and val.__len__()==1 and not val==0:
            show(val,LatexExpr("\\otimes"),e(mu))
            First=false
        elif First and not val==0:
            show(LatexExpr("\\big("),val,LatexExpr("\\big)"),LatexExpr("\\otimes"),e(mu))
            First=false
        elif val.__len__()==1 and not val==0:
            show(LatexExpr("\\qquad +\\ "),val,LatexExpr("\\otimes"),e(mu))
        elif not val==0:
            show(LatexExpr("\\qquad +\\ "),LatexExpr("\\big("),
                 val,LatexExpr("\\big)"),LatexExpr("\\otimes"),e(mu))

@add_method(ClasseTenseurs)
def Pleth(self,g,composante=0,exclude=None):
    k=composante
    return add(c*tensor([s(mu) for mu in v[:k]] 
                        + [s(v[k])(g,exclude=exclude)] 
                        + [s(mu) for mu in v[k+1:]]) 
               for v,c in To_s_tensor(self))


@add_method(ClasseTenseurs)
def Plus(self,k):
    return self.modify(composante=0,transformation=(lambda nu:nu.Plus(Partition([k]).conjugate())))

@add_method(ClasseTenseurs)
def restrict_length(self,k,exact=false,composante=0):
    if exact:
        return add(c*tensor([s(v[0]),s(v[1])]) for v,c in tensor([s,s])(self) 
                   if v[composante].length()==k)
    else:    
        return add(c*tensor([s(v[0]),s(v[1])]) for v,c in tensor([s,s])(self) 
                   if v[composante].length()<=k)


@add_method(ClasseTenseurs)
def joint_restrict_degree(self,K,L={0,1}):
    return add(c*tensor([s(w) for w in v]) for v,c in (self).To_b(b=s) 
                   if add(v[i].size() for i in L)<=K)

@add_method(ClasseTenseurs)
def restrict_support(self,Famille=None,Test=(lambda nu:True),base=s,composante=1):
    if Famille==None:
        Famille=Set([nu for (mu,nu),c in tensor([s,base])(self)])
    if composante==1:    
        return add(c*tensor([s(mu),base(nu)]) for (mu,nu),c in tensor([s,base])(self) 
                   if nu in Famille and Test(nu))
    elif composante==0:    
        return add(c*tensor([s(mu),base(nu)]) for (mu,nu),c in tensor([s,base])(self) 
                   if nu in Famille and Test(mu))


@add_method(ClasseTenseurs)
def skew_by(self,g,k=0):
    if self==0:
        return 0
    elif k==0:
        return add(c*tensor([s(v[0]).skew_by(g)]+[s(nu) for nu in v[1:]]) 
                   for v,c in To_s_tensor(self))
    else:
        return add(c*tensor([s(nu) for nu in v[0:k]] + [s[v[k]].skew_by(g)]+ [s(mu) for mu in v[k+1:]]) 
                   for v,c in To_s_tensor(self))



@add_method(ClasseTenseurs)
def Eval(self,g,composante=0,exclus=None):
    if self==0:
        return self
    else:
        k=composante
        rep=add(c*s(v[k])(g*Un,exclude=exclus).scalar(Un)*tensor([s(mu) for mu in v[:k]] 
                        + [s(mu) for mu in v[k+1:]]) 
               for v,c in To_s_tensor(self))
        if rep.length()==1:
            return add(c*s(v[0]) for v,c in rep)
        else:
            return rep
    
@add_method(ClasseTenseurs)
def pleth(self,g,k=0,exclus={}):
    if k>=self.length():
        return self
    elif k==0:
        return add(c*tensor([s[v[0]](g*Un,exclude=exclus)]+ [s(mu) for mu in v[1:]]) for v,c in To_s_tensor(self))
    else:
        return add(c*tensor([s(nu) for nu in v[0:k]] + [s[v[k]](g*Un,exclude=exclus)]+ [s(mu) for mu in v[k+1:]]) 
                   for v,c in To_s_tensor(self)) 


@add_method(ClasseTenseurs)
def scalar(self,g,k=-1):
    if self==0 or g==0:
        return 0
    elif self.length()==2 and k==-1:
        return add(c*(s(nu).scalar(g))*s(mu) for (mu,nu),c in To_s_tensor(self))
    elif self.length()==2:
        return add(c*(s(mu).scalar(g))*s(nu) for (mu,nu),c in To_s_tensor(self))
    elif k==-1:
        return add(c*(s(v[-1]).scalar(g))*tensor([s(nu) for nu in v[:-1]]) for v,c in To_s_tensor(self))
    else:
        return add(c*(s[v[k]].scalar(g))*tensor([s(nu) for nu in v[0:k]] + [s(mu) for mu in v[k+1:]]) 
                   for v,c in To_s_tensor(self)) 


@add_method(ClasseTenseurs)
def Scalar(self,F):
    if self==0 or F==0:
        return 0
    else: 
        return add(a1*a2*(v==w) for v,a1 in self.To_b(b=s) for w,a2 in F.To_b(b=s))

@add_method(ClasseTenseurs)
def Scalar(self,F,b=s,bdual=s):
    if self==0 or F==0:
        return 0
    else: 
        return add(a1*a2*(v==w) for v,a1 in self.To_b(b=b) for w,a2 in F.To_b(b=bdual))


@add_method(ClasseTenseurs)
def Skew_close_Alt(self):
    n=self.degree()
    A=self.scalar(e[n])
    for k in range(n-3):
        B=A.skew_by(e[n-1-k])
        A=A+add(max(0,c)*s(mu.AjouteColonne(n-1-k)) for mu,c in (self.scalar(s(hook(n,k)))-B))
    return A



@add_method(ClasseTenseurs)
def Skew_Property(self):
    n=self.degree()
    A=self.Skew_close_Alt()
    return self.modify(condition=lambda nu:not Is_Hook(nu))+add(tensor([A.skew_by(e[n-1-k]),s(hook(n,k))]) for k in range(n))

@add_method(ClasseTenseurs)
def SYM(self,b=s):
    return add(c*tensor([b(nu) for nu in v]) for v,c in self.To_b(b) if list(v)==sorted(v))

@add_method(ClasseTenseurs)
def To_b(self,b=s):
    return tensor([b for i in range(self.length())])(self)
    
@add_method(ClasseTenseurs)
def UpFirst(self,k,j=1):
    if j==0:
        return self
    elif j==1:
        return self.modify(composante=0,transformation=lambda nu:nu.AjouteColonne(k))
    else:
        return self.UpFirst(k).UpFirst(k,j-1)
    