from functools import wraps # This convenience func preserves name and docstring
zero=Partition([0])
ClassePartages=zero.__class__


from sage.combinat.diagram import Diagram


def add_method(cls):
    def decorator(func):
        @wraps(func) 
        def wrapper(self, *args, **kwargs): 
            return func(self,*args, **kwargs)
        setattr(cls, func.__name__, wrapper)
        # return func # returning func means func can still be used normally
    return decorator

def Cellules(mu):
    mu=Partition(mu)
    return [(b+1,a+1) for a in range(mu.length()) for b in range(mu[a])]

def grid(k,n):
    return (add(line([(0,i),(k,i)],color='lightgrey') for i in range(n+1))
            +add(line([(i,0),(i,n)],color='lightgrey') for i in range(k+1)))


def carre(c,col='green'):
    e=0.03
    a,b=c
    return polygon([(a-1+e,b-1+e),(a-1+e,b-e),(a-e,b-e),(a-e,b-1+e)],
                   #edgecolor="black"
                   color=col,axes=False,figsize=3,thickness=.5)

def Max(L):
    return max(L+[0])

def Part_to_num(mu):
    mu=Partition(mu)
    if mu==Partition([]):
        return 0
    else:
        n=mu.length()
        return add(10**(n-i-1)*mu[i] for i in range(n))

@add_method(Diagram([(0,0)]).__class__)
def diagram(self,col='green'):
    pp=plot(grid(self.number_of_cols(),self.number_of_rows())
             +add(carre((c[1]+1,c[0]+1),col=col) for c in self.cells()))
    pp.axes(show=False)
    return pp

@add_method(ClassePartages)
def AjouteColonne(self,k):
    if k>=self.length():
        return Partition([k]+list(self.conjugate())).conjugate()
    else:
        return None

@add_method(ClassePartages)
def AjouteLigne(self,k):
    if k>=0:
        return Partition(sorted([k]+list(self),reverse=true))


    

@add_method(ClassePartages)
def B(self):
    return add(q**(i-Integer(1))*t**(j-Integer(1)) for i,j in self.Cells())

@add_method(ClassePartages)
def eta(self):
    return mul(q**(i-Integer(1))*t**(j-Integer(1)) for i,j in self.Cells())



#@add_method(ClassePartages)
#def Bar(self):
#    return Partition(self[1:])

@add_method(ClassePartages)
def Cells(self):
    return Cellules(self)

@add_method(ClassePartages)
def Contained(self):
    return [rho for k in range(self.size()+1) for rho in Partitions(k,outer=self)]


@add_method(ClassePartages)
def Dinv(self,pente,epsilon=1/1000):
    return Set([c for c in self.Cells() 
                if leg(c,self)/(arm(c,self)+1) < pente+epsilon < div(leg(c,self)+1,arm(c,self))])

@add_method(ClassePartages)
def dinv(self,pente,epsilon=1/1000):
    return self.Dinv(pente).cardinality()

@add_method(ClassePartages)
def diagram(self):
    return Diagram(self.cells()).diagram()


@add_method(ClassePartages)
def height(self):
    return self.length()

@add_method(ClassePartages)
def is_dominant(nu):
    nu=Partition(nu)
    return nu.dominates(nu.conjugate())

@add_method(ClassePartages)
def width(self):
    return self.conjugate().length()


@add_method(ClassePartages)
def Moins(self,nu):
    nu=Partition(nu)
    a=self.length()
    b=nu.length()
    try:
        return Partition([self[j]-nu[j] for j in range(b)]
                         +[self[j] for j in range(b,a)])
    except:
        return None

@add_method(ClassePartages)
def nabla_split(self):
    n=self.length()
    rho=staircase(n+1,n+1)
    d=min(self[i]//rho[i] for i in range(n))
    return (d,self.Moins([d*rho[i] for i in range(n)]))


@add_method(ClassePartages)
def Plus(self,nu):
    mu=self
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
    
@add_method(ClassePartages)    
def Pred(self,strict=True):
    if strict:
        first=[]
    else:
        first=[self]
    if self==zero:
        return first
    else:
        L=sorted(Diagonale(self),reverse=true)
        k=L.__len__()
        if Is_Dominant(self):
            return list(reversed(first+[self.RemoveCells(L[0:j]) for j in range(1,k+1)]))
        else:
            return list(reversed(first+[self.RemoveCells(L[j:]) for j in range(0,k)]))
    
@add_method(ClassePartages)
def RemoveCells(self,L):
    if L.__len__()==0:
        return self
    else:
        return enleve(L[0],self).RemoveCells(L[1:])


@add_method(ClassePartages)    
def risers(self):
    nu=list(self)
    return Composition(tuple(nu.count(i) for i in Set(nu)))

@add_method(ClassePartages)    
def contre_marches(self,n=None):
    if n==None:
        n=self.length()+1
    nu=[0 for i in range(n-self.length())]+list(self)
    return Composition(tuple(nu.count(i) for i in range(Max(nu)+1)))


    
@add_method(ClassePartages)
def SchurHatNabla(self,k=0):
    return Formal_Symmetric(1/(q*t)**k*Schur_hat(self).nabla())


    
