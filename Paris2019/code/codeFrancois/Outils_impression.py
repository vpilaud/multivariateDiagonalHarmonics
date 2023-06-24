def Dishout(F,g=None):
    if F==0: 
        show()
    elif g==None:
        Dishout(F,s)
    else:
        n=Deg2(F)
        for mu in sorted(list(Partitions(n))):
            val=Scalar2(F,g(mu))
            if val<>0:
                show((mu,val))

def DishoutLatex(F,g=None):
    if F==0: 
        show()
    elif g==None:
        DishoutLatex(F,s)
    else:
        n=Deg2(F)
        for mu in sorted(list(Partitions(n)),reverse=true):
            val=Scalar2(F,g(mu))
            if val<>0:
                print(latex((val,s(mu))))
                
def Printout(F,g=None):
    if F==0: 
        show()
    elif g==None:
        Printout(F,s)
    else:
        n=Deg2(F)
        for mu in sorted(list(Partitions(n)),reverse=true):
            val=Scalar2(F,g(mu))
            if val<>0:
                print((latex(val),latex(g(mu))))
                
def Printout_buthooks(F,g=None):
    if F==0: 
        show()
    elif g==None:
        Printout(F,s)
    else:
        n=Deg2(F)
        for mu in sorted(list(buthooks(n))):
            val=Scalar2(F,g(mu))
            if val<>0:
                print((latex(val),latex(s(mu))))
                
def Dishout_all(F,N):
    if F==Integer(0): 
        show()
    else:
        for n in range(N+1):
            for mu in sorted(list(Partitions(n))):
                val=Scalar2(F,s(mu))
                show((mu,val))
                
        
def Dishout_hooks(F):
    if F==0: 
        show()
    else:
        n=Deg2(F)
        for mu in hooks(n):
            val=Scalar2(F,s(mu))
            if val<>0:
                show((mu,val))
        
def Dishout_but_hooks(F):
    if F==0: 
        show()
    else:
        n=Deg2(F)
        for mu in buthooks(n):
            val=Scalar2(F,s(mu))
            if val<>0:
                show((mu,val))