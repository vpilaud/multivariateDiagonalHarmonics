Qqt=QQ['q','t']
Rqt=FractionField(Qqt)
Qqtuv=QQ['q','t','u','v','a']
Rqtuv=FractionField(Qqtuv)
Sym = SymmetricFunctions(Rqtuv)
Sym.inject_shorthands(verbose=False)
H = Sym.macdonald().Ht()
t=H.t
q=H.q
Un=s([])

u=(Rqtuv.gens()[2]) 
v=(Rqtuv.gens()[3]) 
a=(Rqtuv.gens()[4]) 




Z=s[1]
X=tensor([Z,Un])
Y=tensor([Un,Z])
M=(1-q)*(1-t)

S = SymmetricFunctions(QQ)


f = e.dual_basis(prefix='f')

def MacCombScalar(mu):
    mu=Partition(mu)
    return (-1)**(mu.size()-mu.length())*(p(mu)(M*Un).scalar(Un))*mu.aut()

def WhittScalar(mu):
    mu=Partition(mu)
    return (p(mu)((1-q)*Un).scalar(Un))*mu.aut()

H_hat = H.dual_basis(MacCombScalar,prefix='H_hat')

def MacScalar(f,g):
    return f.scalar(g,MacCombScalar)

zero=Partition([0])

def mystr(i): 
    if i<10: 
        return str(i) 
    else: 
        return ''.join([str(i),"."])

def compact(mu):
    if mu==Partition([]): return 0
    else: return (''.join(mystr(i) for i in mu))

Partition._latex_= compact
Composition._latex_= compact




f._latex_term = lambda mu: "\\boldsymbol{1}" if mu==zero else "f_{%s}"%(''.join(mystr(i) for i in mu))
s._latex_term = lambda mu: "\\boldsymbol{1}" if mu==zero else "s_{%s}"%(''.join(mystr(i) for i in mu))
p._latex_term = lambda mu: "\\boldsymbol{1}" if mu==zero else "p_{%s}"%(''.join(mystr(i) for i in mu))
h._latex_term = lambda mu: "\\boldsymbol{1}" if mu==zero else "h_{%s}"%(''.join(mystr(i) for i in mu))
e._latex_term = lambda mu: "\\boldsymbol{1}" if mu==zero else "e_{%s}"%(''.join(mystr(i) for i in mu))
m._latex_term = lambda mu: "\\boldsymbol{1}" if mu==zero else "m_{%s}"%(''.join(mystr(i) for i in mu))
H._latex_term = lambda mu: "\\boldsymbol{1}" if mu==zero else "\\widetilde{H}_{%s}"%(''.join(mystr(i) for i in mu))
H_hat._latex_term = lambda mu: "\\boldsymbol{1}" if mu==zero else "\\widehat{H}_{%s}"%(''.join(mystr(i) for i in mu))


W=Sym.macdonald(t=0).P()

W.print_options(prefix="W",latex_prefix="W")
W_dual=W.dual_basis(WhittScalar)
W._latex_term = lambda mu: "\\boldsymbol{1}" if mu==zero else "W_{%s}"%(''.join(mystr(i) for i in mu))
W_dual._latex_term = lambda mu: "\\boldsymbol{1}" if mu==zero else "\\widehat{W}_{%s}"%(''.join(mystr(i) for i in mu))

Wt = Sym.macdonald(q=t,t=0).H()
Wt.print_options(prefix="W")
Wt_dual=Wt.dual_basis(WhittScalar)
Wt._latex_term = lambda mu: "\\boldsymbol{1}" if mu==zero else "W_{%s}"%(''.join(mystr(i) for i in mu))
Wt_dual._latex_term = lambda mu: "\\boldsymbol{1}" if mu==zero else "\\widehat{W}_{%s}"%(''.join(mystr(i) for i in mu))

HL=Sym.macdonald(t=0).Ht()
HL.print_options(prefix="H",latex_prefix="H")
HL._latex_term = lambda mu: "\\boldsymbol{1}" if mu==zero else "H_{%s}"%(''.join(mystr(i) for i in mu))


s.set_print_style('length')
e.set_print_style('length')
H.set_print_style('length')

def homog(k):
    if k<0:
        return 0*Un
    else:
        return s[k]
    


def inv_Pi_op(f):
    return s(add(c*Pi_mu(mu)**(-1)*H(mu) for (mu,c) in H(f)))


def Pi_op(f):
    if f==0:
        return 0
    else:
        return s(add(c*Pi_mu(mu)*H(mu) for (mu,c) in H(f)))


def D_k(k,g):
    if k==0:
        return g-M*Delta(e[1],g)
    elif k>0:
        return s((D_k(0,p[k](Z/M)*g)-p[k](Z/M)*D_k(0,g)))
    else:
        k=-k
        return s((D_k(0,g.skew_by(p[k]))-D_k(0,g).skew_by(p[k])))
    
def Restreint(F,condition=(lambda nu:True)):
    return add(c*s(mu) for mu,c in s(F) if condition(mu))

@cached_function
def VarPlus(g,x=q):
    return add(c*s(mu)*(s(nu.conjugate())(x*Un).scalar(Un)) for (mu,nu),c in tensor([s,s])(g(X+Y)))
    

def Fermeture_Gauche(g_alt,n):
    return add(tensor([g_alt.skew_by(s(Bar(mu.conjugate()))),s(mu)]) for mu in Partitions(n))

def InSchur(F):
    return s(Sym(S.from_polynomial(Rqt(F).numerator()))).restrict_partition_lengths(2,exact=false)

#def Frob_q(D):
#    D=Diagram(D.cells())
#    return q^add(c[1] for c in D)*s(D.specht_module(QQ).frobenius_image())

def PrintLatex(F):
    n=F.degree()
    Famille=Tous_Partages(n)
    First=true
    for mu in sorted(Famille):
        val=F.scalar(f(mu))
        if First and val.__len__()==1 and not val==0:
            print(latex(val),LatexExpr("\\otimes e_{"),latex(mu),LatexExpr("}"))
            First=false
        elif First and not val==0:
            print(LatexExpr("\\big("),latex(val),LatexExpr("\\big) \\otimes e_{"),latex(mu),LatexExpr("}"))
            First=false
        elif val.__len__()==1 and not val==0:
            print(LatexExpr("+"),latex(val),LatexExpr("\\otimes e_{"),latex(mu),LatexExpr("}"))
        elif not val==0:
            print(LatexExpr("+\\big("),
                 latex(val),LatexExpr("\\big) \\otimes e_{"),latex(mu),LatexExpr("}"))

            
def Dishout_F(F):
    if F==0: 
        show(0)
    else:
        Liste=sorted(Set([nu for (mu,nu),c in tensor([s,e])(F)]))
        First=Liste[0]
        for mu in Liste:
            val=F.scalar(f(mu))
            if mu==First and val.__len__()==Integer(1) and not val==Integer(0):
                show(val,LatexExpr("\\otimes"),e(mu))
            elif mu==First and not val==Integer(0):
                show(LatexExpr("\\big("),val,LatexExpr("\\big)"),
                     LatexExpr("\\otimes"),e(mu))
            elif val.__len__()==Integer(1) and not val==Integer(0):
                show(LatexExpr("\\qquad +\\ "),val,LatexExpr("\\otimes"),e(mu))
            elif not val==Integer(0):
                show(LatexExpr("\\qquad +\\ "),LatexExpr("\\big("),
                     val,LatexExpr("\\big)")
                         ,LatexExpr("\\otimes"),e(mu))
                

