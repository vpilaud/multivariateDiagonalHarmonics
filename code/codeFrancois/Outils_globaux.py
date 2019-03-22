def Asympt_CalE_n(n):
    return (Tenseur(Un,Un)+Tenseur(Degree_bound(n,Omega(n+1)),s[1])
              +Tenseur(Degree_bound(n,s((Omega(n+1)*Omega_even(n+1)+Omega_even(n+1)))),s[2])
              +Tenseur(Degree_bound(n+1,s(Omega(n+1)*Omega_odd(n+1)-Omega_even(n+1))),s[1,1]))

def a_mn(m,n,k):
    if n>m:
        return a_mn(n,m,k)
    else:
        return (Scalar2(CalE_mn(m,n),s(hook(n,n-1-k)))
            -Skew(add(addcol(j,a_mn(m,n,j)) 
                      for j in range(k+1,n)),e[k])).restrict_partition_lengths(k,exact=false)

@cached_function
def CalA_mn(m,n=None):
    if n==None:
        return CalA_mn(m,m)
    else: return Scalar2(CalE_mn(m,n),e[n])
    
def CalC(m,mu):
    mu=Partition(mu)
    n=mu.size()
    if m==0 and mu[0]==1:
        return Un
    elif mu[0]>m:
        return 0
    elif mu==[m]:
        return Un
    elif m==2 and mu[0]==1:
        return s([n//2])
    elif m>1 and mu==Partition([m-1,1]):
        return add(s[i] for i in range(1,n))
    elif mu==Partition([1]):
        return Un
    elif mu==Partition([1,1]):
        return s([m//2])
    elif m==n+1:
        return CalC(n,mu)
    elif mu[0]==m:
        return CalC(m,rest(mu))
    elif mu[0]>1 and is_hook(mu):
        return Skew_by_F(CalC(m,colone(n)),e[mu[0]-1])
    elif m<n and mu[0]==2 and mu.length()>=n-2:
        return Skew_by_F(CalC(m,colone(n)),h[n-mu.length()])
    elif m>3 and mu==Partition([1,1,1]):
        return addcol(2,CalC(m-3,mu))+s([3*(m//3)+((m % 3)//2)])
    elif m>=4 and (mu==Partition([1,1,1,1]) or mu==Partition([2,2])):
        return addcol(3,CalC(m-4,mu))+Rec_C_4(m,mu)
    elif rotate(mu,m).size()<n:
        return CalC(m,rotate(mu,m)) 
    elif m<=3 or n<=3:
        return Scalar2(e_mn(m,n),s[mu])
    elif n<7 and m<=n:
        return Scalar2(CalE_mn(m,n),s[mu])
    elif m<n and is_hook(mu):
        return CalC(n,hook(m,m-mu[0]))
    else:
        return q*Un    
    
def CalE(m,n):
    return add(Tenseur(CalC(m,mu),s(mu)) for mu in Partitions(n))

@func_persist
def CalE_Delta(m,n):
    if m>=n:
        rep=e_mn(m,n)
        k=m//n
        for j in range(k*(n-1)):
            for rho in Parkj(k,k*(n-1)-j,n):
                rep=(rep
                     +add(Tenseur(addmu(rho,
                                        Scalar2(restrict_length1(DeltaPrime(s(flip(rho,k,n)),Eval1(e_mn(m-k*n,n),q+t))
                                                                 -Skew1(s(rho),rep),2)
                                                ,s(mu))),s(mu)) 
                          for mu in Partitions(n)))
        return rep
    else: 
        return CalE(m,n)

@cached_function
def CalE_mn(m,n=None):
    if n==None: 
        return CalE_mn(m,m)
    elif m<0: 
        return 0
    elif m==0 and n==0: 
        return Tenseur(Un,Un)
    elif m<=1: 
        return Tenseur(Un,s(e[n]))
    elif m==2:
        return add(Omega2(Tenseur(s[k],s([n//2+k+n%2,n//2-k]))) for k in range(n//2+1))
    elif min(m,n)<=3: 
        return e_mn(m,n)
    elif m==n and n<=7: 
        return Phi[n]
    elif m==n+1: 
        return CalE_mn(n,n)
    elif m==n-1 and n<7: 
        return CalH(n)
    elif m==n-1 and n==7:
        return (Make_all(e_mn(m,n))
                -Tenseur(s[1,1,1,1]+s[5,1,1]+2*s[4,1,1]+2*s[3,1,1]+s[2,1,1]+s[1,1,1],s[2,2,2,1])
                -Tenseur(s[4,1,1]+2*s[3,1,1]+2*s[2,1,1]+2*s[1,1,1],s[3,2,2]))
    elif n==4: 
        return UpArrow1(3,CalE_mn(m-n,n))+e_mn(m,n)
    elif m==7 and n==5: 
        return CalE_Delta(7,5)
    elif m==9 and n==5:
        return (CalE_Delta(9,5)
                +Addmu([2,2,2],DeltaPrime(e[1],e[5])
                       -Tenseur(s[1],e[5]))+Tenseur(s[2,2,2,1],e[5]))
    elif m==4 and n==6:
        return (add(Tenseur(Premier(mu,m),s(mu)) for mu in hooks(n))+
                add(Tenseur(Scalar2(E_mn(m,n),s(nu)),s(nu)) 
                    for nu in Set(Partitions(n)).difference(hooks(n)))
                +Tenseur(s[1,1,1],s[2,2,1,1]))
    elif m>=n:
        return CalE_Delta(m,n)
    else: 
        return Make_all(e_mn(m,n))
    
@func_persist
def CalF_mu(mu):
    return tensor([s,e])(Pleth1(N_Schur(mu),X+1))
    
@func_persist
def CalF_mn(m,n=None):
    if n==None: 
        return CalF_mn(m,m)
    else: return tensor([s,e])(Pleth1(CalE_mn(m,n),X+Un))

def CalF_mn_k(m,n):
    return ((e(Eval1(CalF_mn(m,n),q,{q})).map_coefficients(tobinom))).map_coefficients(ToBin)

def En_zeta(F):
    return ((e(Eval1(F,q,{q})).map_coefficients(tobinom))).map_coefficients(ToBin)
    
@cached_function
def CalR_n(n,N):
    return restrict_degree1(s([n])(Tenseur(add(h[k] for k in range(N)),X)),N)

def CalE_star(n):
    return q*Strip2(CalE_mn(1,1))+add(q**k*(Strip2(CalE_mn(k,k))-Strip2(CalE_mn(k-1,k-1))) for k in range(2,n+1))

def CalE_mn_evalk_plus_un(m,n):
    return (((e(Eval1(CalE_mn(m,n),1+q,{q}))).map_coefficients(tobinom))).map_coefficients(ToBin)

def CalE_mn_evalk(m,n):
    return ((Eval1(CalE_mn(m,n),q,{q}))).map_coefficients(tobinom).map_coefficients(ToBin)

@cached_function
def CalGen_mn(m,n,F_mn):
    if gcd(m,n)==1: return CalE_mn(m,n)
    elif m<0: return 0
    elif m==0 and n==0: return Tenseur(Un,Un)
    elif min(m,n)<=3: return F_mn(m,n)
    elif n==4 and m>n: return UpArrow1(3,CalF_mn(m-n,n,F_mn))+F_mn(m,n)
    else: 
        return Make_all(F_mn(m,n))


def CalH(n):
    if n<=4: 
        return Nabla(Schur_hat([n]))
    elif n==5:
        return Nabla(Schur_hat([5]))+Tenseur(s[1,1,1],e[5])
    elif n==6:
        return (Nabla(Schur_hat([6]))
                +Tenseur(s[1,1,1,1],s[1,1,1,1,1,1])
                +Tenseur(s[2,1,1],s[2,2,2])
                +Tenseur(s[3,1,1]+s[4,1,1]+s[5,1,1],e[6])
                +Tenseur(s[1,1,1]+s[2,1,1]+s[3,1,1],s[2,2,1,1])
                +Tenseur(s[1,1,1]+s[2,1,1]+s[3,1,1]+s[4,1,1],s[2,1,1,1,1]))
    else: return 0


def CalH_mn(m,n=None):
    if n==None: return CalH_mn(m,m)
    elif gcd(m,n)==1: return CalE_mn(m,n)
    elif m<0: return 0
    elif m==0 and n==0: return Tenseur(Un,Un)
    elif min(m,n)<=3: return H_mn(m,n)
    elif m==n: return CalH(n)
    elif n==4: return UpArrow1(3,CalH_mn(m-n,n))+H_mn(m,n)
    else: 
        return Make_all(H_mn(m,n))
    
def CalM_mu(mu,a=None,b=None):
    if a==None or b==None: return CalM_mu(mu,1,1)
    else: 
        mu=Partition(mu)
        d=mu.size()
        return Make_all(M_mu(mu,a,b))   
    
def CalP_mn(m,n=None):
    if n==None: return CalP_mn(m,m)
    else: 
        return CalF_mn(m,n,P_mn)
    
def CalQ_mn(m,n=None):
    if n==None: return CalQ_mn(m,m)
    elif gcd(m,n)==1: return CalE_mn(m,n)
    elif m<0: return 0
    elif m==0 and n==0: return Tenseur(Un,Un)
    elif min(m,n)<=3: return Q_mn(m,n)
    elif m==n: return Skew1(Un+e[1],CalE_mn(m,n))
    elif n==4: return UpArrow1(3,CalQ_mn(m-n,n))+Q_mn(m,n)
    else: 
        return Make_all(Q_mn(m,n))

@cached_function
def CalQ_n(n):
    return add(N_Schur(hook(n,k)) for k in range(n))

@cached_function
def Decal_CalE_mn(m,n):
    return tensor([s,e])(Pleth1(CalE_mn(m,n),Un+X))

def Diff_Cal(m,n):
    d=Deg1(E_mn(m,n))-Deg1(E_mn(m,n-1))
    return (CalE_mn(m,n)
            -E_mn(m,n)
            -odot(Tenseur(s[d],s[1]),CalE_mn(m,n-1)-E_mn(m,n-1))
            -odot(Tenseur(e[n-1],Un),CalE_mn(m-n,n)))

@cached_function
def Direct_equerres(m,n):
    if m==1 or n-m<=1:
        correction=-x**(n-m)/(1-q)+x**(n-m)
    else:
        correction=0
    return (add(factor(c)*x**(mu.length()-1) for mu,c in (Eval1(CalE_mn(m,n),t*(1-q))/(1-q)) 
                 if mu in hooks(n))+correction).substitute({q:-q})

def epsilon(k,m,n=None):
    if n==None: return Skew1(e[k],length_component(CalE_mn(m),k))
    return Skew1(e[k],length_component(CalE_mn(m,n),k))

def Epsilon(mu,k):
    return Skew1(e[k],length_component(N_Schur(mu),k))

def E_Hooks(m,n):
    if m<n:
        return add_rows(colone(n-m),E_Hooks(n,m))
    else:
        return add(Strip(restrict_to_hooks(CalC(m,hook(n,k))))*s(hook(n,k)) for k in range(n))
    

def fff(n,k,F):
    return s(1/k*formule_F_f(n,F)(X*k))

def Formal(G):
    n=G.degree()
    return add(Tenseur(Formal_coeff(G,mu),s(mu)) for mu in Partitions(n))

def Formal_P(mu,k=None):
    if k==None:
        return Formal_P(mu,0)
    else:
        return add(Tenseur(Formal_q((q**k*c).numerator()),s(mu)) for (mu,c) in s(P(mu)))
    
def Formule_CalE_nn_1(n):
    return add(add(Tenseur(s[tau.cocharge()],s(mu)) for tau in StandardTableaux(mu)) for mu in Partitions(n))

def Formula_CalE_mn_hooks(m,n):
    if m==0:
        return Tenseur(Un,e[n])
    return add(Tenseur(Ombre(Hook_shape_formula(m,n,k)),s_hook(n,k)) for k in range(n))

def formule_carre_trois(n):
    return ppp(n,2*n+1,central_moins)

def Formule_equerres(n):
    return add(factor(mul(t^i+q for i in range(1,k))*q_binomial(n-1,k,t))*(x*t)^k for k in range(n))

def formule_F_f(n,F):
    return add(mul(F(k) for k in mu)*f(mu) for mu in Partitions(n))

def formule_F_p(n,F):
    return (-1)**n*add(mul(F(k) for k in mu)*p(mu)(-X)/zee(mu) for mu in Partitions(n))

def formule_hook(k,n):
    return q**((k-1)*binomial(n,2)+n-1)*mul(q**i+t for i in range(1,n-1))

def formule_hook_component(n,k=None):
    if k==None: return formule_hook_component(n,n-1)
    else: return Ombre(q**k*gaussian_binomial(n-1,k)*mul(q**j+t for j in range(1,k)))

def formule_nabh_trois(n):
    return 1/(2*n-1)*ppp(n,2*n-1,central_moins)

def formule_nabh_three(n):
    return 1/(2*n-1)*fff(n,2*n-1,catalan_number)

def formule_nabe_three(n):
    return 1/(n+1)*fff(n,n+1,catalan_plus)

def formule_nabe_trois(n):
    return 1/(n+1)*ppp(n,n+1,central)

def formule_squarre_three(n):
    return fff(n,2*n+1,catalan_number)

def HH_mn(n):
    return (odot(Tenseur(s[0],s[1]),CalE_mn(n-1))
            +Addmu2([1,1],CalE_mn(n-2))
            -Tenseur(Un,s[n-1,1]))

def hook_component(n,mu):
    return restrict_to_hooks(Scalar2(CalE_mn(n),s(mu)))

def Hook_shape_formula(m,n,k):
    if k<0 or k>n-1:
        return 0
    elif m==1 and k==n-1:
        return 1
    elif m<n: 
        return Hook_shape_formula(n,m,m-n+k)
    elif m==n or m==n+1: 
        return q**(binomial(n,2)-binomial(n-k,2))*gaussian_binomial(n-1,k).substitute({q:1/q})*mul(1+t/q**j 
                                                                                                   for j in range(1,k))
    elif m<2*n: 
        return Strip(Scalar2(RestrictHooks(E_mn(m,n)),s_hook(n,k)))
    else:
        return q**(degmn(m,n)-binomial(n-k,2))*gaussian_binomial(n-1,k).substitute({q:1/q})*mul(1+t/q**j 
                                                                                                for j in range(1,n-1))

@cached_function
def Make_alt(F):
    n=Deg2(F)
    res=Scalar2(F,e[n])
    for k in range(n): 
        diff=addmu(colone(n-1-k),Scalar2(F,s(hook(n,k)))-Skew(res,e[n-1-k]).restrict_partition_lengths(2,exact=false))
        res=res+diff
    return res

@cached_function
def Make_hooks(F):
    n=Deg2(F)
    F_e=Make_alt(F)
    return add(Tenseur(Skew(F_e,e[n-1-k]),s(hook(n,k))) for k in range(n))

@cached_function
def Make_all(F):
    n=Deg2(F)
    F_e=Make_alt(F)
    return F+add(Tenseur(restreint_plus_trois(Skew(F_e,s(conj_skew(mu)))),s(mu)) for mu in Partitions(n))

def min_comp_hook(n,k):
    if k==n-1:
        return Tenseur(Un,s(e[n]))
    elif k>n//2:
        return add(Tenseur(s[k].skew_by(s[i]),s(Partition([n-i,i]).conjugate())) for i in range(n//2+1))
    else:
        return add(Tenseur(s[k].skew_by(s[i]),s(Partition([n-i,i]).conjugate())) for i in range(k+1))

def N_Schur(mu):
    a=0
    b=0
    n=Partition(mu).size()
    if n<=3: 
        return Nabla(Schur_hat(mu))
    elif mu[0]==1:
        return CalE_mn(n,n)
    elif mu[0]==n:
        return CalE_mn(n-1,n)
    elif n==4 and mu==[2,2]: 
        return Nabla(Schur_hat(mu)/(q*t))
    elif n==4: 
        return Nabla(Schur_hat(mu))
    elif n==5 and mu==[2,1,1,1]:
        return (Nabla(Schur_hat(mu))
                +Tenseur(s[4,1,1],e[5])
                +Tenseur(s[3,1,1],s[2,1,1,1])
                +Tenseur(s[2,1,1],s[2,2,1]))
    elif n==5 and mu==[2,2,1]:
        return (Nabla(Schur_hat(mu)/(q*t))
                +Tenseur(s[2,1,1],e[5])
                +Tenseur(s[1,1,1],s[2,1,1,1]))
    elif n==5 and mu==[3,2]:
        return (Nabla(Schur_hat([3,2])/(q*t))
                +Tenseur(s[1,1,1],e[5]))
    elif n==5 and mu==[3,1,1]:
        return (Nabla(Schur_hat([3,1,1]))
                +Tenseur(s[3,1,1],e[5])
                +Tenseur(s[2,1,1],s[2,1,1,1])
                +Tenseur(s[1,1,1],s[2,2,1]))
    elif n==5 and mu==[4,1]:
        return (Nabla(Schur_hat([4,1]))
                +Tenseur(s[2,1,1],e[5])
                +Tenseur(s[1,1,1],s[2,1,1,1]))
    elif n==6 and mu==[5,1]:
        return (Nabla(Schur_hat([5,1]))
                +(Tenseur(s[2,1,1,1]+s[3,2,1]+s[4,1,1]
                          +s[4,2,1]+s[5,1,1]+s[6,1,1],e[6])
                  +Tenseur(s[1,1,1]+s[2,1,1]+s[3,1,1],s[3,1,1,1])
                  +Tenseur(s[1,1,1,1]+s[2,1,1]+s[2,2,1]
                           +2*s[3,1,1]+s[3,2,1]+2*s[4,1,1]+s[5,1,1],s[2,1,1,1,1])
                  +Tenseur(2*s[2,1,1]+s[2,2,1]+2*s[3,1,1]+s[4,1,1],s[2,2,1,1])
                  +Tenseur(s[1,1,1]+s[2,1,1],s[3,2,1])
                  +Tenseur(s[2,1,1]+s[3,1,1],s[2,2,2])))
    elif n==6 and mu==[4,1,1]:
        return (Nabla(Schur_hat(mu))
                +(Tenseur(s[3,1,1,1]+s[3,3,1]+s[4,2,1]+s[5,1,1]
                          +s[5,2,1]+s[6,1,1]+s[7,1,1],e[6])
                  +Tenseur(s[2,1,1,1]+s[3,1,1]+2*s[3,2,1]+2*s[4,1,1]
                           +s[4,2,1]+2*s[5,1,1]+s[6,1,1],s[2,1,1,1,1])
                  +Tenseur(s[1,1,1,1]+s[2,2,1]+3*s[3,1,1]+s[3,2,1]
                           +2*s[4,1,1]+s[5,1,1],s[2,2,1,1])
                  +Tenseur(2*s[2,1,1]+s[3,1,1],s[3,2,1])
                  +Tenseur(a*s[1,1,1],s[3,3])
                  +Tenseur(s[2,1,1]+s[2,2,1]+s[3,1,1]+s[4,1,1],s[3,1,1,1])
                  +Tenseur(s[2,1,1]+s[3,1,1]+s[4,1,1],s[2,2,2])))
    elif n==6 and mu==[3,1,1,1]:
        return (Nabla(Schur_hat(mu))
                +(Tenseur(s[4,1,1,1]+s[4,3,1]+s[5,2,1]+s[6,1,1]
                          +s[6,2,1]+s[7,1,1]+s[8,1,1], s[1,1,1,1,1,1])
                  +Tenseur(s[3,1,1,1]+s[3,3,1]+s[4,1,1]+2*s[4,2,1]
                           +2*s[5,1,1]+s[5,2,1]+2*s[6,1,1]+s[7,1,1], s[2,1,1,1,1])
                  +Tenseur(s[2,1,1,1]+2*s[3,2,1]+3*s[4,1,1]+s[4,2,1]
                           +2*s[5,1,1]+s[6,1,1], s[2,2,1,1])
                  +Tenseur(s[2,2,1]+2*s[3,1,1]+s[4,1,1], s[3,2,1])
                  +Tenseur(b*s[2,1,1],s[3,3])
                  +Tenseur(s[3,1,1]+s[3,2,1]+s[4,1,1]+s[5,1,1],s[3,1,1,1])
                  +Tenseur(a*s[1,1,1,1]+a*s[2,2,1]+s[3,2,1]+s[3,1,1]+s[4,1,1]+s[5,1,1],s[2,2,2])))
    elif n==6 and mu==[2,1,1,1,1]:
        return (Nabla(Schur_hat(mu))
                +(Tenseur(s[5,1,1,1]+s[4,3,1]+s[5,3,1]+s[6,2,1]
                          +s[7,1,1]+s[7,2,1]+s[8,1,1]+s[9,1,1],e[6])
                  +Tenseur(s[3,2,1]+s[4,1,1]+s[4,2,1]+s[5,1,1]+s[6,1,1],s[3,1,1,1])
                  +Tenseur(s[3,3,1]+s[4,1,1,1]+s[4,2,1]+s[4,3,1]+s[5,1,1]+2*s[5,2,1]
                           +2*s[6,1,1]+s[6,2,1]+2*s[7,1,1]+s[8,1,1],s[2,1,1,1,1])
                  +Tenseur(s[3,1,1,1]+s[3,2,1]+s[3,3,1]+s[4,1,1]+2*s[4,2,1]
                           +3*s[5,1,1]+s[5,2,1]+2*s[6,1,1]+s[7,1,1],s[2,2,1,1])
                  +Tenseur(s[2,2,1]+s[3,1,1]+s[3,2,1]+2*s[4,1,1]+s[5,1,1],s[3,2,1])
                  +Tenseur(s[3,1,1],s[3,3])
                  +Tenseur(s[2,1,1,1]+s[3,2,1]+s[4,2,1]
                           +s[4,1,1]+s[5,1,1]+s[6,1,1], s[2,2,2])))
    elif n==6 and mu==[4,2]:
        return (Nabla(Schur_hat(mu)/(q*t))
                +Tenseur(s[1,1,1,1]+s[3,1,1]+s[4,1,1],e[6])
                +Tenseur(s[1,1,1]+s[2,1,1]+s[3,1,1],s[2,1,1,1,1])
                +Tenseur(s[1,1,1]+s[2,1,1],s[2,2,1,1]))
    else:
        return Nabla(Schur_hat(mu))
    
def Ombre(F):
    if F==0:
        return 0
    elif F==1:
        return Un
    return  add(c*s([mqt.degrees()[0]-mqt.degrees()[1]]+[1 for j in range(mqt.degrees()[1])]) 
                for c,mqt in F.numerator())

def omega_star(F):
    if F==0: 
        return 0
    else:
        return (F.map_coefficients(lambda c: c.substitute(q=1/q,t=1/t))).omega()

def Omega_star(F):
    if F==0: return 0
    else:
        d=F.degree()-1
        return ((-q*t)**d*F).map_coefficients(lambda c: c.substitute(q=1/q,t=1/t)).omega()
    
def oper_p1_perp(F):
    n=F.degree()
    g=(e[1,1]-e[2](M*X)/M)
    return Formal_F(s(add((g(B_mu(mu)*Un)/B_mu(mu)*e[n](B_mu(mu)*Un)).scalar(Un)*c*H(mu)
               for mu,c in H(F))))

def Oper_p1_perp(n):
    F=(-1)^(n-1)*p[n](X/M)
    g=(M*e[1,1]-e[2](M*X))
    return Delta(g,F.nabla())

def Phi_mnk(m,n,k):
    return gcd(n,m*n+k)/(m*n+k)*varphi(m,n)((m*n+k)*X)

def Phi_q(r,n):
    return 1/q*varphi(r,n)(q*X,exclude={q})

def ppp(n,k,F):
    return s(1/k*formule_F_p(n,F)(X*k))

def Rec_C_4(m,mu):
    k=m//4
    j=m % 4
    if m<=14:
        return Scalar2(e_mn(m,4),s(mu))
    else:
        return (addmu([6],Rec_C_4(4*(k-1)+j,mu))
                +addmu([3,1],Rec_C_4(4*(k-1)+j,mu))
                +addmu([4,4],Rec_C_4(4*(k-2)+j,mu))
                -addmu([9,1],Rec_C_4(4*(k-2)+j,mu))
                -addmu([7,5],Rec_C_4(4*(k-3)+j,mu))
                -addmu([10,4],Rec_C_4(4*(k-3)+j,mu))
                +addmu([13,5],Rec_C_4(4*(k-4)+j,mu)))

@cached_function
def Restreint_mn_Un(m,n,j):
    if j>n-1:
        return 0
    elif m%n==1:
        return Restreint_mn_Un(m-1,n,j)
    else:
        return P(Eval1(DeltaPrime(e[n-1-j],Eval1(E_mn(m,n),q+t)),q))

@cached_function
def Rest_mn_Un(m,n,j):
    if j>n:
        return 0
    elif m%n==1:
        return Rest_mn_Un(m-1,n,j)
    else:
        return P(Eval1(Delta(e[j],Eval1(E_mn(m,n),q+t)),q))
    
@cached_function
def squarre_trois(n):
    return Eval1(CalQ_n(n),3)

@cached_function
def StripHooks(m,n):
    if m<n:
        return add_rows(colone(n-m),StripHooks(n,m))
    else:
        return add(Strip(Scalar2(RestrictHooks(CalE_mn(m,n)),s(mu)))*s(mu) for mu in hooks(n))

def StripHooks_formnula(m,n):
    if m<n:
        return add_rows(colone(n-m),StripHooks_formnula(n,m))
    else:
        return add(Hook_shape_formula(m,n,k)*s(hook(n,k)) for k in range(n))

@cached_function
def test_delta_f_emn(f,m,n):
    k=f.degree()
    return (q^(k*(n-1))*Eval1(delta_f_emn(f,m,n),q+1/q)/q_eval(f,n))

def trivarie_E(r,n):
    m=r*n+1
    return 1/m**2*varphi(r,n)(((r+1)*m)*X)

def trivarie_H(r,n):
    m=r*n-1
    return 1/(m+n)**2*varphi(r,n)(((r+1)*m+1)*X)

def trivarie_Q(r,n):
    m=r*n
    return 1/((r+1)*m+1)*varphi(r,n)(((r+1)*m+1)*X) 

def Undefined(m,n):
    res=Scalar1(CalE(m,n),Un)-Scalar1(CalE(m,n),Un).map_coefficients(lambda c: c.substitute({q:0}))
    return Set([mu for mu,c in res])

def varphi(r,n):
    return add(C_r_mu(r,mu)*f(mu) for mu in Partitions(n))

def VdM(A):
    n=A.__len__()
