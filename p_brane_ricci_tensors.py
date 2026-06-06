#!/usr/bin/env cadabra2
{\mu,\nu,\rho,\sigma,\lambda,\alpha,\beta}::Indices(position=independent).
{m,n,r,s,t,u,v,q,i,j,k,l}::Indices(position=independent).

{\partial{#}}::PartialDerivative.
\delta^{a}_{b}::KroneckerDelta.

\hat{\eta}_{\mu\nu}::Symmetric.
\hat{\eta}^{\mu\nu}::Symmetric.
\tilde{g}_{mn}::Symmetric.
\tilde{g}^{mn}::Symmetric.
\tilde{\Gamma}^{m}_{n r}::TableauSymmetry(shape={2}, indices={1,2}).

# ------------------------------------------------------------
# Nonzero Christoffel rules
# ------------------------------------------------------------

GWrule := GW^{\rho}_{\mu n} -> a H**(-1) \partial_{n}{H} \delta^{\rho}_{\mu};

GIrule := GI^{i}_{\alpha\beta} ->
 - a H**(2 a - 2 b - 1) \tilde{g}^{i j} \partial_{j}{H} \hat{\eta}_{\alpha\beta};

TWrule := TW_{j} -> a (p+1) H**(-1) \partial_{j}{H};

TIrule := TI_{j} ->
 \tilde{\Gamma}^{i}_{i j}
 + b (D-p-1) H**(-1) \partial_{j}{H};

Ccoef := a (p+1) + b (D-p-1);

# ------------------------------------------------------------
# Useful derivative rules
# ------------------------------------------------------------

dHm1_r := \partial_{r}{H**(-1)} -> - H**(-2) \partial_{r}{H};
dHm1_n := \partial_{n}{H**(-1)} -> - H**(-2) \partial_{n}{H};

dHGI := \partial_{r}{H**(2 a - 2 b - 1)} ->
        (2 a - 2 b - 1) H**(2 a - 2 b - 2) \partial_{r}{H};

# ------------------------------------------------------------
# R_{mu nu}
# ------------------------------------------------------------

A1 := \partial_{r}{GI^{r}_{\mu\nu}};
A2 := TW_{r} GI^{r}_{\mu\nu};
A3 := TI_{r} GI^{r}_{\mu\nu};
A4 := - GW^{\rho}_{\nu r} GI^{r}_{\mu \rho};
A5 := - GI^{r}_{\nu \rho} GW^{\rho}_{\mu r};

substitute(A1, GIrule);
substitute(A2, TWrule);
substitute(A2, GIrule);
substitute(A3, TIrule);
substitute(A3, GIrule);
substitute(A4, GWrule);
substitute(A4, GIrule);
substitute(A5, GWrule);
substitute(A5, GIrule);

product_rule(A1);
product_rule(A2);
product_rule(A3);
product_rule(A4);
product_rule(A5);

distribute(A1);
distribute(A2);
distribute(A3);
distribute(A4);
distribute(A5);

substitute(A1, dHGI);
substitute(A1, $\partial_{r}{\hat{\eta}_{\mu\nu}} -> 0$);
substitute(A1, $\partial_{r}{a} -> 0$);
substitute(A1, $\partial_{r}{b} -> 0$);

eliminate_kronecker(A1);
eliminate_kronecker(A2);
eliminate_kronecker(A3);
eliminate_kronecker(A4);
eliminate_kronecker(A5);

canonicalise(A1);
canonicalise(A2);
canonicalise(A3);
canonicalise(A4);
canonicalise(A5);

collect_terms(A1);
collect_terms(A2);
collect_terms(A3);
collect_terms(A4);
collect_terms(A5);

Rmunu := @(A1) + @(A2) + @(A3) + @(A4) + @(A5);

canonicalise(Rmunu);
collect_terms(Rmunu);

Rmunu;

# ------------------------------------------------------------
# R_{mn}
# ------------------------------------------------------------

B1a := \partial_{q}{\tilde{\Gamma}^{q}_{m n}};
B1b := b \partial_{q}{H**(-1) \delta^{q}_{m} \partial_{n}{H}};
B1c := b \partial_{q}{H**(-1) \delta^{q}_{n} \partial_{m}{H}};
B1d := - b \partial_{q}{H**(-1) \tilde{g}_{m n} \tilde{g}^{q s} \partial_{s}{H}};

B2a := - \partial_{n}{\tilde{\Gamma}^{q}_{q m}};
B2b := - Ccoef \partial_{n}{H**(-1) \partial_{m}{H}};

B3a := \tilde{\Gamma}^{q}_{q r} \tilde{\Gamma}^{r}_{m n};
B3b := b H**(-1) \tilde{\Gamma}^{q}_{q r} \delta^{r}_{m} \partial_{n}{H};
B3c := b H**(-1) \tilde{\Gamma}^{q}_{q r} \delta^{r}_{n} \partial_{m}{H};
B3d := - b H**(-1) \tilde{\Gamma}^{q}_{q r} \tilde{g}_{m n} \tilde{g}^{r s} \partial_{s}{H};
B3e := Ccoef H**(-1) \partial_{r}{H} \tilde{\Gamma}^{r}_{m n};
B3f := b Ccoef H**(-2) \partial_{r}{H} \delta^{r}_{m} \partial_{n}{H};
B3g := b Ccoef H**(-2) \partial_{r}{H} \delta^{r}_{n} \partial_{m}{H};
B3h := - b Ccoef H**(-2) \partial_{r}{H} \tilde{g}_{m n} \tilde{g}^{r s} \partial_{s}{H};

B4 := - (p+1) a**2 H**(-2) \partial_{m}{H} \partial_{n}{H};

B5a := - \tilde{\Gamma}^{r}_{n s} \tilde{\Gamma}^{s}_{m r};

B5b := - b H**(-1) \tilde{\Gamma}^{r}_{n s} \delta^{s}_{m} \partial_{r}{H};
B5c := - b H**(-1) \tilde{\Gamma}^{r}_{n s} \delta^{s}_{r} \partial_{m}{H};
B5d := + b H**(-1) \tilde{\Gamma}^{r}_{n s} \tilde{g}_{m r} \tilde{g}^{s j} \partial_{j}{H};

B5e := - b H**(-1) \delta^{r}_{n} \partial_{s}{H} \tilde{\Gamma}^{s}_{m r};
B5f := - b H**(-1) \delta^{r}_{s} \partial_{n}{H} \tilde{\Gamma}^{s}_{m r};
B5g := + b H**(-1) \tilde{g}_{n s} \tilde{g}^{r i} \partial_{i}{H} \tilde{\Gamma}^{s}_{m r};

B5h := - b**2 H**(-2) \delta^{r}_{n} \partial_{s}{H} \delta^{s}_{m} \partial_{r}{H};
B5i := - b**2 H**(-2) \delta^{r}_{n} \partial_{s}{H} \delta^{s}_{r} \partial_{m}{H};
B5j := + b**2 H**(-2) \delta^{r}_{n} \partial_{s}{H} \tilde{g}_{m r} \tilde{g}^{s j} \partial_{j}{H};

B5k := - b**2 H**(-2) \delta^{r}_{s} \partial_{n}{H} \delta^{s}_{m} \partial_{r}{H};
B5l := - b**2 H**(-2) \delta^{r}_{s} \partial_{n}{H} \delta^{s}_{r} \partial_{m}{H};
B5m := + b**2 H**(-2) \delta^{r}_{s} \partial_{n}{H} \tilde{g}_{m r} \tilde{g}^{s j} \partial_{j}{H};

B5n := + b**2 H**(-2) \tilde{g}_{n s} \tilde{g}^{r i} \partial_{i}{H} \delta^{s}_{m} \partial_{r}{H};
B5o := + b**2 H**(-2) \tilde{g}_{n s} \tilde{g}^{r i} \partial_{i}{H} \delta^{s}_{r} \partial_{m}{H};
B5p := - b**2 H**(-2) \tilde{g}_{n s} \tilde{g}^{r i} \partial_{i}{H} \tilde{g}_{m r} \tilde{g}^{s j} \partial_{j}{H};

product_rule(B1b);
product_rule(B1c);
product_rule(B1d);
product_rule(B2b);

distribute(B1a);
distribute(B1b);
distribute(B1c);
distribute(B1d);
distribute(B2a);
distribute(B2b);

substitute(B1b, $\partial_{q}{\delta^{q}_{m}} -> 0$);
substitute(B1c, $\partial_{q}{\delta^{q}_{n}} -> 0$);

substitute(B1b, dHm1_r);
substitute(B1c, dHm1_r);
substitute(B1d, dHm1_r);
substitute(B2b, dHm1_n);

eliminate_kronecker(B1a);
eliminate_kronecker(B1b);
eliminate_kronecker(B1c);
eliminate_kronecker(B1d);
eliminate_kronecker(B2a);
eliminate_kronecker(B2b);

eliminate_kronecker(B3a);
eliminate_kronecker(B3b);
eliminate_kronecker(B3c);
eliminate_kronecker(B3d);
eliminate_kronecker(B3e);
eliminate_kronecker(B3f);
eliminate_kronecker(B3g);
eliminate_kronecker(B3h);

eliminate_kronecker(B5a);
eliminate_kronecker(B5b);
eliminate_kronecker(B5c);
eliminate_kronecker(B5d);
eliminate_kronecker(B5e);
eliminate_kronecker(B5f);
eliminate_kronecker(B5g);
eliminate_kronecker(B5h);
eliminate_kronecker(B5i);
eliminate_kronecker(B5j);
eliminate_kronecker(B5k);
eliminate_kronecker(B5l);
eliminate_kronecker(B5m);
eliminate_kronecker(B5n);
eliminate_kronecker(B5o);
eliminate_kronecker(B5p);

substitute(B5l, $\delta^{r}_{r} -> D-p-1$);

canonicalise(B1a);
canonicalise(B1b);
canonicalise(B1c);
canonicalise(B1d);
canonicalise(B2a);
canonicalise(B2b);

canonicalise(B3a);
canonicalise(B3b);
canonicalise(B3c);
canonicalise(B3d);
canonicalise(B3e);
canonicalise(B3f);
canonicalise(B3g);
canonicalise(B3h);

canonicalise(B4);

canonicalise(B5a);
canonicalise(B5b);
canonicalise(B5c);
canonicalise(B5d);
canonicalise(B5e);
canonicalise(B5f);
canonicalise(B5g);
canonicalise(B5h);
canonicalise(B5i);
canonicalise(B5j);
canonicalise(B5k);
canonicalise(B5l);
canonicalise(B5m);
canonicalise(B5n);
canonicalise(B5o);
canonicalise(B5p);

collect_terms(B1a);
collect_terms(B1b);
collect_terms(B1c);
collect_terms(B1d);
collect_terms(B2a);
collect_terms(B2b);

collect_terms(B3a);
collect_terms(B3b);
collect_terms(B3c);
collect_terms(B3d);
collect_terms(B3e);
collect_terms(B3f);
collect_terms(B3g);
collect_terms(B3h);

collect_terms(B4);

collect_terms(B5a);
collect_terms(B5b);
collect_terms(B5c);
collect_terms(B5d);
collect_terms(B5e);
collect_terms(B5f);
collect_terms(B5g);
collect_terms(B5h);
collect_terms(B5i);
collect_terms(B5j);
collect_terms(B5k);
collect_terms(B5l);
collect_terms(B5m);
collect_terms(B5n);
collect_terms(B5o);
collect_terms(B5p);

Rmn :=
  @(B1a) + @(B1b) + @(B1c) + @(B1d)
+ @(B2a) + @(B2b)
+ @(B3a) + @(B3b) + @(B3c) + @(B3d) + @(B3e) + @(B3f) + @(B3g) + @(B3h)
+ @(B4)
+ @(B5a) + @(B5b) + @(B5c) + @(B5d) + @(B5e) + @(B5f) + @(B5g)
+ @(B5h) + @(B5i) + @(B5j) + @(B5k) + @(B5l) + @(B5m) + @(B5n) + @(B5o) + @(B5p);

canonicalise(Rmn);
collect_terms(Rmn);

Rmn;

# ------------------------------------------------------------
# R_{mu n}
# ------------------------------------------------------------

C1 := \partial_{\rho}{GW^{\rho}_{\mu n}};

substitute(C1, GWrule);

product_rule(C1);
distribute(C1);

substitute(C1, $\partial_{\rho}{a} -> 0$);
substitute(C1, $\partial_{\rho}{H**(-1)} -> 0$);
substitute(C1, $\partial_{\rho}{\partial_{n}{H}} -> 0$);
substitute(C1, $\partial_{\rho}{\delta^{\rho}_{\mu}} -> 0$);

eliminate_kronecker(C1);
canonicalise(C1);
collect_terms(C1);

Rmun := @(C1);

canonicalise(Rmun);
collect_terms(Rmun);

Rmun;
