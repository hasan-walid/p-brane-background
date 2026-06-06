#!/usr/bin/env cadabra2
{\mu,\nu,\rho,\sigma,\lambda}::Indices(position=independent).
{m,n,r,s,t,u}::Indices(position=independent).

{\partial{#}}::PartialDerivative.
\delta^{a}_{b}::KroneckerDelta.

\hat{\eta}_{\mu\nu}::Symmetric.
\hat{\eta}^{\mu\nu}::Symmetric.
\tilde{g}_{mn}::Symmetric.
\tilde{g}^{mn}::Symmetric.

# ----------------------------------------------------------------
# 1. Gamma^rho_{mu nu}
# Do NOT write: \Gamma^{\rho}_{\mu\nu} = ...
# Just define the RHS first.
# ----------------------------------------------------------------

GammaWW := 1/2 H**{-2 a}\hat{\eta}^{\rho\sigma}(
    \partial_{\mu}{(H**{2 a}\hat{\eta}_{\sigma\nu})}
  + \partial_{\nu}{(H**{2 a}\hat{\eta}_{\sigma\mu})}
  - \partial_{\sigma}{(H**{2 a}\hat{\eta}_{\mu\nu})}
);

product_rule(GammaWW);
distribute(GammaWW);

substitute(GammaWW, $\partial_{\mu}{H**{2 a}} -> 0$);
substitute(GammaWW, $\partial_{\nu}{H**{2 a}} -> 0$);
substitute(GammaWW, $\partial_{\sigma}{H**{2 a}} -> 0$);

substitute(GammaWW, $\partial_{\mu}{\hat{\eta}_{\sigma\nu}} -> 0$);
substitute(GammaWW, $\partial_{\nu}{\hat{\eta}_{\sigma\mu}} -> 0$);
substitute(GammaWW, $\partial_{\sigma}{\hat{\eta}_{\mu\nu}} -> 0$);

canonicalise(GammaWW);
collect_terms(GammaWW);

GammaWW;
# ----------------------------------------------------------------
# 2. Gamma^rho_{mu n}
# ----------------------------------------------------------------

GammaWI := 1/2 H**{-2 a}\hat{\eta}^{\rho\sigma}
           \partial_{n}{(H**{2 a}\hat{\eta}_{\sigma\mu})};

product_rule(GammaWI);
distribute(GammaWI);

substitute(GammaWI, $\partial_{n}{H**{2 a}} -> 2 a H**{2 a-1}\partial_{n}{H}$);
substitute(GammaWI, $\partial_{n}{\hat{\eta}_{\sigma\mu}} -> 0$);

substitute(GammaWI, $\hat{\eta}^{\rho\sigma}\hat{\eta}_{\sigma\mu} -> \delta^{\rho}_{\mu}$);

canonicalise(GammaWI);
eliminate_kronecker(GammaWI);
collect_terms(GammaWI);

GammaWI;
# ----------------------------------------------------------------
# 3. Gamma^m_{mu nu}
# ----------------------------------------------------------------

GammaIWW := -1/2 H**{-2 b}\tilde{g}^{m r}
            \partial_{r}{(H**{2 a}\hat{\eta}_{\mu\nu})};

product_rule(GammaIWW);
distribute(GammaIWW);

substitute(GammaIWW, $\partial_{r}{H**{2 a}} -> 2 a H**{2 a-1}\partial_{r}{H}$);
substitute(GammaIWW, $\partial_{r}{\hat{\eta}_{\mu\nu}} -> 0$);

canonicalise(GammaIWW);
collect_terms(GammaIWW);

GammaIWW;
# ----------------------------------------------------------------
# 4. Gamma^m_{nr}
# ----------------------------------------------------------------

GammaIII := 1/2 H**{-2 b}\tilde{g}^{m s}(
    \partial_{n}{(H**{2 b}\tilde{g}_{s r})}
  + \partial_{r}{(H**{2 b}\tilde{g}_{s n})}
  - \partial_{s}{(H**{2 b}\tilde{g}_{n r})}
);

product_rule(GammaIII);
distribute(GammaIII);

substitute(GammaIII, $\partial_{n}{H**{2 b}} -> 2 b H**{2 b-1}\partial_{n}{H}$);
substitute(GammaIII, $\partial_{r}{H**{2 b}} -> 2 b H**{2 b-1}\partial_{r}{H}$);
substitute(GammaIII, $\partial_{s}{H**{2 b}} -> 2 b H**{2 b-1}\partial_{s}{H}$);

canonicalise(GammaIII);
collect_terms(GammaIII);

GammaIII;


