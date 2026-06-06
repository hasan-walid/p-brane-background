#!/usr/bin/env cadabra2
# ----------------------------------------------------------
# Index sets
# ----------------------------------------------------------
{M,N,P,Q,R,S,T,U#}::Indices(full, position=fixed).
{\mu,\nu,\rho,\sigma,\lambda,\kappa#}::Indices(ext, parent=full, position=fixed).
{m,n,p,q,r,s#}::Indices(int, parent=full, position=fixed).

# ----------------------------------------------------------
# Basic tensors
# ----------------------------------------------------------
g_{M N}::Symmetric.
g^{M N}::Symmetric.

\eta_{\mu\nu}::Symmetric.
\eta^{\mu\nu}::Symmetric.

\hat{g}_{m n}::Symmetric.
\hat{g}^{m n}::Symmetric.

\delta^{M}_{N}::KroneckerDelta.

\partial_{#}::PartialDerivative.

# ----------------------------------------------------------
# Metric ansatz and inverse metric
# ----------------------------------------------------------
metric_rules := {
  g_{\mu\nu} -> \Omega^{2}\eta_{\mu\nu},
  g_{m n}    -> \Omega^{2}\hat{g}_{m n},
  g_{\mu m}  -> 0,
  g_{m \mu}  -> 0,

  g^{\mu\nu} -> \Omega^{-2}\eta^{\mu\nu},
  g^{m n}    -> \Omega^{-2}\hat{g}^{m n},
  g^{\mu m}  -> 0,
  g^{m \mu}  -> 0
};

# ----------------------------------------------------------
# Dependence rules:
# Omega and \hat g depend only on y,
# eta depends only on x
# ----------------------------------------------------------
dep_rules := {
  \partial_{\mu}{\Omega} -> 0,

  \partial_{\mu}{\hat{g}_{m n}} -> 0,
  \partial_{\mu}{\hat{g}^{m n}} -> 0,

  \partial_{m}{\eta_{\mu\nu}} -> 0,
  \partial_{m}{\eta^{\mu\nu}} -> 0
};

# ----------------------------------------------------------
# Inverse-metric contractions
# ----------------------------------------------------------
inv_rules := {
  \eta^{\mu\rho}\eta_{\rho\nu} -> \delta^{\mu}_{\nu},
  \eta_{\mu\rho}\eta^{\rho\nu} -> \delta_{\mu}^{\nu},

  \hat{g}^{m p}\hat{g}_{p n} -> \delta^{m}_{n},
  \hat{g}_{m p}\hat{g}^{p n} -> \delta_{m}^{n}
};

# ----------------------------------------------------------
# Define the Christoffel symbol from the full metric
# ----------------------------------------------------------
GammaRule := \Gamma^{P}_{M N} ->
  (1/2) g^{P Q}(
      \partial_{M}{g_{N Q}}
    + \partial_{N}{g_{M Q}}
    - \partial_{Q}{g_{M N}}
  );

# ----------------------------------------------------------
# Compute only the nonzero Christoffels needed for R_{mu nu}
# ----------------------------------------------------------

# Gamma^rho_{mu nu}
G1 := \Gamma^{\rho}_{\mu\nu};
substitute(G1, GammaRule);
substitute(G1, metric_rules);
distribute(G1);
product_rule(G1);
substitute(G1, dep_rules);
substitute(G1, inv_rules);
canonicalise(G1);
rename_dummies(G1);
collect_terms(G1);

# Gamma^m_{mu nu}
G2 := \Gamma^{m}_{\mu\nu};
substitute(G2, GammaRule);
substitute(G2, metric_rules);
distribute(G2);
product_rule(G2);
substitute(G2, dep_rules);
substitute(G2, inv_rules);
canonicalise(G2);
rename_dummies(G2);
collect_terms(G2);

# Gamma^rho_{mu n}
G3 := \Gamma^{\rho}_{\mu n};
substitute(G3, GammaRule);
substitute(G3, metric_rules);
distribute(G3);
product_rule(G3);
substitute(G3, dep_rules);
substitute(G3, inv_rules);
canonicalise(G3);
rename_dummies(G3);
collect_terms(G3);

# Gamma^m_{n p}
G4 := \Gamma^{m}_{n p};
substitute(G4, GammaRule);
substitute(G4, metric_rules);
distribute(G4);
product_rule(G4);
substitute(G4, dep_rules);
substitute(G4, inv_rules);
canonicalise(G4);
rename_dummies(G4);
collect_terms(G4);

# ----------------------------------------------------------
# For convenience, define symbols for the external and internal
# Levi-Civita connections
# ----------------------------------------------------------
GammaEta := \Gammaeta^{\rho}_{\mu\nu} ->
  (1/2)\eta^{\rho\sigma}(
      \partial_{\mu}{\eta_{\nu\sigma}}
    + \partial_{\nu}{\eta_{\mu\sigma}}
    - \partial_{\sigma}{\eta_{\mu\nu}}
  );

GammaHat := \Gammahat^{m}_{n p} ->
  (1/2)\hat{g}^{m q}(
      \partial_{n}{\hat{g}_{p q}}
    + \partial_{p}{\hat{g}_{n q}}
    - \partial_{q}{\hat{g}_{n p}}
  );

# Replace the computed Christoffels by simpler named objects
G1rule := \Gamma^{\rho}_{\mu\nu} -> \Gammaeta^{\rho}_{\mu\nu};
G2rule := \Gamma^{m}_{\mu\nu} -> -\hat{g}^{m n}\eta_{\mu\nu}\Omega^{-1}\partial_{n}{\Omega};
G3rule := \Gamma^{\rho}_{\mu n} -> \delta^{\rho}_{\mu}\Omega^{-1}\partial_{n}{\Omega};
G4rule := \Gamma^{m}_{n p} ->
   \Gammahat^{m}_{n p}
 + \delta^{m}_{n}\Omega^{-1}\partial_{p}{\Omega}
 + \delta^{m}_{p}\Omega^{-1}\partial_{n}{\Omega}
 - \hat{g}_{n p}\hat{g}^{m q}\Omega^{-1}\partial_{q}{\Omega};

# ----------------------------------------------------------
# Ricci tensor R_{mu nu}
# Split the summed index P into external rho and internal m
# ----------------------------------------------------------
Rmunu := 
    \partial_{\rho}{\Gamma^{\rho}_{\mu\nu}}
  + \partial_{m}{\Gamma^{m}_{\mu\nu}}
  - \partial_{\nu}{\Gamma^{\rho}_{\mu\rho}}
  - \partial_{\nu}{\Gamma^{m}_{\mu m}}
  + \Gamma^{\rho}_{\rho\lambda}\Gamma^{\lambda}_{\mu\nu}
  + \Gamma^{\rho}_{\rho n}\Gamma^{n}_{\mu\nu}
  + \Gamma^{m}_{m\lambda}\Gamma^{\lambda}_{\mu\nu}
  + \Gamma^{m}_{m n}\Gamma^{n}_{\mu\nu}
  - \Gamma^{\rho}_{\nu\lambda}\Gamma^{\lambda}_{\mu\rho}
  - \Gamma^{\rho}_{\nu n}\Gamma^{n}_{\mu\rho}
  - \Gamma^{m}_{\nu\lambda}\Gamma^{\lambda}_{\mu m}
  - \Gamma^{m}_{\nu n}\Gamma^{n}_{\mu m};

substitute(Rmunu, G1rule);
substitute(Rmunu, G2rule);
substitute(Rmunu, G3rule);
substitute(Rmunu, G4rule);

distribute(Rmunu);
product_rule(Rmunu);
substitute(Rmunu, dep_rules);
substitute(Rmunu, inv_rules);
substitute(Rmunu, GammaEta);
substitute(Rmunu, GammaHat);
distribute(Rmunu);
canonicalise(Rmunu);
eliminate_kronecker(Rmunu);
rename_dummies(Rmunu);
collect_terms(Rmunu);

Rmunu;
{M,N,P,Q,R,S,T,U#}::Indices(full, position=fixed).
{\mu,\nu,\rho,\sigma,\lambda,\kappa#}::Indices(ext, parent=full, position=fixed).
{m,n,p,q,r,s#}::Indices(int, parent=full, position=fixed).

\eta_{\mu\nu}::Symmetric.
\eta^{\mu\nu}::Symmetric.
\hat{g}_{m n}::Symmetric.
\hat{g}^{m n}::Symmetric.
\delta^{M}_{N}::KroneckerDelta.

D::Symbol.
d::Symbol.
qdim::Symbol.

Delta := \partial_{m}{\Gamma^{m}_{\mu\nu}}
       + \Gamma^{\rho}_{\rho n}\Gamma^{n}_{\mu\nu}
       + \Gamma^{m}_{m n}\Gamma^{n}_{\mu\nu}
       - \Gamma^{\rho}_{\nu n}\Gamma^{n}_{\mu\rho}
       - \Gamma^{m}_{\nu\lambda}\Gamma^{\lambda}_{\mu m};

rule1 := \Gamma^{\rho}_{\mu\nu} -> \Gammaeta^{\rho}_{\mu\nu};
rule2 := \Gamma^{m}_{\mu\nu} -> -\eta_{\mu\nu} H^{m};
rule3 := \Gamma^{\rho}_{\mu n} -> \delta^{\rho}_{\mu} A_{n};
rule4 := \Gamma^{\rho}_{n\mu} -> \delta^{\rho}_{\mu} A_{n};
rule5 := \Gamma^{m}_{\mu n} -> 0;
rule6 := \Gamma^{m}_{n\mu} -> 0;
rule7 := \Gamma^{m}_{n p} -> \Gammahat^{m}_{n p}
                           + \delta^{m}_{n} A_{p}
                           + \delta^{m}_{p} A_{n}
                           - \hat{g}_{n p} H^{m};

rule8 := A_{m} -> \Omega^{-1}\partial_{m}{\Omega};
rule9 := H^{m} -> \hat{g}^{m n} A_{n};

rule10 := \delta^{\rho}_{\rho} -> d;
rule11 := \delta^{m}_{m} -> qdim;

substitute(Delta, rule1);
substitute(Delta, rule2);
substitute(Delta, rule3);
substitute(Delta, rule4);
substitute(Delta, rule5);
substitute(Delta, rule6);
substitute(Delta, rule7);
substitute(Delta, rule8);
substitute(Delta, rule9);
substitute(Delta, rule10);
substitute(Delta, rule11);

eliminate_kronecker(Delta);
canonicalise(Delta);
rename_dummies(Delta);
collect_terms(Delta);

Delta;
