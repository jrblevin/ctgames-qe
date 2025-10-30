# example.ctl -- example quality ladder control file  -*-conf-*-

# Model
nn          3
nw          4
we          2
mktsize     0.15
rho         0.09
maxiter     10000
vf_tol      1.0e-13
Q_tol       1.0e-13
delta       1.00
discrete    T
nvf         100
prmod       0

# Model Parameters:
#            lam_l   lam_h   gam   kap   eta   fc
theta        0.8     1.2     0.4   0.6   4.0   0.8

# Monte Carlo parameters
nmc         25
nm          400
nt          200

# L-BFGS-B Optimization
eps_lbfgs     1.0e-12
pgtol_lbfgs   1.0e-9
iprint_lbfgs  1
