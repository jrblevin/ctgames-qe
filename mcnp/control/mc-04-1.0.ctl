# Model
nn           4
nw           7
we           4
mktsize      0.60
rho          0.05
maxiter      20000
vf_tol       1.0e-13
Q_tol        1.0e-13
delta        1.00
discrete     T
nvf          100

# L-BFGS-B Optimization
eps_lbfgs    1.0e-8
pgtol_lbfgs  1.0e-12

# Model Parameters:
#            lam_l   lam_h   gam   kap   eta   fc
theta        1.0     1.2     0.4   0.8   4.0   0.9

# Monte Carlo parameters
nmc          100
nm           10000
restart      mc-04-1.0-restart.ctl
