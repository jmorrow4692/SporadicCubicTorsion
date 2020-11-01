from maartens_sage_functions import (congruence_groups_between_gamma0_and_gamma1,
                                     count_points_J_H,
                                     counts,
                                     gonality_lower_bound,
                                     cuspidal_integral_structure_matrix,
					galois_action,
					galois_orbit,
					galois_orbits,
                                     generators_of_subgroups_of_unit_group,
                                     diamond_operator,
                                     positive_part,
                                     rational_cuspidal_classgroup,
                                     tate_normal_form,
				     integral_period_mapping)

from cuspidal_classgroup import (intersection,
                                 cuspidal_rational_subgroup_mod_rational_cuspidal_subgroup,
                                 upper_bound_index_cusps_in_JG_torsion,
                                 JG_torsion_upperbound,
  				 kill_torsion_coprime_to_q,
				 cuspidal_integral_structure)

_sage_const_1 = Integer(1);

# In the below computations, we actually work with the modular curve $X_{\mu}(N)$ (and its quotients) rather than $X_1(N)$. 
# Note that there is an isomorphism $X_{\mu}(N) \to X_1(N)$ that sends the cusps to the cusps; isomorphism is $(E,f: \mu_N -> E[N]) to (E/im(f),f^{\vee}: Z/NZ -> E[N]/im(f))$. 
# The reason for working with $X_{\mu}(N)$ rather than $X_1(N)$ is that when using the model $X_{\mu}(N)$ the cusp $\infty$ in $\H^*$ becomes rational, and hence many computations become easier.
# We also note that the isomorphism $X_{\mu}(N) \to X_1(N)$ does not preserve Hecke operators, and so we use the Hecke operators $T_q - \langle q \rangle - q$ to kill the prime to $q$ torsion in $J_{\mu}(N)(\Q)$ (cf.~\cite[page 87]{DiamondIm:ModularCurves}).

def do_cusps_generate_J_H_torsion(G,bound=60):
    N = G.level()
    M=ModularSymbols(G);
    S=M.cuspidal_subspace()
    Sint=cuspidal_integral_structure(M)
    kill_mat=(M.star_involution().matrix().restrict(Sint)-_sage_const_1 )
    kill=kill_mat.transpose().change_ring(ZZ).row_module()
    for p in prime_range(3 ,bound):
        if not N % p ==0 :
            kill+=kill_torsion_coprime_to_q(p,M).restrict(Sint).change_ring(ZZ).transpose().row_module()
        #if kill.matrix().is_square() and kill.matrix().determinant()==d:
        #    #print p
        #    break
    kill_mat = kill.matrix().transpose()
    clnr1 = kill_mat.determinant()
    clnr2 = rational_cuspidal_classgroup(G).cardinality()
    bad_primes = ZZ(clnr1/clnr2).prime_divisors()

    d_kill_mat = kill_mat.smith_form()[0][-1,-1]
    pm=integral_period_mapping(M)
    period_images1=[sum([M.coordinate_vector(M([c,infinity])) for c in cusps])*pm for cusps in galois_orbits(G)]
    period_images2=[M.coordinate_vector(M([c,infinity]))*pm for c in G.cusps() if c != Cusp(oo)]
    #period_matrix1,d1=Matrix(period_images1)._clear_denom()
    period_matrix2,d2=Matrix(period_images2)._clear_denom()

    for p in bad_primes:
        pn = p**max(d_kill_mat.valuation(p),d2.valuation(p))
        kill_p = kill+pn*kill.ambient_module()
        kill_mat_p=kill_p.matrix().transpose()
        dn = pn*d2.prime_to_m_part(p)
        Sint2 = ZZ**Sint.rank()
        Vp = Sint2/(pn*Sint2)
        tors_bound = Vp.hom([i.lift()*kill_mat_p for i in Vp.gens()],codomain=Vp).kernel()
        rat_csp_grp = tors_bound.submodule([dn*i for i in period_images1])
        csp_grp = Vp.submodule([dn*i for i in period_images2])
        appr_csp_grp = tors_bound.submodule([i for i in csp_grp.V().intersection(tors_bound.V()).gens()])
        if len((tors_bound/rat_csp_grp).invariants()) != len((appr_csp_grp/rat_csp_grp).invariants()):
            return False,(tors_bound/rat_csp_grp).invariants(),(appr_csp_grp/rat_csp_grp).invariants()
    return True, rational_cuspidal_classgroup(G).invariants()


#######################################################
# Proof of Proposition 4.11
#######################################################

[ [N,cuspidal_rational_subgroup_mod_rational_cuspidal_subgroup(Gamma1(N))] for N in [10..55] ]

#######################################################
# Output
#######################################################

[[10, ()],
 [11, ()],
 [12, ()],
 [13, ()],
 [14, ()],
 [15, ()],
 [16, ()],
 [17, ()],
 [18, ()],
 [19, ()],
 [20, ()],
 [21, ()],
 [22, ()],
 [23, ()],
 [24, ()],
 [25, ()],
 [26, ()],
 [27, ()],
 [28, ()],
 [29, ()],
 [30, ()],
 [31, ()],
 [32, ()],
 [33, ()],
 [34, ()],
 [35, ()],
 [36, ()],
 [37, ()],
 [38, ()],
 [39, ()],
 [40, ()],
 [41, ()],
 [42, ()],
 [43, ()],
 [44, ()],
 [45, ()],
 [46, ()],
 [47, ()],
 [48, ()],
 [49, ()],
 [50, ()],
 [51, ()], 
 [52, ()], 
 [53, ()], 
 [54, ()], 
 [55, ()]
]



#######################################################
# Cuspidal part of J1(N)(Q)
#######################################################

[ [N,rational_cuspidal_classgroup(Gamma1(N)).invariants()] for N in [13..55] ]

#######################################################
# Output
#######################################################

 [
 [[13, (19,)],
 [14, (6,)],
 [15, (4,)],
 [16, (2, 10)],
 [17, (584,)],
 [18, (21,)],
 [19, (4383,)],
 [20, (60,)],
 [21, (364,)],
 [22, (5, 775)],
 [23, (408991,)],
 [24, (2, 2, 120)],
 [25, (227555,)],
 [26, (133, 1995)],
 [27, (3, 3, 52497)],
 [28, (2, 4, 12, 936)],
 [29, (4, 4, 64427244)],
 [30, (4, 8160)],
 [31, (10, 1772833370)],
 [32, (2, 2, 2, 4, 120, 11640)],
 [33, (5, 42373650)],
 [34, (8760, 595680)],
 [35, (13, 109148520)],
 [36, (12, 252, 7812)],
 [37, (160516686697605,)],
 [38, (9, 4383, 33595695)],
 [39, (7, 31122, 3236688)],
 [40, (2, 2, 2, 8, 120, 895440)],
 [41, (107768799408099440,)],
 [42, (182, 1092, 131040)],
 [43, (2, 1563552532984879906)],
 [44, (4, 620, 3100, 6575100)]],
 [45, (3, 9, 36, 16592750496)],
 [46, (408991, 546949390174)],
 [47, (3279937688802933030787,)],
 [48, (2, 2, 2, 2, 4, 40, 40, 240, 1436640)],
 [49, (7, 52367710906884085342)],
 [50, (5, 1137775, 47721696825)],
 [51, (8, 1168, 7211322610146240)],
 [52, (4, 28, 532, 7980, 17470957140)],
 [53, (182427302879183759829891277,)],
 [54, (3, 3, 3, 9, 9, 1102437, 1529080119)],
 [55, (5, 550, 8972396739917886000)]]
]

#######################################################
# Hecke bound for J1(N) for 13 <= N <= 55
#######################################################

[ [N,upper_bound_index_cusps_in_JG_torsion(Gamma1(N),rational_cuspidal_classgroup(Gamma1(N)).cardinality())] for N in [13..55] ]

#######################################################
# Output
#######################################################

[[13, 1],
 [14, 1],
 [15, 1],
 [16, 1],
 [17, 1],
 [18, 1],
 [19, 1],
 [20, 1],
 [21, 1],
 [22, 1],
 [23, 1],
 [24, 2],
 [25, 1],
 [26, 1],
 [27, 1],
 [28, 1],
 [29, 1],
 [30, 1],
 [31, 1],
 [32, 2],
 [33, 2],
 [34, 1],
 [35, 1],
 [36, 1],
 [37, 1],
 [38, 1],
 [39, 1],
 [40, 4],
 [41, 1],
 [42, 1],
 [43, 1],
 [44, 1],
 [45, 1],
 [46, 1],
 [47, 1],
 [48, 16],
 [49, 1],
 [50, 1],
 [51, 1],
 [52, 1],
 [53, 1],
 [54, 3],
 [55, 1]
]


##################################################################################
# Computing p-torsion of (M/Cl_Q^c X1(N)) and (M′/Cl_Q^c X1(N)) in Theorem 4.13 
##################################################################################

bound = 24
[ [N,do_cusps_generate_J_H_torsion(Gamma1(N),bound)] for N in [24,32,33,40,43,48,54] ]

####################################################### 
# Output
#######################################################

[
 [24, (True, (2, 2, 120))],
 [32, (True, (2, 2, 2, 4, 120, 11640))],
 [33, (True, (5, 42373650))]]
 [40, (True, (2, 2, 2, 8, 120, 895440))],
 [43, (True, (2, 1563552532984879906))],
 [48, (True, (2, 2, 2, 2, 4, 40, 40, 240, 1436640))],
 [54, (False, (3,), ())]
]

######################################################
# Cuspidal torsion for J0(N) for N = 30,33,35,39
######################################################

[ [N,rational_cuspidal_classgroup(Gamma0(N)).invariants()] for N in [30,33,35,39] ]

####################################################### 
# Output
#######################################################

[[30, (2, 4, 24)], [33, (10, 10)], [35, (2, 24)], [39, (2, 28)]] 

#######################################################
# Hecke bounds for J0(N) for N = 30,33,35,39
#######################################################

[ [N,upper_bound_index_cusps_in_JG_torsion(Gamma0(N),rational_cuspidal_classgroup(Gamma0(N)).cardinality())] for N in [30,33,35,39] ]

####################################################### 
#Output
#######################################################

[[30, 4], [33, 2], [35, 2], [39, 2]]

#################################################################
# Additional argument from Theorem 4.13 does not improve bounds
#################################################################

[ [N,do_cusps_generate_J_H_torsion(Gamma0(N),bound)] for N in [30,33,35,39] ]

####################################################### 
#Output
#######################################################

[[30, (False, (2, 2), ())],
 [33, (False, (2,), ())],
 [35, (False, (2,), ())],
 [39, (False, (2,), ())]]


#######################################################
# Cupsidal torsion for X_H(45)
#######################################################

rational_cuspidal_classgroup(GammaH(45,[ 1, 4, 16, 19, 31, 34 ])).invariants()

#######################################################
# Output
#######################################################

(2, 4, 48) 

#######################################################
# Hecke bound for XH(45)
#######################################################

upper_bound_index_cusps_in_JG_torsion(GammaH(45,[ 1, 4, 16, 19, 31, 34 ]),rational_cuspidal_classgroup(GammaH(45,[ 1, 4, 16, 19, 31, 34 ])).cardinality())

#######################################################
# Output
#######################################################

1


#########################################################################
# Proof of Proposition 4.16 
# Cuspidal torsion, Hecke bounds, and additional arguments on J1(2,2N)
# Note that X_{Delta}(4N) \cong X1(2,2N), with Delta := {±1, ±(2N + 1)} 
#########################################################################

G210 = GammaH(20,[1,-1,11,-11]);
J210 = rational_cuspidal_classgroup(G210)
J210.invariants()  # (6)
d = J210.cardinality() # 6
cuspidal_rational_subgroup_mod_rational_cuspidal_subgroup(G210) 
upper_bound_index_cusps_in_JG_torsion(G210,d) # 1

G212 = GammaH(24,[1,-1,13,-13]);
J212 = rational_cuspidal_classgroup(G212)
J212.invariants()  # (4)
d = J212.cardinality() # 4
upper_bound_index_cusps_in_JG_torsion(G212,d) #1

G214 = GammaH(28,[1,-1,15,-15]);
J214 = rational_cuspidal_classgroup(G214)
J214.invariants()  # (2, 2, 6, 18)
d = J214.cardinality() # 432
upper_bound_index_cusps_in_JG_torsion(G214,d) # 1

G216 = GammaH(32,[1,-1,17,-17]);
J216 = rational_cuspidal_classgroup(G216)
J216.invariants()  # (2, 20, 20)
d = J216.cardinality() # 800
upper_bound_index_cusps_in_JG_torsion(G216,d) # 1

G218 = GammaH(36,[1,-1,19,-19]);
J218 = rational_cuspidal_classgroup(G218)
J218.invariants()  # (2, 42, 126)
d = J218.cardinality() # 10584
upper_bound_index_cusps_in_JG_torsion(G218,d) # 1

G220 = GammaH(40,[1,-1,21,-21]);
J220 = rational_cuspidal_classgroup(G220)
J220.invariants()  # (4, 60, 120)
d = J220.cardinality() # 28800
upper_bound_index_cusps_in_JG_torsion(G220,d) # 1

G222 = GammaH(44,[1,-1,23,-23]);
J222 = rational_cuspidal_classgroup(G222)
J222.invariants()  # (2, 10, 1550, 4650)
d = J222.cardinality() # 144150000
upper_bound_index_cusps_in_JG_torsion(G222,d) # 1

G224 = GammaH(48,[1,-1,25,-25]);
J224 = rational_cuspidal_classgroup(G224)
J224.invariants()  # (2, 4, 4, 120, 240)
J224.cardinality() # 921600
d = upper_bound_index_cusps_in_JG_torsion(G224,d) # 4
do_cusps_generate_J_H_torsion(G224,bound = 24) # (True,(2, 4, 4, 120, 240))

G226 = GammaH(52,[1,-1,27,-27]);
J226 = rational_cuspidal_classgroup(G226)
J226.invariants()  # (2, 14, 266, 3990, 11970)
d = J226.cardinality() # 355718714400
upper_bound_index_cusps_in_JG_torsion(G226,d) # 1

G228 = GammaH(56,[1,-1,29,-29]);
J228 = rational_cuspidal_classgroup(G228)
J228.invariants()  # (2, 4, 4, 4, 8, 24, 936, 936)
d = J228.cardinality() # 21530935296
upper_bound_index_cusps_in_JG_torsion(G228,d) # 1

G230 = GammaH(60,[1,-1,31,-31]);
J230 = rational_cuspidal_classgroup(G230)
J230.invariants()  # (2, 2, 2, 2, 2, 24, 8160, 8160)
d = J230.cardinality() # 51137740800
upper_bound_index_cusps_in_JG_torsion(G230,d) # 2
do_cusps_generate_J_H_torsion(G230,bound = 24) # (True, (2, 2, 2, 2, 2, 24, 8160, 8160))

G232 = GammaH(64,[1,-1,33,-33]);
J232 = rational_cuspidal_classgroup(G232)
J232.invariants()  # (2, 2, 2, 4, 4, 24, 120, 23280, 23280)
d = J232.cardinality() # 199787544576000
upper_bound_index_cusps_in_JG_torsion(G232,d) # 4
do_cusps_generate_J_H_torsion(G232,bound = 24) # (True, (2, 2, 2, 4, 4, 24, 120, 23280, 23280))

