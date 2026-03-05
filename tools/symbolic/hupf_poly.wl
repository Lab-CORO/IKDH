(* ================================================================ *)
(*  IKDH — HuPf Symbolic Polynomial Reduction                       *)
(*  Robot : ABB GoFa CRB 15000-5                                    *)
(*                                                                   *)
(*  Pipeline :                                                       *)
(*   1. Hyperplan matrix H (8x8)                                    *)
(*   2. Kinematic surface r[0..7] via 7x7 determinants              *)
(*   3. e1 = Study quadric,  e2 = 8th hyperplane constraint         *)
(*   4. realPol = Resultant[e1, e2, u]   (eliminates left variable) *)
(*   5. Factor[realPol]  ->  extract degree-16 factor in v          *)
(*   6. Export coefficient expressions as C++ / Python code         *)
(*                                                                   *)
(*  Usage :  Open in Mathematica, then Evaluation > Run All         *)
(*           Or from terminal :                                      *)
(*             math -script hupf_poly.wl                            *)
(* ================================================================ *)

Print["=== IKDH HuPf Symbolic Reduction — GoFa5 ==="];
Print["Started: ", DateString[]];

(* Set working directory to the folder containing this script *)
(* Works both from command line (math -script) and Mathematica notebook *)
SetDirectory[If[ValueQ[$InputFileName] && $InputFileName =!= "",
  DirectoryName[$InputFileName],
  NotebookDirectory[]]];
If[!DirectoryQ["checkpoints"], CreateDirectory["checkpoints"]];
Print["Working directory: ", Directory[]];

(* ---------------------------------------------------------------- *)
(* 1. DH PARAMETERS (exact rationals)                               *)
(*    al_i = tan(alpha_i / 2)                                       *)
(*    alpha = {-Pi/2, 0, -Pi/2, Pi/2, -Pi/2, 0}                    *)
(* ---------------------------------------------------------------- *)

al1 = -1;  al2 = 0;  al3 = -1;
al4 =  1;  al5 = -1; al6 =  0;

a1 = 0;    a2 = 444/1000;  a3 = 110/1000;
a4 = 0;    a5 =  80/1000;  a6 = 0;

d2 = 0;    d3 = 0;
d4 = 470/1000;  d5 = 0;  d6 = 101/1000;

(* ---------------------------------------------------------------- *)
(* 2. SYMBOLIC VARIABLES                                            *)
(*    t0..t7 : Study parameters of TCP pose                         *)
(*    u      : left  chain half-angle  v2 = tan(theta2/2)          *)
(*    v      : right chain half-angle  v4 = tan(theta4/2)          *)
(* ---------------------------------------------------------------- *)

(* Declare as symbols — Mathematica treats undefined names as symbols *)
(* t0,t1,t2,t3,t4,t5,t6,t7,u,v are automatically symbolic          *)

(* ---------------------------------------------------------------- *)
(* 3. LEFT CHAIN HYPERPLANES : h_tc_v2                              *)
(*    Polynomials in u with rational coefficients.                  *)
(*    C++ Polynomial(c0, c1)  ->  c0 + c1*u                        *)
(* ---------------------------------------------------------------- *)

h1tcv2 = {
  (-2*al3*al1*d2 + 2*al3*al2*d2 - 2*al3*d3*al2 - 2*al3*d3*al1) +
  ( 2*al3*a2 + 2*al3*al1*al2*a2 - 2*al3*al1*a1*al2 - 2*al3*a1 + 2*a3*al1 - 2*al2*a3)*u,

  ( 2*al1*d2 - 2*al2*d2 + 2*d3*al2 + 2*d3*al1) +
  (-2*a2 - 2*al1*al2*a2 + 2*al1*a1*al2 + 2*a1 + 2*al3*a3*al1 - 2*a3*al3*al2)*u,

  (-2*a1 - 2*a2 + 2*a2*al1*al2 + 2*a1*al1*al2 - 2*a3*al3*al2 - 2*a3*al3*al1) +
  ( 2*d2*al1 + 2*d2*al2 + 2*d3*al1 - 2*d3*al2)*u,

  ( 2*al3*a1 + 2*al3*a2 - 2*al3*al1*al2*a2 - 2*al3*al1*al2*a1 - 2*al2*a3 - 2*al1*a3) +
  (-2*al1*d2*al3 - 2*d2*al3*al2 - 2*d3*al1*al3 + 2*al3*d3*al2)*u,

  0 + (-4*al3*(al1 - al2))*u,
  0 + ( 4*al1 - 4*al2)*u,
  (-4*al2 - 4*al1),
  ( 4*al3*(al1 + al2))
};

h2tcv2 = {
  (-2*d2 - 2*al1*al2*d2 + 2*d3*al1*al2 - 2*d3) +
  ( 2*al2*a2 + 2*al1*a1 - 2*al1*a2 - 2*al2*a1 - 2*al3*a3*al1*al2 - 2*al3*a3)*u,

  (-2*al3*d2 - 2*al3*al1*al2*d2 + 2*al3*d3*al1*al2 - 2*al3*d3) +
  ( 2*a2*al3*al2 + 2*al3*al1*a1 - 2*al3*al1*a2 - 2*a1*al3*al2 + 2*a3*al1*al2 + 2*a3)*u,

  ( 2*a2*al3*al1 + 2*a1*al3*al1 + 2*a1*al3*al2 + 2*a2*al3*al2 - 2*a3*al1*al2 + 2*a3) +
  ( 2*d2*al3 - 2*al1*d2*al3*al2 + 2*al3*d3*al1*al2 + 2*d3*al3)*u,

  ( 2*al1*a2 + 2*al1*a1 + 2*al2*a1 + 2*al2*a2 + 2*al3*al1*al2*a3 - 2*al3*a3) +
  ( 2*d2 - 2*al1*d2*al2 + 2*d3*al1*al2 + 2*d3)*u,

  0 + (-4*al1*al2 - 4)*u,
  0 + (-4*al3*(1 + al1*al2))*u,
  4*al3*(-1 + al1*al2),
  4*al1*al2 - 4
};

h3tcv2 = {
  ( 2*al3*a1 + 2*al3*a2 - 2*al3*al1*al2*a2 - 2*al3*al1*al2*a1 - 2*al2*a3 - 2*al1*a3) +
  (-2*al1*d2*al3 - 2*d2*al3*al2 - 2*d3*al1*al3 + 2*al3*d3*al2)*u,

  (-2*a1 - 2*a2 + 2*a2*al1*al2 + 2*a1*al1*al2 - 2*a3*al3*al2 - 2*a3*al3*al1) +
  ( 2*d2*al1 + 2*d2*al2 + 2*d3*al1 - 2*d3*al2)*u,

  ( 2*al2*d2 - 2*al1*d2 - 2*d3*al2 - 2*d3*al1) +
  ( 2*a2 - 2*a1 + 2*al1*al2*a2 - 2*al1*a1*al2 + 2*a3*al3*al2 - 2*al3*a3*al1)*u,

  (-2*al3*al2*d2 + 2*al3*al1*d2 + 2*al3*d3*al2 + 2*al3*d3*al1) +
  (-2*al3*a2 + 2*al3*a1 - 2*al3*al1*al2*a2 + 2*al3*al1*a1*al2 + 2*al2*a3 - 2*a3*al1)*u,

  4*al3*(al1 + al2),
  (-4*al2 - 4*al1),
  0 + ( 4*al2 - 4*al1)*u,
  0 + ( 4*al3*(al1 - al2))*u
};

h4tcv2 = {
  (-2*al2*a1 - 2*al1*a1 - 2*al2*a2 - 2*al1*a2 + 2*al3*a3 - 2*al3*al1*al2*a3) +
  (-2*d2 + 2*al1*d2*al2 - 2*d3*al1*al2 - 2*d3)*u,

  (-2*a1*al3*al2 - 2*a1*al3*al1 - 2*a2*al3*al2 - 2*a2*al3*al1 - 2*a3 + 2*a3*al1*al2) +
  (-2*d2*al3 + 2*al1*d2*al3*al2 - 2*al3*d3*al1*al2 - 2*d3*al3)*u,

  (-2*al3*d2 - 2*al3*al1*al2*d2 + 2*al3*d3*al1*al2 - 2*al3*d3) +
  ( 2*a2*al3*al2 + 2*al3*al1*a1 - 2*al3*al1*a2 - 2*a1*al3*al2 + 2*a3*al1*al2 + 2*a3)*u,

  (-2*d2 - 2*al1*al2*d2 + 2*d3*al1*al2 - 2*d3) +
  ( 2*al2*a2 + 2*al1*a1 - 2*al1*a2 - 2*al2*a1 - 2*al3*a3*al1*al2 - 2*al3*a3)*u,

  4 - 4*al1*al2,
  -(4*al3*(-1 + al1*al2)),
  0 + (-(4*al3*(1 + al1*al2)))*u,
  0 + (-4*al1*al2 - 4)*u
};

(* ---------------------------------------------------------------- *)
(* 4. RIGHT CHAIN HYPERPLANES : h_v4q                              *)
(*    Polynomials in v with coefficients linear in t0..t7.          *)
(*    C++ Polynomial(c0_expr, c1_expr)  ->  c0_expr + c1_expr*v    *)
(* ---------------------------------------------------------------- *)

h1v4q = {
  (-t0*al6*a6 - 2*t4 + t2*al4*d5 - 2*t4*al4*al6 + t3*d4 + t3*d6 + t0*al4*a6
   - t2*al4*d4 - t0*al4*a4 + t0*a5*al5 + t0*al6*a4 + t2*al6*d5 + t2*al6*d4
   + t2*al6*d6 + t2*d6*al4 - t1*a5*al5*al4 - t1*al6*al4*a4 + t1*al6*a5*al5
   + t1*al6*a6*al4 - t3*al6*al4*d5 + t3*al6*al4*d4 - t3*al6*al4*d6
   - t1*a4 + t1*a6 + t3*d5 + t0*al6*a5*al5*al4 + 2*t5*al4 - 2*t5*al6) +
  (t1*al6*d4 - t3*a5*al5 - 2*t6*al4 + 2*t7 + t0*d5 - t3*al6*a5*al5*al4
   - t2*al6*a6*al4 + t2*a4 - t0*al6*al4*d5 + t1*al4*d5 + t1*al6*d5
   + t1*al6*d6 + t3*al4*a4 - t3*al6*a4 + t3*al6*a6 - t3*a6*al4
   + 2*t7*al6*al4 + t2*a5*al5*al4 - t2*al6*a5*al5 + t2*al6*al4*a4
   - t1*al4*d4 + t1*al4*d6 + t0*al6*al4*d4 - t0*al6*al4*d6
   + 2*t6*al6 + t0*d6 + t0*d4 - t2*a6)*v,

  (t1*a5*al5 + t1*al4*a6 - 2*t5 - 2*t4*al4 - t2*d6 - t1*al6*a6
   - t0*al6*a5*al5 + t3*al6*d5 - 2*t5*al4*al6 - t3*al4*d4 + t3*al6*d4
   + t3*al4*d5 + t0*a5*al5*al4 + t0*al6*al4*a4 - t0*al6*a6*al4
   + t2*al6*al4*d5 - t2*al6*al4*d4 - t1*al4*a4 + t2*al6*al4*d6
   + t1*al6*a4 + t3*al6*d6 + t3*d6*al4 + t0*a4 - t0*a6 - t2*d5
   - t2*d4 + 2*t4*al6 + t1*al6*a5*al5*al4) +
  (-t3*al6*a5*al5 - 2*t6 - t2*al6*a6 - t3*a6 + 2*t7*al6 + t2*a5*al5
   + t3*a4 + t1*d5 + t2*al6*a5*al5*al4 - t0*al4*d5 - t0*al6*d5
   - t0*al6*d6 - t2*al4*a4 + t2*al6*a4 + t2*a6*al4 - 2*t6*al6*al4
   - 2*t7*al4 - t1*al6*al4*d5 + t3*a5*al5*al4 + t3*al6*al4*a4
   - t3*al6*a6*al4 + t1*d6 + t1*d4 + t1*al6*al4*d4 - t1*al6*al4*d6
   + t0*al4*d4 - t0*al6*d4 - t0*al4*d6)*v,

  (-t3*al6*a5*al5 - 2*t6 - t2*al6*a6 - t3*a6 + 2*t7*al6 + t2*a5*al5
   + t3*a4 + t1*d5 + t2*al6*a5*al5*al4 - t0*al4*d5 - t0*al6*d5
   - t0*al6*d6 - t2*al4*a4 + t2*al6*a4 + t2*a6*al4 - 2*t6*al6*al4
   - 2*t7*al4 - t1*al6*al4*d5 + t3*a5*al5*al4 + t3*al6*al4*a4
   - t3*al6*a6*al4 + t1*d6 + t1*d4 + t1*al6*al4*d4 - t1*al6*al4*d6
   + t0*al4*d4 - t0*al6*d4 - t0*al4*d6) +
  (-t1*a5*al5 - t1*al4*a6 + 2*t5 + 2*t4*al4 + t2*d6 + t1*al6*a6
   + t0*al6*a5*al5 - t3*al6*d5 + 2*t5*al4*al6 + t3*al4*d4 - t3*al6*d4
   - t3*al4*d5 - t0*a5*al5*al4 - t0*al6*al4*a4 + t0*al6*a6*al4
   - t2*al6*al4*d5 + t2*al6*al4*d4 + t1*al4*a4 - t2*al6*al4*d6
   - t1*al6*a4 - t3*al6*d6 - t3*d6*al4 - t0*a4 + t0*a6 + t2*d5
   + t2*d4 - 2*t4*al6 - t1*al6*a5*al5*al4)*v,

  (-t1*al6*d4 + t3*a5*al5 + 2*t6*al4 - 2*t7 - t0*d5 + t3*al6*a5*al5*al4
   + t2*al6*a6*al4 - t2*a4 + t0*al6*al4*d5 - t1*al4*d5 - t1*al6*d5
   - t1*al6*d6 - t3*al4*a4 + t3*al6*a4 - t3*al6*a6 + t3*a6*al4
   - 2*t7*al6*al4 - t2*a5*al5*al4 + t2*al6*a5*al5 - t2*al6*al4*a4
   + t1*al4*d4 - t1*al4*d6 - t0*al6*al4*d4 + t0*al6*al4*d6
   - 2*t6*al6 - t0*d6 - t0*d4 + t2*a6) +
  (-t0*al6*a6 - 2*t4 + t2*al4*d5 - 2*t4*al4*al6 + t3*d4 + t3*d6
   + t0*al4*a6 - t2*al4*d4 - t0*al4*a4 + t0*a5*al5 + t0*al6*a4
   + t2*al6*d5 + t2*al6*d4 + t2*al6*d6 + t2*d6*al4 - t1*a5*al5*al4
   - t1*al6*al4*a4 + t1*al6*a5*al5 + t1*al6*a6*al4 - t3*al6*al4*d5
   + t3*al6*al4*d4 - t3*al6*al4*d6 - t1*a4 + t1*a6 + t3*d5
   + t0*al6*a5*al5*al4 + 2*t5*al4 - 2*t5*al6)*v,

  (-2*t0 - 2*t0*al6*al4 + 2*t1*al4 - 2*t1*al6) +
  (-2*t2*al4 + 2*t2*al6 + 2*t3 + 2*t3*al6*al4)*v,

  (-2*t1 - 2*t1*al6*al4 - 2*t0*al4 + 2*t0*al6) +
  (-2*t3*al4 + 2*t3*al6 - 2*t2 - 2*t2*al6*al4)*v,

  (-2*t3*al4 + 2*t3*al6 - 2*t2 - 2*t2*al6*al4) +
  ( 2*t0*al4 - 2*t0*al6 + 2*t1 + 2*t1*al6*al4)*v,

  (-2*t3 - 2*t3*al6*al4 + 2*t2*al4 - 2*t2*al6) +
  (-2*t0 - 2*t0*al6*al4 + 2*t1*al4 - 2*t1*al6)*v
};

h2v4q = {
  (2*t4*al5*al6 - t1*a5 - 2*t5*al5 - t3*al5*al4*d5 + t1*al6*al5*a4
   - t1*al6*a6*al5 - t1*al5*al4*a4 + t3*al5*al4*d4 - t2*al6*al5*al4*d5
   + t2*al6*al5*al4*d4 - t0*al6*al5*al4*a6 + t2*al5*d5
   + t2*al6*al5*al4*d6 - t3*al6*al5*d5 + t1*a6*al5*al4
   + t0*al6*al5*al4*a4 - 2*t5*al6*al5*al4 - t3*al6*al5*d4
   + t3*d6*al5*al4 + t0*al6*a5 - t1*al6*a5*al4 + t3*al6*al5*d6
   - 2*t4*al5*al4 - t0*al5*a6 + t2*al5*d4 - t2*al5*d6
   - t0*a5*al4 + t0*al5*a4) +
  (2*t6*al5 + t3*a6*al5 + 2*t7*al5*al4 - t3*al5*a4 - 2*t7*al6*al5
   - t1*al5*d6 + t3*a5*al4 - t0*al6*al5*d5 + t1*al6*d6*al5*al4
   - t3*al6*a5 + t2*al5*al4*a4 - t0*al6*al5*d4 + t0*al5*al4*d6
   - t2*al5*al6*a6 + t2*al5*al4*a6 + 2*t6*al5*al4*al6
   + t0*al5*al4*d4 + t2*a5 + t0*al6*al5*d6 - t0*al5*al4*d5
   + t3*al6*a6*al5*al4 - t1*al6*al5*al4*d5 + t1*al5*d5
   + t1*al5*d4 + t2*al6*a5*al4 + t1*al6*al5*al4*d4
   - t3*al6*al5*al4*a4 - t2*al6*al5*a4)*v,

  (2*t4*al5 + t3*al5*d4 + t1*al6*a5 - t1*a6*al5 + t0*a5
   + 2*t5*al6*al5 - t2*al6*al5*d6 - t2*al5*d6*al4 - t3*d6*al5
   + t2*al6*al5*d5 - t3*al6*al5*al4*d5 + t0*al5*al4*a4
   - t0*al6*al5*a4 - t1*a5*al4 + t3*al6*al5*al4*d4 + t0*al6*a5*al4
   + t0*al6*a6*al5 + t1*al5*a4 + t2*al6*al5*d4
   + 2*t4*al5*al6*al4 + t3*al6*al5*al4*d6 - t0*al5*a6*al4
   - 2*t5*al5*al4 + t3*al5*d5 + t1*al6*al5*al4*a4
   + t2*al5*al4*d5 - t2*al5*al4*d4 - t1*al6*al5*al4*a6) +
  (t3*a5 + t1*al6*d6*al5 + t2*al6*a5 + 2*t7*al5 - 2*t6*al5*al4
   + t3*al6*a5*al4 + t3*al5*al4*a4 + 2*t6*al6*al5
   - t2*al5*al6*a6*al4 - t0*al5*d4 + t0*al5*d6
   + t2*al6*al5*al4*a4 - t1*al5*al4*d5 + t0*al6*al5*al4*d5
   - t3*al6*al5*a4 - t0*al6*al5*d6*al4 - t1*al6*al5*d5
   - t2*a6*al5 - t0*al6*al5*al4*d4 + t3*al5*al6*a6
   - t3*al5*al4*a6 + t1*al5*al4*d6 + 2*t7*al5*al4*al6
   - t1*al5*al6*d4 + t1*al5*al4*d4 - t0*al5*d5
   - t2*a5*al4 + t2*al5*a4)*v,

  (t3*a5 + t1*al6*d6*al5 + t2*al6*a5 + 2*t7*al5 - 2*t6*al5*al4
   + t3*al6*a5*al4 + t3*al5*al4*a4 + 2*t6*al6*al5
   - t2*al5*al6*a6*al4 - t0*al5*d4 + t0*al5*d6
   + t2*al6*al5*al4*a4 - t1*al5*al4*d5 + t0*al6*al5*al4*d5
   - t3*al6*al5*a4 - t0*al6*al5*d6*al4 - t1*al6*al5*d5
   - t2*a6*al5 - t0*al6*al5*al4*d4 + t3*al5*al6*a6
   - t3*al5*al4*a6 + t1*al5*al4*d6 + 2*t7*al5*al4*al6
   - t1*al5*al6*d4 + t1*al5*al4*d4 - t0*al5*d5
   - t2*a5*al4 + t2*al5*a4) +
  (-2*t4*al5 - t3*al5*d4 - t1*al6*a5 + t1*a6*al5 - t0*a5
   - 2*t5*al6*al5 + t2*al6*al5*d6 + t2*al5*d6*al4 + t3*d6*al5
   - t2*al6*al5*d5 + t3*al6*al5*al4*d5 - t0*al5*al4*a4
   + t0*al6*al5*a4 + t1*a5*al4 - t3*al6*al5*al4*d4 - t0*al6*a5*al4
   - t0*al6*a6*al5 - t1*al5*a4 - t2*al6*al5*d4
   - 2*t4*al5*al6*al4 - t3*al6*al5*al4*d6 + t0*al5*a6*al4
   + 2*t5*al5*al4 - t3*al5*d5 - t1*al6*al5*al4*a4
   - t2*al5*al4*d5 + t2*al5*al4*d4 + t1*al6*al5*al4*a6)*v,

  (-2*t6*al5 - t3*a6*al5 - 2*t7*al5*al4 + t3*al5*a4 + 2*t7*al6*al5
   + t1*al5*d6 - t3*a5*al4 + t0*al6*al5*d5 - t1*al6*d6*al5*al4
   + t3*al6*a5 - t2*al5*al4*a4 + t0*al6*al5*d4 - t0*al5*al4*d6
   - t2*al5*al6*a6 + t2*al5*al4*a6 - 2*t6*al5*al4*al6
   - t0*al5*al4*d4 - t2*a5 - t0*al6*al5*d6 + t0*al5*al4*d5
   - t3*al6*a6*al5*al4 + t1*al6*al5*al4*d5 - t1*al5*d5
   - t1*al5*d4 - t2*al6*a5*al4 - t1*al6*al5*al4*d4
   + t3*al6*al5*al4*a4 + t2*al6*al5*a4) +
  (2*t4*al5*al6 - t1*a5 - 2*t5*al5 - t3*al5*al4*d5 + t1*al6*al5*a4
   - t1*al6*a6*al5 - t1*al5*al4*a4 + t3*al5*al4*d4
   - t2*al6*al5*al4*d5 + t2*al6*al5*al4*d4 - t0*al6*al5*al4*a6
   + t2*al5*d5 + t2*al6*al5*al4*d6 - t3*al6*al5*d5
   + t1*a6*al5*al4 + t0*al6*al5*al4*a4 - 2*t5*al6*al5*al4
   - t3*al6*al5*d4 + t3*d6*al5*al4 + t0*al6*a5 - t1*al6*a5*al4
   + t3*al6*al5*d6 - 2*t4*al5*al4 - t0*al5*a6 + t2*al5*d4
   - t2*al5*d6 - t0*a5*al4 + t0*al5*a4)*v,

  (-2*t0*al5*al4 + 2*t0*al6*al5 - 2*t1*al5 - 2*t1*al6*al5*al4) +
  ( 2*t2*al5 + 2*t2*al6*al5*al4 + 2*t3*al5*al4 - 2*t3*al6*al5)*v,

  (-2*t1*al5*al4 + 2*t1*al6*al5 + 2*t0*al5 + 2*t0*al6*al5*al4) +
  ( 2*t3*al5 + 2*t3*al6*al5*al4 - 2*t2*al5*al4 + 2*t2*al6*al5)*v,

  ( 2*t3*al5 + 2*t3*al6*al5*al4 - 2*t2*al5*al4 + 2*t2*al6*al5) +
  (-2*t0*al5 - 2*t0*al6*al5*al4 + 2*t1*al5*al4 - 2*t1*al6*al5)*v,

  (-2*t3*al5*al4 + 2*t3*al6*al5 - 2*t2*al5 - 2*t2*al6*al5*al4) +
  (-2*t0*al5*al4 + 2*t0*al6*al5 - 2*t1*al5 - 2*t1*al6*al5*al4)*v
};

h3v4q = {
  (-2*t6*al5 - t3*a6*al5 + 2*t7*al5*al4 - t3*al5*a4 + 2*t7*al6*al5
   + t1*al5*d6 + t3*a5*al4 + t0*al6*al5*d5 + t1*al6*d6*al5*al4
   + t3*al6*a5 - t2*al5*al4*a4 + t0*al6*al5*d4 + t0*al5*al4*d6
   - t2*al5*al6*a6 - t2*al5*al4*a6 + 2*t6*al5*al4*al6
   + t0*al5*al4*d4 - t2*a5 - t0*al6*al5*d6 - t0*al5*al4*d5
   + t3*al6*a6*al5*al4 - t1*al6*al5*al4*d5 - t1*al5*d5
   - t1*al5*d4 + t2*al6*a5*al4 + t1*al6*al5*al4*d4
   + t3*al6*al5*al4*a4 - t2*al6*al5*a4) +
  (2*t4*al5*al6 - t1*a5 - 2*t5*al5 + t3*al5*al4*d5 - t1*al6*al5*a4
   - t1*al6*a6*al5 - t1*al5*al4*a4 - t3*al5*al4*d4
   + t2*al6*al5*al4*d5 - t2*al6*al5*al4*d4 + t0*al6*al5*al4*a6
   + t2*al5*d5 - t2*al6*al5*al4*d6 - t3*al6*al5*d5
   - t1*a6*al5*al4 + t0*al6*al5*al4*a4 + 2*t5*al6*al5*al4
   - t3*al6*al5*d4 - t3*d6*al5*al4 + t0*al6*a5 + t1*al6*a5*al4
   + t3*al6*al5*d6 + 2*t4*al5*al4 - t0*al5*a6 + t2*al5*d4
   - t2*al5*d6 + t0*a5*al4 - t0*al5*a4)*v,

  (-t3*a5 - t1*al6*d6*al5 - t2*al6*a5 - 2*t7*al5 - 2*t6*al5*al4
   + t3*al6*a5*al4 - t3*al5*al4*a4 - 2*t6*al6*al5
   - t2*al5*al6*a6*al4 + t0*al5*d4 - t0*al5*d6
   - t2*al6*al5*al4*a4 - t1*al5*al4*d5 + t0*al6*al5*al4*d5
   - t3*al6*al5*a4 - t0*al6*al5*d6*al4 + t1*al6*al5*d5
   + t2*a6*al5 - t0*al6*al5*al4*d4 - t3*al5*al6*a6
   - t3*al5*al4*a6 + t1*al5*al4*d6 + 2*t7*al5*al4*al6
   + t1*al5*al6*d4 + t1*al5*al4*d4 + t0*al5*d5 - t2*a5*al4 + t2*al5*a4) +
  (2*t4*al5 + t3*al5*d4 + t1*al6*a5 - t1*a6*al5 + t0*a5
   + 2*t5*al6*al5 - t2*al6*al5*d6 + t2*al5*d6*al4 - t3*d6*al5
   + t2*al6*al5*d5 + t3*al6*al5*al4*d5 + t0*al5*al4*a4
   + t0*al6*al5*a4 + t1*a5*al4 - t3*al6*al5*al4*d4 - t0*al6*a5*al4
   + t0*al6*a6*al5 - t1*al5*a4 + t2*al6*al5*d4 - 2*t4*al5*al6*al4
   - t3*al6*al5*al4*d6 + t0*al5*a6*al4 + 2*t5*al5*al4 + t3*al5*d5
   + t1*al6*al5*al4*a4 - t2*al5*al4*d5 + t2*al5*al4*d4
   + t1*al6*al5*al4*a6)*v,

  (t3*al5*d4 + 2*t4*al5 + t2*al5*al4*d4 + t1*a5*al4 + t1*al6*al5*al4*a4
   - t1*al5*a6 + t3*al6*al5*al4*d5 + t1*al6*al5*al4*a6 - t1*al5*a4
   + t2*d6*al5*al4 + t1*al6*a5 - t2*al5*al4*d5 + t2*al6*al5*d5
   - t3*al6*al5*al4*d6 + 2*t5*al6*al5 + t0*a5 - t3*d6*al5
   + t0*al5*al4*a4 + t0*al6*al5*a4 - t3*al6*al5*al4*d4
   - t0*al6*a5*al4 + t0*al6*a6*al5 + t2*al6*al5*d4
   - t2*al6*al5*d6 - 2*t4*al6*al5*al4 + t0*al5*a6*al4
   + 2*t5*al5*al4 + t3*al5*d5) +
  (t3*a5 + t1*al6*d6*al5 + t2*al6*a5 + 2*t7*al5 + 2*t6*al5*al4
   - t3*al6*a5*al4 + t3*al5*al4*a4 + 2*t6*al6*al5
   + t2*al5*al6*a6*al4 - t0*al5*d4 + t0*al5*d6
   + t2*al6*al5*al4*a4 + t1*al5*al4*d5 - t0*al6*al5*al4*d5
   + t3*al6*al5*a4 + t0*al6*al5*d6*al4 - t1*al6*al5*d5
   - t2*a6*al5 + t0*al6*al5*al4*d4 + t3*al5*al6*a6 + t3*al5*al4*a6
   - t1*al5*al4*d6 - 2*t7*al5*al4*al6 - t1*al5*al6*d4
   - t1*al5*al4*d4 - t0*al5*d5 + t2*a5*al4 - t2*al5*a4)*v,

  (-t0*al6*a5 - t0*al6*al5*al4*a4 + t1*a5 + t2*d6*al5 + t1*al5*al4*a4
   - t2*al5*d5 - 2*t5*al6*al5*al4 + 2*t5*al5 + t2*al6*al5*al4*d6
   + t2*al6*al5*al4*d4 + t3*al6*al5*d4 + t1*al5*a6*al4
   + t0*al5*a4 - t3*al5*al4*d5 - t2*al6*al5*al4*d5
   - t1*al6*a5*al4 + t3*al5*al4*d4 + t1*al6*al5*a4
   + t1*al6*a6*al5 - t2*al5*d4 + t3*d6*al5*al4 + t3*al6*al5*d5
   - 2*t4*al5*al4 - t0*a5*al4 + t0*al5*a6 - t3*al6*al5*d6
   - 2*t4*al6*al5 - t0*al6*al5*al4*a6) +
  (2*t7*al6*al5 + 2*t7*al5*al4 - t2*a5 - t2*al6*a6*al5
   + t1*al6*al5*al4*d4 + t3*al6*a5 - t0*al6*d6*al5 - 2*t6*al5
   - t1*al6*al5*al4*d5 + t3*a5*al4 - t1*al5*d5 + t1*al6*d6*al5*al4
   + t0*al5*al4*d4 + t0*al6*al5*d4 - t2*al5*al4*a6
   + t0*al5*al4*d6 + 2*t6*al5*al4*al6 - t2*al6*al5*a4
   + t2*al6*a5*al4 + t3*al5*al6*a6*al4 - t2*al5*al4*a4
   - t3*a6*al5 + t3*al6*al5*al4*a4 - t3*al5*a4
   - t0*al5*al4*d5 - t1*al5*d4 + t1*al5*d6 + t0*al6*al5*d5)*v,

  (-2*t2*al5 + 2*t2*al6*al5*al4 + 2*t3*al5*al4 + 2*t3*al6*al5) +
  ( 2*t0*al5*al4 + 2*t0*al6*al5 - 2*t1*al5 + 2*t1*al6*al5*al4)*v,

  (-2*t3*al5 + 2*t3*al6*al5*al4 - 2*t2*al5*al4 - 2*t2*al6*al5) +
  ( 2*t1*al5*al4 + 2*t1*al6*al5 + 2*t0*al5 - 2*t0*al6*al5*al4)*v,

  ( 2*t1*al5*al4 + 2*t1*al6*al5 + 2*t0*al5 - 2*t0*al6*al5*al4) +
  ( 2*t2*al5*al4 + 2*t2*al6*al5 + 2*t3*al5 - 2*t3*al6*al5*al4)*v,

  ( 2*t1*al5 - 2*t1*al6*al5*al4 - 2*t0*al5*al4 - 2*t0*al6*al5) +
  (-2*t2*al5 + 2*t2*al6*al5*al4 + 2*t3*al5*al4 + 2*t3*al6*al5)*v
};

h4v4q = {
  (2*t6*al6 + 2*t7 - t2*a6 + t1*al6*d6 + t0*d4 - t1*al4*d5 + t3*al6*a6
   - t2*a4 + t3*al6*a4 + t1*al4*d4 + 2*t6*al4 + t0*d6 + t3*al4*a6
   + t0*d5 + t0*al6*al4*d6 - t1*d6*al4 + t1*al6*d5 - t3*a5*al5
   + t0*al6*al4*d5 - 2*t7*al4*al6 + t1*al6*d4 - t0*al6*al4*d4
   + t2*al6*a6*al4 + t2*al6*al4*a4 - t2*al6*a5*al5 + t3*al4*a4
   + t3*al6*a5*al5*al4 - t2*a5*al5*al4) +
  (t2*al4*d5 - 2*t4*al6*al4 + t0*al4*a4 - t1*a6 + 2*t4 + t0*al6*a4
   - t0*a5*al5 - t3*d4 - t3*d5 - t3*d6 + t3*al6*al4*d4
   - t3*al6*al4*d6 - t2*al6*d5 + 2*t5*al6 + t1*al6*a6*al4
   + t2*al4*d6 - t2*al6*d4 - t2*al4*d4 + t1*al6*al4*a4
   - t1*a5*al5*al4 - t1*a4 + 2*t5*al4 + t0*al6*a5*al5*al4
   + t0*al6*a6 - t3*al6*al4*d5 + t0*a6*al4
   - t1*al6*a5*al5 - t2*al6*d6)*v,

  (-t2*al6*a6 + 2*t7*al6 - t0*al6*d6 - 2*t6 - t2*al6*a4 - t2*al4*a6
   - t0*al6*d4 - t2*al4*a4 + t1*al6*al4*d5 - t3*a6 + t2*a5*al5
   + t1*d5 - t1*al6*al4*d4 + t0*d6*al4 + t1*d4 - t0*al4*d4
   + 2*t6*al4*al6 - t0*al6*d5 - t3*a4 + t1*d6
   - t2*al6*a5*al5*al4 + t1*al6*al4*d6 - t3*a5*al5*al4
   - 2*t7*al4 - t3*al6*a5*al5 + t3*al6*al4*a4
   + t0*al4*d5 + t3*al6*a6*al4) +
  (2*t5 + t2*d6 - 2*t4*al6 + t1*al6*a5*al5*al4 - t3*al6*d4
   - t3*al4*d4 + t3*al4*d6 + t2*d5 + t0*a6 + t0*al6*a5*al5
   + t0*a5*al5*al4 - t3*al6*d5 - t3*al6*d6 + t1*al4*a4
   + t1*al6*a4 + t2*d4 - 2*t5*al6*al4 + t3*al4*d5
   - t0*al6*al4*a4 - t0*al6*a6*al4 - 2*t4*al4 + t0*a4
   + t1*al6*a6 + t2*al6*al4*d5 + t1*a6*al4 + t2*al6*al4*d6
   - t2*al6*al4*d4 - t1*a5*al5)*v,

  (2*t5 + t2*d6 - 2*t4*al6 + t1*al6*a5*al5*al4 - t3*al6*d4
   - t3*al4*d4 + t3*al4*d6 + t2*d5 + t0*a6 + t0*al6*a5*al5
   + t0*a5*al5*al4 - t3*al6*d5 - t3*al6*d6 + t1*al4*a4
   + t1*al6*a4 + t2*d4 - 2*t5*al6*al4 + t3*al4*d5
   - t0*al6*al4*a4 - t0*al6*a6*al4 - 2*t4*al4 + t0*a4
   + t1*al6*a6 + t2*al6*al4*d5 + t1*a6*al4 + t2*al6*al4*d6
   - t2*al6*al4*d4 - t1*a5*al5) +
  (t2*al6*a6 - 2*t7*al6 + t0*al6*d6 + 2*t6 + t2*al6*a4 + t2*al4*a6
   + t0*al6*d4 + t2*al4*a4 - t1*al6*al4*d5 + t3*a6 - t2*a5*al5
   - t1*d5 + t1*al6*al4*d4 - t0*d6*al4 - t1*d4 + t0*al4*d4
   - 2*t6*al4*al6 + t0*al6*d5 + t3*a4 - t1*d6 + t2*al6*a5*al5*al4
   - t1*al6*al4*d6 + t3*a5*al5*al4 + 2*t7*al4 + t3*al6*a5*al5
   - t3*al6*al4*a4 - t0*al4*d5 - t3*al6*a6*al4)*v,

  (-t2*al4*d5 + 2*t4*al6*al4 - t0*al4*a4 + t1*a6 - 2*t4 - t0*al6*a4
   + t0*a5*al5 + t3*d4 + t3*d5 + t3*d6 - t3*al6*al4*d4
   + t3*al6*al4*d6 + t2*al6*d5 - 2*t5*al6 - t1*al6*a6*al4
   - t2*al4*d6 + t2*al6*d4 + t2*al4*d4 - t1*al6*al4*a4
   + t1*a5*al5*al4 + t1*a4 - 2*t5*al4 - t0*al6*a5*al5*al4
   - t0*al6*a6 + t3*al6*al4*d5 - t0*a6*al4
   + t1*al6*a5*al5 + t2*al6*d6) +
  (2*t6*al6 + 2*t7 - t2*a6 + t1*al6*d6 + t0*d4 - t1*al4*d5 + t3*al6*a6
   - t2*a4 + t3*al6*a4 + t1*al4*d4 + 2*t6*al4 + t0*d6 + t3*al4*a6
   + t0*d5 + t0*al6*al4*d6 - t1*d6*al4 + t1*al6*d5 - t3*a5*al5
   + t0*al6*al4*d5 - 2*t7*al4*al6 + t1*al6*d4 - t0*al6*al4*d4
   + t2*al6*a6*al4 + t2*al6*al4*a4 - t2*al6*a5*al5 + t3*al4*a4
   + t3*al6*a5*al5*al4 - t2*a5*al5*al4)*v,

  (2*t2*al4 + 2*t2*al6 + 2*t3 - 2*t3*al6*al4) +
  (2*t0 - 2*t0*al6*al4 + 2*t1*al4 + 2*t1*al6)*v,

  (2*t3*al4 + 2*t3*al6 - 2*t2 + 2*t2*al6*al4) +
  (2*t1 - 2*t1*al6*al4 - 2*t0*al4 - 2*t0*al6)*v,

  (2*t1 - 2*t1*al6*al4 - 2*t0*al4 - 2*t0*al6) +
  (2*t2 - 2*t2*al6*al4 - 2*t3*al4 - 2*t3*al6)*v,

  (-2*t1*al4 - 2*t1*al6 - 2*t0 + 2*t0*al6*al4) +
  (2*t2*al4 + 2*t2*al6 + 2*t3 - 2*t3*al6*al4)*v
};

(* ---------------------------------------------------------------- *)
(* 5. ASSEMBLE 8x8 HYPERPLANE MATRIX                               *)
(*    H[[i]] = row i  (1-based in Mathematica)                     *)
(*    Rows 1-4 : left chain,  Rows 5-8 : right chain               *)
(* ---------------------------------------------------------------- *)

H = {h1tcv2, h2tcv2, h3tcv2, h4tcv2,
     h1v4q,  h2v4q,  h3v4q,  h4v4q};

(* Evaluate all DH numeric constants *)
H = Expand[H];

Print["[1] Hyperplane matrix H built (8x8). Expanding..."];

(* ---------------------------------------------------------------- *)
(* 6. KINEMATIC SURFACE r[0..7]  via  Cramer's rule                *)
(*    A[[i,j]] = H[[i, j+1]]   for i=1..7, j=1..7   (x0 = 1)     *)
(*    b[[i]]   = -H[[i, 1]]    (RHS)                               *)
(*    r[[1]]   = Det[A]         (common denominator)                *)
(*    r[[k+1]] = Det[A with col k replaced by b]  for k=1..7       *)
(* ---------------------------------------------------------------- *)

Print["\n[2] Computing kinematic surface r[1..8] via 7x7 determinants..."];

(* Build 7x7 matrix A and RHS b *)
A7 = Table[Expand[H[[i, j+1]]], {i, 1, 7}, {j, 1, 7}];
b7 = Table[Expand[-H[[i, 1]]], {i, 1, 7}];

(* Launch parallel kernels (one per physical core) *)
LaunchKernels[];
Print["  Parallel kernels available: ", $KernelCount];
DistributeDefinitions[A7, b7, u, v];

(* r[1] = det(A)  and  r[2..8] = Cramer numerators — all in parallel *)
{tAll, rVec} = Timing[
  ParallelTable[
    Module[{M = A7},
      If[k > 1,
        (* Replace column k-1 with b7 for Cramer numerators *)
        Do[M[[i, k-1]] = b7[[i]], {i, 1, 7}]
      ];
      Expand[Det[M]]
    ],
    {k, 1, 8}  (* k=1: plain det; k=2..8: Cramer for x1..x7 *)
  ]
];
Print["  All 8 determinants done in ", tAll, "s"];
Do[
  Print["  r[", k, "]: deg(u)=", Exponent[rVec[[k]], u],
        "  deg(v)=", Exponent[rVec[[k]], v]],
  {k, 1, 8}
];

(* Unpack: r1=rVec[[1]], ..., r8=rVec[[8]] *)
{r1,r2,r3,r4,r5,r6,r7,r8} = rVec;

(* ---------------------------------------------------------------- *)
(* 7. BUILD e1 AND e2                                               *)
(*    e1 = Study quadric:  r0*r4 + r1*r5 + r2*r6 + r3*r7          *)
(*         (0-indexed: rVec[[1]]*rVec[[5]] + ...)                   *)
(*    e2 = H[[8]] . rVec   (8th hyperplane constraint)              *)
(* ---------------------------------------------------------------- *)

Print["\n[3] Building e1 (Study quadric) and e2 (8th hyperplane)..."];

{te1, e1} = Timing[Expand[
  r1*r5 + r2*r6 + r3*r7 + r4*r8
]];
Print["  e1 in ", te1, "s  |  deg(u)=", Exponent[e1, u],
      "  deg(v)=", Exponent[e1, v]];

{te2, e2} = Timing[Expand[
  Sum[H[[8, j]] * rVec[[j]], {j, 1, 8}]
]];
Print["  e2 in ", te2, "s  |  deg(u)=", Exponent[e2, u],
      "  deg(v)=", Exponent[e2, v]];

(* ---------------------------------------------------------------- *)
(* 8. RESULTANT via PRS — with checkpoints at each step            *)
(*    Polynomial Remainder Sequence eliminates u step by step.     *)
(*    Each step is saved so computation can be resumed if stopped. *)
(*    Resume logic: loads last completed step automatically.        *)
(*                                                                  *)
(*    e1: deg(u)=8,  e2: deg(u)=4  ->  4-5 PRS steps expected     *)
(* ---------------------------------------------------------------- *)

Print["\n[4] Computing Resultant via PRS (with checkpoints)..."];

(* Save e1 and e2 as starting point *)
Export["checkpoints/prs_f0.mx", e1];
Export["checkpoints/prs_g0.mx", e2];
Print["  Saved starting polynomials prs_f0.mx, prs_g0.mx"];

(* Resume logic: find the last completed PRS step *)
resumeStep = 0;
While[FileExistsQ["checkpoints/prs_step" <> ToString[resumeStep+1] <> ".mx"],
  resumeStep++
];

If[resumeStep == 0,
  (* Start fresh *)
  prsF = e1;
  prsG = e2;
  Print["  Starting PRS from beginning"],
  (* Resume from last saved step *)
  Print["  Resuming from step ", resumeStep];
  If[resumeStep == 1,
    prsF = e2;
    prsG = Import["checkpoints/prs_step1.mx"],
    prsF = Import["checkpoints/prs_step" <> ToString[resumeStep-1] <> ".mx"];
    prsG = Import["checkpoints/prs_step" <> ToString[resumeStep] <> ".mx"]
  ]
];

(* PRS loop: each step reduces degree in u *)
prsStep = resumeStep;
While[
  Module[{degG = Exponent[prsG, u, -Infinity]},
    Print["\n  -- PRS step ", prsStep+1,
          ": deg(f,u)=", Exponent[prsF, u],
          "  deg(g,u)=", degG, " --"];
    degG > 0
  ],

  prsStep++;
  {tStep, prsR} = Timing[
    Expand[PolynomialRemainder[prsF, prsG, u]]
  ];
  degR = Exponent[prsR, u, -Infinity];
  Print["  Step ", prsStep, " done in ", tStep, "s",
        " | deg(rem,u)=", degR,
        " deg(rem,v)=", Exponent[prsR, v]];

  Export["checkpoints/prs_step" <> ToString[prsStep] <> ".mx", prsR];
  Print["  Saved checkpoints/prs_step", prsStep, ".mx"];

  prsF = prsG;
  prsG = prsR;
];

(* Last non-zero remainder = resultant (up to scalar in v,t0..t7) *)
realPol = Expand[prsG];
Print["\n  PRS complete after ", prsStep, " steps."];
Print["  realPol: deg(u)=", Exponent[realPol, u, -Infinity],
      "  deg(v)=", Exponent[realPol, v]];
Print["  Degree in v = ", Exponent[realPol, v]];

(* Save intermediate result *)
Export["checkpoints/realPol.mx", realPol];
Print["  Saved to checkpoints/realPol.mx"];

(* ---------------------------------------------------------------- *)
(* 9. FACTOR  ->  extract degree-16 factor in v                    *)
(* ---------------------------------------------------------------- *)

Print["\n[5] Factoring the resultant polynomial ..."];

{tFac, factors} = Timing[FactorList[realPol]];
Print["  Factoring done in ", tFac, "s"];
Print["  Number of factors: ", Length[factors]];

Do[
  {fac, mult} = factors[[i]];
  dv = Exponent[fac, v];
  Print["  Factor ", i, ": degree_v=", dv, "  multiplicity=", mult],
  {i, 1, Length[factors]}
];

(* ---------------------------------------------------------------- *)
(* 10. EXTRACT AND EXPORT degree-16 factor                         *)
(* ---------------------------------------------------------------- *)

Print["\n[6] Extracting degree-16 factor ..."];

poly16 = Select[factors, Exponent[#[[1]], v] == 16 &];

If[Length[poly16] == 0,
  Print["  WARNING: No degree-16 factor found!"],

  poly16 = poly16[[1, 1]];  (* the polynomial itself *)
  Print["  Found! degree_v=", Exponent[poly16, v]];

  (* Extract coefficients c[0]..c[16] where poly = Sum[c[k]*v^k] *)
  coeffs = CoefficientList[poly16, v];  (* coeffs[[k+1]] = c[k] *)
  Print["  Coefficients c[0]..c[16]:"];
  Do[
    c = Simplify[coeffs[[k+1]]];
    Print["    c[", k, "] = ", c],
    {k, 0, 16}
  ];

  (* Save *)
  Export["checkpoints/poly16.mx", poly16];
  Export["checkpoints/coeffs16.mx", coeffs];

  (* Export as text for use in C++ or Python *)
  outLines = {"(* Degree-16 IK polynomial for GoFa5 *)",
              "(* p(v) = Sum[c[k]*v^k, {k,0,16}] *)",
              "(* v = tan(theta4/2), t0..t7 = Study parameters of TCP *)",
              ""};
  Do[
    c = Simplify[coeffs[[k+1]]];
    AppendTo[outLines, "c[" <> ToString[k] <> "] = " <> ToString[c]],
    {k, 0, 16}
  ];
  Export["poly16_coeffs.txt", StringJoin[Riffle[outLines, "\n"]], "Text"];
  Print["  Coefficients saved to poly16_coeffs.txt"];
];

Print["\n=== Done! Finished: ", DateString[], " ==="];