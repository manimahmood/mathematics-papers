# Projector and elimination data

Exact finite data for the accompanying manuscript.

## S2. Exact six-projection seed and whitening

The columns of the following 4-by-6 matrices are u_j and v_j:

    U = [
      [1, 0, 0, 0,  1, 1],
      [0, 1, 0, 0,  0, 1],
      [0, 0, 1, 0, -4, 2],
      [0, 0, 0, 1, -3, 3]
    ]

    V = [
      [ 0,  5, -1, -1, -17, -5],
      [-5,  0,  0,  5,   4,  4],
      [ 2, -2,  0, -3,  -2,  2],
      [-2,  1,  0,  0,  -3, -1]
    ]

For j=1,...,6,

    P_j = u_j u_j^T/(u_j^T u_j) + v_j v_j^T/(v_j^T v_j),
    S   = sum_j P_j.

Every u_j^T v_j is zero. Hence every P_j is symmetric, idempotent,
positive semidefinite, and has rank two. The exact factorization S=L D L^T
uses

    L = [
      [1,                   0,                       0,                    0],
      [-1206036/7353455,    1,                       0,                    0],
      [-629974/7353455,    -620743491/3418884916,    1,                    0],
      [865837/7353455,      562593553/3418884916,     39888940419/121436249827, 1]
    ]

    diag(D) = [
      1470691/332787,
      3418884916/1213320075,
      121436249827/51283273740,
      213445555084/121436249827
    ].

An independent positive-definiteness check used the four leading Sylvester
minors

    1470691/332787
    3418884916/274549275
    5279836949/179053875
    213445555084/4118239125

and found each strictly positive. Define

    C_w = D^(-1/2) L^(-1),
    B_j = C_w P_j C_w^T.

Then sum_j B_j=I_4 exactly. The order D^(-1/2)L^(-1) is part of the
construction; swapping it is a rejected mutation.

## S3. Trace split and five fixed functionals

For j=1,...,5, set

    alpha_j = tr(B_j)/4,
    D_j     = B_j - alpha_j I_4.

The exact offsets, in order, are

    alpha_1 = 62457323133/426891110168
    alpha_2 = 135091765759/853782220336
    alpha_3 = 79415386985/426891110168
    alpha_4 = 23248840457/106722777542
    alpha_5 = 7361309426/53361388771

Each D_j is trace-free. The inverse dictionary in
lem:quaternion-dictionary defines the exact 3-by-3 coefficient matrix K_j
without a numerical fit. For a unit quaternion q and C=R(q),

    y_j = q^T B_j q = alpha_j + <K_j,C>_F,
    y_6 = 1 - sum_{j=1}^5 y_j.

The exact determinants of K_1,...,K_5 are

    det(K_1) = 3473223797711049750536848194297
               /77794936864983961888480902470661632

    det(K_2) = 167547333464197895643225749749879
               /622359494919871695107847219765293056

    det(K_3) = 69943244850495814877402080561925
               /77794936864983961888480902470661632

    det(K_4) = 4706840730252752914553772389333
               /4862183554061497618030056404416352

    det(K_5) = 34895673885154660754615970997
               /303886472128843601126878525276022

All five determinants are positive and all five coefficient matrices have
rank three. Their rank sum 15 is only an algebraic rank-one coefficient-term
upper bound. It is not a primitive Green-shot count, multiplexing theorem,
device cost, or implementation optimum.

The large multiquadratic expanded entries of the five K_j are intentionally
not duplicated here. Equations eq:exact-uv through eq:functionals are the
exact constructive definition. Readers can reconstruct the entries from
the supplied projector columns and exact whitening factors; no external
expanded-entry file is required by that definition.

## S4. Measurement-kernel coordinates and exact basis

Use the symmetric-coordinate order

    [q00,q01,q02,q03,q11,q12,q13,q22,q23,q33].

The exact 6-by-10 measurement map has rank 6 and kernel dimension 4. In the
coordinate order above, its row-reduced kernel basis is

    Q_a = [
      0,
      791515/238017,
      -568640/238017,
      3294920/238017,
      -1477400/238017,
      -1608733/238017,
      1,0,0,0
    ]

    Q_b = [
      -1,
      -23164577/3570255,
      2506819/1428102,
      -177432509/7140510,
      7249859/714051,
      40134749/3570255,
      0,1,0,0
    ]

    Q_c = [
      0,
      1384846/396695,
      66939/79339,
      2527256/396695,
      -566396/396695,
      -866673/396695,
      0,0,1,0
    ]

    Q_d = [
      0,
      -5944223/3570255,
      1028125/1428102,
      -76065521/7140510,
      2854424/714051,
      18554201/3570255,
      0,0,0,1
    ].

Thus a general kernel element is Q(a,b,c,d)=a Q_a+b Q_b+c Q_c+d Q_d.
A symmetric 4-by-4 matrix has rank at most two exactly when its adjugate is
zero, equivalently when its ten upper-triangular cofactors vanish. Those ten
cofactors are the complete determinantal system used here.

## S5. Projective elimination and Sturm certificate

Homogeneity gives the four-patch projective cover

    d=1
    d=0,c=1
    d=0,c=0,b=1
    d=0,c=0,b=0,a=1.

On d=1, lexicographic elimination has basis length three and produces a
degree-ten polynomial f(c). With coefficients listed in descending degree,

    f(c) = a10*c^10 + a9*c^9 + ... + a1*c + a0,

where

    a10 = 519966562263643554945384191703616165395119112637573956248783182163623419737152512
    a9  = 328579249789044588040378180884308246920612283125532941525689039614516382331655616
    a8  = -1488937659336445964244382640314269130820570919591539086664937634093065289567722168
    a7  = -1233940048680917718690405606336150050552321029972511785557378326694535160186545716
    a6  = 560830252887171704842938122056614428129006430576409737097486332815580249915632343
    a5  = 862820282355455834964156023668088334391251839891193308388567471700918807326409512
    a4  = 775754320104988082038883422728415511210018447819392679061447643613437650180479290
    a3  = 470277354383117587803315722906687832056957450220672159363137099743786021781725000
    a2  = 200124465435786576051259835347973959273830283051981724273474893893577713039639000
    a1  = 61607571777035859344852037093432432693206774944581030976959966534882027371910000
    a0  = 8636626929016108140668606241999544716256996255890976455134579423688060066275000.

The signed-remainder Sturm chain has degrees

    [10,9,8,7,6,5,4,3,2,1,0].

Its signs at negative infinity are

    [1,-1,1,-1,-1,-1,-1,1,1,-1,-1],

and its signs at positive infinity are

    [1,1,1,1,-1,1,-1,-1,1,1,-1].

Both sign-variation counts are 5, so f has zero real roots. On the next two
patches the exact Groebner bases are respectively [1] and [1]. On the final
a=1 patch, clear denominators from each cofactor polynomial, divide its
integer coefficient content, and fix its sign. This nonzero rational
rescaling preserves the common determinantal zero locus. Seven primitive
cofactor polynomials are nonzero there, and the first declared nonzero
primitive denominator-cleared cofactor polynomial evaluates to

    2324385247266000.

Therefore the projective cofactor system has no real point and the symmetric
measurement kernel contains no nonzero real rank-at-most-two matrix.

The exact certificate identities are

    cofactor_polynomial_sha256 = 81cfd7c8a004cd667a60c6df3e8f397ff22c32b43350ebbad6eba76705b8a49c
    canonical_d1_sturm_patch_sha256 = 071f38cebd57219fbef4c6fd43f850cc6caf327af1de46fba4f3eb55326395f1
    certificate_sha256 = 580b8fa19de802dad2c3fcfb15019a80128d2028f27d9d395c6b6535cfda7be2

The d=1 hash is the SHA-256 of canonical compact sorted JSON for the complete
G3.d_nonzero_patch object. The certificate route is
RREF_KERNEL__TEN_COFACTORS__PROJECTIVE_GROEBNER__STURM. No imported Maple
verdict is consumed.
