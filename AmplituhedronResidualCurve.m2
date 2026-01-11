needsPackage "Resultants"

-- params for amplituhedron
n = 10;
k = 3;
m = 2;

----------- setup for ideal -------------
R = QQ[gens Grass(k-1,k+m-1)];
G = sub(ideal Grass(k-1,k+m-1),R);
Gr = R/G;

use R;

-- make n x (k+m) (= d+1) matrix
Z = transpose matrix {toList(n:1)} | matrix apply(toList(0..n-1), i -> apply(k+m-1, j -> (i)^(j+1)));

-- function for extracting Plueckers/minors of size n
plueckers = (M) -> (
    S := subsets(0..(numcols M)-1,numrows M)/toList;
    tmp := {};
    for i to #S-1 do (tmp=tmp|{det M_(S#i)});
    return tmp
);

-- chow forms
CH = (j) -> (
    M := Z^{j,(j+1)%n};
    sgns := reverse {1,-1,1,-1,1,-1,1,1,-1,1};
    pl := plueckers(M);
    G := reverse gens Grass(k-1,k+m-1);
    sub(sum for i to #G-1 list G#i*pl#i*sgns#i,R)
);


-- make residual curve
E = ideal(CH(0)) + ideal(CH(2)) + ideal(CH(4)) + ideal(CH(6)) + ideal(CH(8)) + G;

-- check irreducibility
# decompose E

-- check that is a complete intersection inside Grassmannian
codim sub(E,Gr), numgens trim sub(E,Gr)

-- check smoothness
(decompose ideal singularLocus E)#0 : sub(ideal(gens R),Gr)

-- compute genus
genus(Proj R/E)
