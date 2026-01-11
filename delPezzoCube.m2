-- setup the projective plane
R = QQ[x,y,z];

-- collect all degree 3 monomials in 3 variable
B = basis(3,R); -- the parameter space of cubics in the plane

-- ring with params corresponding to the basis of cubics
E = QQ[{x,y,z} | toList(c_0..c_(numcols B -1)) | {t1 , t2, t}]; -- the t's parameterize the families of cubics we will find later on; t is for plotting purposes t1,t2 for algebra purpose in here
B = sub(B,E);
C = toList(c_0..c_(numcols B -1)); -- makes a list of the coefficients 


-- fix 6 points in the plane
P = {{0,1,1}, {-1,0,1}, {0,-1,1}, {1, 1/2,1}, {3/4,-1/2,1}, {-1/4, -1/2,1}}; -- these correspond to A,B,C,D,E,F in geogebra

-- the qudrilateral Q is bounded by
L = ideal(y-x-z); -- the line through A and B
S = ideal(x^2 + y^2 -1/2*x*y -z^2); -- interpolates P#0--P#4

-- generic linear combination of basis elements of cubics
g = sum( for i to numcols B -1 list B_(0,i)*c_i);


-------------------------------------------------------------
--- collect conditions for making red nodal cubic at A ---

-- RED gradient vanishes at (0,1) = P#0 
dxRed = (coefficients( sub(diff(x,g), {x=>0, y=>1, z=> 1}), Monomials => C))#1;
dyRed = (coefficients( sub(diff(y,g), {x=>0, y=>1, z=> 1}), Monomials => C))#1;


-- these are the interpolation conditions of the 6 points
-- built function for containment
mCon = P -> (
    tmpCon := for i from 0 to #P-1 list sub(g, {x => (P#i)#0 , y => (P#i)#1 , z => 1});
    matrix{apply(tmpCon, p -> (coefficients(p, Monomials => C))#1)}
);

-------------------------------------------------------------
---- make faces corresponding to the exceptional divisor ----


---------- collect all linear constraints for RED ----------
mRed = transpose sub(dxRed | dyRed | mCon(P),QQ); -- all conditions for containing the points and having a node in A

-- get the solution space
kRed = kernel mRed;

-- compute cubic family
red = (B*(t1*(gens kRed)_{0} + t2*(gens kRed)_{1} ))_(0,0);

-- for plotting choose a generic t1 value and substitute t2 with t, also choice of affine patch U_z
planeRed = toString sub(red, {t1=>1, t2=>t, z=>1});

---------- make blue nodal cubic ----------

-- this obtained as a pullback of 2-plane containing the exceptional divisor which makes the cube nicely bounded
blue = -(5/3)*x^3+(5/6)*x^2*y+(115/12)*x*y^2-(29/12)*x^2*z+(5/6)*x*y*z-(11/12)*y^2*z+(1/6)*x*z^2+(11/12)*z^3;

-- for plotting we choose the affine chart z=1
planeBlue = sub(blue, z => 1);


-------------------------------------------------------------
------- make the other 2 (pencils of) linear faces ----------

---- 1: cubic containing L and conic containing remaining 4 pts
Bconic = sub(basis(2,R),E);
Pconic  = drop(P,2);

gConic = sum( for i to numcols Bconic -1 list Bconic_(0,i)*c_i);

-- only need take up to degree 2 elements of basis
mConic = for p in Pconic list sub(gConic, {x => p#0, y=> p#1, z => 1})
MConic  = matrix{apply( mConic, e -> (coefficients(e, Monomials => {c_0,c_1,c_2,c_3,c_4,c_5}))#1)}

-- compute solution space
kConic = kernel transpose MConic;

-- compute cubic family
conic = (L_0)*((Bconic*(t1*(gens kConic)_{0} + t2*(gens kConic)_{1} ))_(0,0));

-- for plotting choose a generic t1 value and substitute t2 with t, also choice of affine patch U_z
planeConic = toString sub(conic, {t1=>1, t2=>t, z=>1});


-------------------------
---- 2: cubic containing the conic S and the line through the remaining point P#5
lineEq = y + 1/4*x + 9/16*z; 

-- compute cubic family by multiplying with line
line = lineEq*(S_0);

-- for plotting choose affine patch U_z
planeLine =  toString sub(line, z=>1);

--------------------------------------------------------
------- make the remaning opposing BLACK face ----------

-- this obtained as a pullback of 2-plane very close to our quadrilateral on the cubic surface
-- old oppo = -(144281346333049/62393779574120)*x^3-(75702941361939/124787559148240)*x^2*y+(328368662103031/49915023659296)*x*y^2-(2009304096539/1177241124040)*y^3-(710931951169261/249575118296480)*x^2*z+(709818682343/588620562020)*x*y*z-(8376116223613/4708964496160)*y^2*z+(731433004751/588620562020)*x*z^2+(2009304096539/1177241124040)*y*z^2+(8376116223613/4708964496160)*z^3;
oppo = (28557/2650)*x^3-(12727/5300)*x^2*y-(37547/10600)*x*y^2+y^3+(55451/10600)*x^2*z-(391/100)*x*y*z+(691/200)*y^2*z-9*x*z^2-y*z^2-(691/200)*z^3;

-- for plotting on the affine chart z=1
planeOppo = sub(oppo, z => 1);

-----------------------------------------------------------------------------
--------- compute the basis of cubics interpolating the 6 points ------------

-- get the conditions on the c_i that interpolating 6 points puts on them
basisCubics = transpose sub( mCon(P),QQ);

-- compute the 4 dimensional solution space of c_i
bCubics =  gens kernel basisCubics;

-- replace it with a nice plane at infinity in the last column which is our oppo face
bCubics = bCubics_{0,1,2} | matrix {{-54/53}, {-3/53}, {-161/106}, {717/106}, {1}, {0}, {0}, {-1/2}, {0}, {1/2}}; -- this equaled 1 col + 1/10 4th col

-- get explicit equations of the 4 basis cubics
gCubics = for i to 3 list sub(g, matrix {{x,y,z}} | transpose bCubics_{i} | matrix {{1,1,1}});

-- for plotting choose affine patch U_z
toString sub(gCubics#0 + t*gCubics#3, z=>1);

------------------------------------------------------------------------------
---------- express cubics in the basis cub1, cub2, cub3, cub4  ---------------

-- compute the coefficients of the chosen cubics in the standard monomial basis, that is with a choice for the params t1, t2
redt = sub((coefficients(sub(red, {t1=>1, t2=>-1}), Monomials => flatten entries B))#1,QQ);
conict = sub((coefficients(sub(conic, {t1=>1, t2=>-1/10}), Monomials => flatten entries B))#1,QQ);
linet = sub((coefficients(sub(line, {t1=>1, t2=>-1/10}), Monomials => flatten entries B))#1,QQ);

-- setup ring with coefficients corresponding to basis elements
cubicR = QQ[u_0..u_3, x,y,z];
M = matrix{ drop(gens cubicR,-3)};

-- compute normal vector of the RED hyperplane in PP^3 = CC[cub1, cub2, cub3, cub4]
redSurface = (M*solve(bCubics, redt))_(0,0);
plotRed = sub(M*solve(bCubics, redt), {u_0=>x, u_1=>y, u_2=>z, u_3=>1});

-- compute normal vector of the BLUE hyperplane in PP^3 = CC[cub1, cub2, cub3, cub4]
blueSurface = -u_0+(1/6)*u_1+(11/6)*u_3;
plotBlue = sub(blueSurface, {u_0=>x, u_1=>y, u_2=>z, u_3=>1});

-- compute normal vector of the CONIC (GREEN) hyperplane in PP^3 = CC[cub1, cub2, cub3, cub4]
conicSurface = (M*solve(bCubics, conict))_(0,0);
plotConic = sub(M*solve(bCubics, conict), {u_0=>x, u_1=>y, u_2=>z, u_3=>1});

-- compute normal vector of the LINE (YELLOW) hyperplane in PP^3 = CC[cub1, cub2, cub3, cub4]
lineSurface = (M*solve(bCubics, linet))_(0,0);
plotLine = sub(M*solve(bCubics, linet), {u_0=>x, u_1=>y, u_2=>z, u_3=>1});

-- simpler oppo that also makes one connected component and nicer picture
oppoSurface = -(691/100)*u_3 + 3*u_0 - 9*u_1 - u_2;
plotOppo = sub(oppoSurface, {u_0=>x, u_1=>y, u_2=>z, u_3=>1});

--------------------------------------------------------------------------------
------------------------- compute the blow up ----------------------------------

-- function for maximal ideals of the chosen points
maxP = p -> (
    ideal((x - p#0*z, y - p#1*z))
);

-- compute all of the blow-ups by eliminating the graph ideal for the basis cubics and the new coordinates u_i
graphI = minors(2,matrix {toList(u_0..u_3)} || sub(matrix{gCubics}, cubicR)) :  intersect(ideal(x,y,z), ideal(u_0..u_3)); -- we saturate the irrelevant ideals of projective space(s)

-- saturate the exceptional divisors away
satI = graphI : intersect( ( for p in P list maxP(p))); -- these are the exceptional divisors over the 6 points

-- get equation of the blown-up surface in the space of cubics
cubicSurface = (eliminate({x,y,z}, satI))_0;

-- for plotting we change variables to x,y,z and affine patch u_3=1
plotCubic = sub(cubicSurface, {u_0=>x, u_1=>y, u_2=>z, u_3=>1});

--------------------------------------------------------------------------
---------------------------- compute adjoint  ----------------------------

-- collect all residual intersections AFTER fixing semi-algebraic conditions
L1 = ideal(redSurface) + ideal(blueSurface);
L2 = ideal(conicSurface) + ideal(lineSurface);
C1 = (decompose(ideal(cubicSurface) + ideal(redSurface)))#1;
C2 = (decompose(ideal(cubicSurface) + ideal(blueSurface)))#1;
C3 = (decompose(ideal(cubicSurface) + ideal(conicSurface)))#1;
C4 = (decompose(ideal(cubicSurface) + ideal(lineSurface)))#1;
El = ideal(cubicSurface) + ideal(oppoSurface); -- is irreducible


-- compute space of adjoints
adjI = intersect( L1, L2, C1, C2, C3, C4, El);

-- select UNIQUE generator of degree 4
adj = ideal((select(adjI_*, f -> degree f == {4}))#0);

-- plot adjoint (use x,y,z coordinates for the PP^3 of cubics)
plotAdj = toString (sub( adj, {u_3=>1, u_0=>x, u_1=>y, u_2=>z}))_0;


--------------------------------------------------------------------------------
------- sample points in plane regions ------

for g in gCubics list sub(g, {x=>0,y=>0,z=>1});
projPt = for g in gCubics list sub(g, {x=>-6/10,y=>6/10,z=>1});
affPt = toString ((projPt_{0,1,2})*(projPt#3)^(-1));


--------------------------------------------------------------------------------
---- to get equations of the YELLOW (LINE) cubic in the plane call 
planeLine

-- to get the hyperplane corresponding to YELLOW (LINE) PP^3 with coordinates [cub1,cub2,cub3,cub4] call
plotLine

----------------------------------------------------------------------
-------------------------- compute vertices --------------------------

VOLG = decompose ideal(oppoSurface, lineSurface, conicSurface);
VOLR = decompose ideal(oppoSurface, lineSurface, redSurface);
VOLB = decompose ideal(oppoSurface, lineSurface, blueSurface);
VCLR = decompose ideal(cubicSurface, lineSurface, redSurface);
VCGR = decompose ideal(cubicSurface, conicSurface, redSurface);
VCLB = decompose ideal(cubicSurface, lineSurface, blueSurface);
VOCB = decompose ideal(oppoSurface, conicSurface, blueSurface);


----------------------------------------------------------------------
------------------------ resiude computations ------------------------

-- compute the boundary i.e. the denominator of the canonical form
boundary = redSurface * blueSurface * lineSurface * conicSurface * oppoSurface * cubicSurface;

-- compute the canonical form
canForm = adj_0/boundary;


----------------------------------------------------------
------ residue of the red line on the cubic surface ------

use cubicR;

-- compute differential
d0Cubic = diff(u_0, cubicSurface);
d1Cubic = diff(u_1, cubicSurface);
d1Red = diff(u_1, redSurface);
d0Red = diff(u_0, redSurface);

-- red line resiude (before substituting params)
redRes = canForm * cubicSurface * redSurface/(d0Cubic*d1Red-d0Red*d1Cubic);

-- setup subsitution (i.e. solutions for the params)
g = (decompose ideal(cubicSurface,redSurface))#0;
u1Sub = sub((-g_0 + diff(u_1,g_0)*u_1)/diff(u_1,g_0), cubicR);
u0Sub = sub((-g_1 + diff(u_0,g_1)*u_0)/diff(u_0,g_1), cubicR);
redRes = sub(sub(redRes, {u_1=>u1Sub}), {u_0=> u0Sub});

-- show boundaries
factor redRes

-- compute residue along vertex of cubic--red--yellow/line (vcry)
vcry = ((factor denominator redRes)#0)#0 -- get vertex cooridantes
vcrySub = sub((-vcry + diff(u_2,vcry)*u_2)/diff(u_2,vcry), cubicR); -- substitution rule for local param
d2Yellow = diff(u_2,vcry); -- differential wrt. to local param u_2

-- residue at vcry
resCRY = sub(sub(redRes * sub(vcry, frac(cubicR)) * 1/d2Yellow, {u_2=>vcrySub}), {u_3=>1})


---------- same vertex but in order cubic--yellow/line--red

-- compute differential
d0Cubic = diff(u_0, cubicSurface);
d1Cubic = diff(u_1, cubicSurface);
d1Yellow = diff(u_1, lineSurface);
d0Yellow = diff(u_0, lineSurface);

-- yellow/line line resiude (before substituting params)
yellowRes = canForm * cubicSurface * lineSurface/(d0Cubic*d1Yellow-d0Yellow*d1Cubic);

-- setup subsitution (i.e. solutions for the params)
g = (decompose ideal(cubicSurface,lineSurface))#0;
u1Sub = sub((-g_0 + diff(u_1,g_0)*u_1)/diff(u_1,g_0), cubicR);
u0Sub = sub((-g_1 + diff(u_0,g_1)*u_0)/diff(u_0,g_1), cubicR);
yellowRes = sub(sub(yellowRes, {u_1=>u1Sub}), {u_0=> u0Sub});

-- show boundaries
factor yellowRes

-- compute residue along vertex of cubic--yellow--red (vcyr)
vcyr = ((factor denominator yellowRes)#0)#0; -- get vertex cooridantes
vcyrSub = sub((-vcyr + diff(u_2,vcyr)*u_2)/diff(u_2,vcyr), cubicR); -- substitution rule for local param
d2Red = diff(u_2,vcyr); -- differential wrt. to local param u_2

-- residue at vcyr
resCYR = sub(sub(yellowRes * sub(vcyr, frac(cubicR)) * 1/d2Red, {u_2=>vcyrSub}), {u_3=>1})



-------------------------------------------------------
------ residue of the red line on the oppo facet ------

use cubicR;

-- compute differential
d0Oppo = diff(u_0, oppoSurface);
d1Oppo = diff(u_1, oppoSurface);
d1Red = diff(u_1, redSurface);
d0Red = diff(u_0, redSurface);

-- red line resiude on oppo facet (before substituting params)
redRes = canForm * oppoSurface * redSurface/(d0Oppo*d1Red-d0Red*d1Oppo);

-- setup subsitution (i.e. solutions for the params)
g = (decompose ideal(oppoSurface,redSurface))#0;
u1Sub = sub((-g_0 + diff(u_1,g_0)*u_1)/diff(u_1,g_0), cubicR);
u0Sub = sub((-g_1 + diff(u_0,g_1)*u_0)/diff(u_0,g_1), cubicR);
redRes = sub(sub(redRes, {u_1=>u1Sub}), {u_0=> u0Sub});

-- show boundaries
factor redRes

-- compute residue along vertex of oppo--red--yellow/line (vory)
vory = ((factor denominator redRes)#0)#0 -- get vertex cooridantes
vorySub = sub((-vory + diff(u_2,vory)*u_2)/diff(u_2,vory), cubicR); -- substitution rule for local param
d2Yellow = diff(u_2,vory); -- differential wrt. to local param u_2

-- residue at vcry
resORY = sub(sub(redRes * sub(vory, frac(cubicR)) * 1/d2Yellow, {u_2=>vorySub}), {u_3=>1})
