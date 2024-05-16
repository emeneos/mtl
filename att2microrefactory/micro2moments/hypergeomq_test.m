m = 5;
n = 202;
x = logspace(-154,0,n);

 

h1 = hypergeom([0.5 0.5],[1.5 1.5],x);


  

h2 = hypgeo_mex(0.5, 0.5, 1.5, x);

