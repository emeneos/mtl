function options = create_atti2micro_options(mask,M,N,P,G)
options.lambda = 0.006;    
options.tl = 1.0e-6;
options.tu = 1-1.0e-6;
options.ADC0 = 3.0e-3;
options.usef = false;
options.mask = mask;
options.bth = 100;
options.mlperp = 0.0;
options.mu = 0.0;
options.verbose = true;
options.L = 6;

options.tl = 1.0e-7;        
options.tu = 1-options.tl;      
options.ADC0 = 3.0e-3;      
options.usef = false;      
% -------------------------------------------------------------------------
options.mlperp = 0.01e-3;  
options.mu = 5.0e-5;        
options.nu = 0.0;           
options.nmax = 100;         
options.dl = 1.0e-8;        
options.dC = 1.0e-6;        
options.regf = true;        
options.lpar = ones(M,N,P);    
options.forcelpar = false;  
options.nolpwarn = false;   
% -------------------------------------------------------------------------
options.fmin = 0.01;        
options.fmax = 1;          
% -------------------------------------------------------------------------
options.mask = true(M,N,P); 
options.bth = 1;            
options.chunksz = 10000;    
options.verbose = false;  
end