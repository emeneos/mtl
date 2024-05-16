mtlroot = matlabroot;
if(ispc)
    LPATH = fullfile(matlabroot,'extern','lib',computer('arch'),'microsoft');
    mex( '-R2018a', 'CXXFLAGS="$CXXFLAGS -fopenmp -D_USE_OMP_THREAD_CONTROL"','LDFLAGS="$LDFLAGS -fopenmp"', ...
    ['-L',LPATH], ...
    'dmri_2F1mex.cpp', 'hypergeom2F1.cxx', 'threadHelper.cpp' );
else
    libsdir = sprintf('%s/bin/glnxa64',mtlroot);
    mkllib  = sprintf('%s/mkl.so',libsdir);
    mex( '-R2018a', 'CXXFLAGS="$CXXFLAGS -D_USE_MKL_THREAD_CONTROL"', ...
    'dmri_2F1mex.cpp', 'hypergeom2F1.cxx', 'threadHelper.cpp', ...
    '-lpthread', mkllib );
end

fprintf(1,'Running test...');
test_dmri_2F1mex;
fprintf(1,' Done\n');
zoom('on');
fprintf(1,'Zoom in so that you can check circles are aligned with stars far from x=1\n');
