clear; 
clc;
close all
%Make all radial cases compatible with power models
files = dir()
for j=1:length(files)
    f = files(j)
    if ~f.isdir
        name = f.name
        mpc = loadcase(name)
        savecase(name,mpc)
    end
end