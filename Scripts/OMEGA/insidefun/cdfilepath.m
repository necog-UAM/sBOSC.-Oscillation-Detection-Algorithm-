function outfile = cdfilepath(fpath,filename)

cd(fpath)
a = dir;
for i = 1:length(a)
    if strfind(a(i).name,filename) > 0
        outfile = a(i).name;
    end
end
cd(outfile)