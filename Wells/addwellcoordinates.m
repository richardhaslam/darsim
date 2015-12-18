function [Well] = addwellcoordinates(Istart, Iend, Jstart, Jend)
Iindexes = Istart:1:Iend;
Jindexes = Jstart:1:Jend;
perforations = max(length(Iindexes), length(Jindexes));
if length(Iindexes) == 1
    Well.x = Iindexes * ones(perforations,1);
    Well.y = Jindexes';
elseif length(Jindexes) == 1
    Well.x = Iindexes';
    Well.y = Jindexes * ones(perforations,1);
else
    disp('ERROR: Check the coordinates of your wells. Wells must be either horizontal or vertical!!');
    return
end
end