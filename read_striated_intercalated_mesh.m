%
% read in a striated cell mesh file
%
% J.Rugis
% 27.09.21
%
%

function [verts, faces, tri_type, duct_idx, dist_from_duct, dist_along_duct] = read_striated_intercalated_mesh(fname)
  % get the data counts
  pfile = fopen(fname,'r');
  nverts = get_count(pfile, 'vertex');           
  nfaces = get_count(pfile, 'face');
  skip_header(pfile);

  % get the vertex data
  verts = zeros(nverts,3);
  for i = 1:nverts
    verts(i,:) = str2double(split(fgetl(pfile)));
  end
  
  % get the face data
  faces = zeros(nfaces,3);
  tri_type = zeros(nfaces,1);
  duct_idx = zeros(nfaces,1);
  dist_from_duct = zeros(nfaces,1);
  dist_along_duct = zeros(nfaces,1);
  for i = 1:nfaces
    tokens = str2double(split(fgetl(pfile)));
    faces(i,:) = tokens(2:4)+1;
    tri_type(i) = tokens(5);
    duct_idx(i) = tokens(6);
    dist_from_duct(i) = tokens(7);
    dist_along_duct(i) = tokens(8);
  end

  fclose(pfile);
end

function [count] = get_count(pfile, value)
  while 1
    tokens = split(fgetl(pfile));
    if strcmp(tokens{1},'element') && strcmp(tokens{2},value)
      count = str2double(tokens{3});
      break;
    end
  end
end

function [] = skip_header(pfile)
  while 1
    if(strcmp(fgetl(pfile),'end_header'))
      break;
    end
  end
end

