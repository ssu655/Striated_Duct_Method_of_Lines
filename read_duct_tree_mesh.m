%
% read in the duct tree file
%
% J.Rugis
% 27.09.21
%
%

function [nodes, segments, radii, seg_types] = read_duct_tree_mesh(fname)
  % get the data counts
  pfile = fopen(fname,'r');
  nnodes = get_count(pfile, 'duct_node');           
  nsegs = get_count(pfile, 'duct_segment');
  skip_header(pfile);

  % get lumen node data
  nodes = zeros(nnodes,3);
  radii = zeros(nnodes,1);
  for i = 1:nnodes
    tokens = str2double(split(fgetl(pfile)));
    nodes(i,:) = tokens(1:3);
    radii(i) = tokens(4);
  end
  
  % get lumen segment data
  segments = zeros(nsegs,2);
  seg_types = zeros(nsegs,1);
  for i = 1:nsegs
    tokens = str2double(split(fgetl(pfile)));
    segments(i,:) = tokens(1:2);
    seg_types(i,:) = tokens(3);
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

