% the length of each disc by construction
disc_length = s_lumen_prop.disc_length;

% the 3D coordinates of each disc node, node 1 is root node
coordinates = s_lumen_prop.disc_nodes_coord;

X = coordinates(:,1);
Y = coordinates(:,2);
Z = coordinates(:,3);
scatter3(X,Y,Z)

% the [input node, output node] vector of each disc, disc 1 is root disc
nodes_ind = s_lumen_prop.disc_2_nodes;

% n is number of discs in the duct
n = size(nodes_ind,1);

fprintf("disc no.| in node no.| out node no.| recorded length | measured length \n",i,in_node,out_node,norm(distance),disc_length(i))

for i = 1:n
    % input node index of disc i
    in_node = nodes_ind(i,1);
    in = [X(in_node),Y(in_node),Z(in_node)];

    % output node index of disc i
    out_node = nodes_ind(i,2);
    out = [X(out_node),Y(out_node),Z(out_node)];

    % check the distance between input and output nodes
    distance = in - out;
    fprintf("%7g |  %9g |  %10g |  %14.2f |  %14.2f \n",i,in_node,out_node,disc_length(i),norm(distance))
end
