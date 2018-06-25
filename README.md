# LineThruPoly
will check if node connections pass through certain obstacles and assign associated penalties to those lines


function [ dist_in_poly, index_bus_in_poly, index_line_in_poly] ...
    = LineThruPoly(polygonx, polygony, node_from_index, node_to_index, node_from_location, node_to_location)
%function will accept vector of transmission lines from, and to, determine
%their Euclidean distance in the plane, and determine whether their
%Euclidean distance passes through certain input polygons, such as lakes
%and mountains. furthermore, the function will determine whether the bus is
%contained within such an input polygon. The function will determine if it
%is cheaper to go 'around' the polygon or through it. 
%   outputs are the vectors of the coordinates where they first enter the
%   polygon, then when they leave the polygon, the distance given by that,
%   and the indices of buses in the polygon as well as indices of
%   connections in the polygon. INPUT LATITUDE THEN LONGITUDE FOR POLYGON
%   AND NODE LOCATION. polygonx is longitude, polygony is latitude. 
node_from = [node_from_index, node_from_location];
node_to = [node_to_index, node_to_location];
new_node_from = zeros(length(node_from_index),3);
for i = 1:length(node_from) %switch all node_froms to be the minimum node index in the pair
    if node_from(i,1) > node_to(i,1)
        new_node_from(i,:) = node_to(i,:);
        node_to(i,:) = node_from(i,:);
        node_from(i,:) = new_node_from(i,:);
    end
end
connections = [node_from, node_to];
connections = sortrows(connections,1); %will sort rows based on node_froms        
all_nodes = unique([node_from; node_to],'rows'); %will return the unique rows and reject the ones that are identical i.e. same index and location
all_nodes_x = all_nodes(:,3); %collects longitudes for y data
all_nodes_y = all_nodes(:,2); %collects latitudes for x data

index_bus_in_poly = zeros(length(all_nodes),length(polygonx));

[row_location,~] = size(node_from_location);
[row_nodes,~] = size(all_nodes);
[row_polygon,~] = size(polygonx);
%firstly, can find nodes in polygon
for i = 1:row_nodes
    for j = 1:row_polygon
        if inpolygon(all_nodes_x(i),all_nodes_y(i),polygonx(j,:), polygony(j,:))
            index_bus_in_poly(i,j) = 1;
        end
    end
end


index_bus_in_poly = sparse(index_bus_in_poly); %returns the node index i in polygon j


[I,J] = find(index_bus_in_poly); %finding the I nodes and their corresponding J polygons. ith row corresponds to the ith node, jth column corresp to jth polygon
connections_from_size = zeros(length(I),1); %preallocate for speed

%for i = 1:length(I) %checking the num of connections for each node i in its corresp polygon
    %connection_indices(i) = find(node_from(:,1)==I(i));
 %   connections_from_size(i) = length(find(node_from(:,1)==I(i))); 
    
%end
%dist_node_in_poly = zeros(length(I),max(connections_from_size)); %preallocate matrix size based on i nodes in corresp polys, and the one with the most connections
for i =1:length(I)%knowing which nodes are already in a polygon, find the distance of the connection inside the polygon
    node_de = find(node_from(:,1)==I(i)); %find index of node inside corresp polygon and all its connections
    node_a = node_to(node_de); %will return node indices of node i´s connections
    node_de_lugar = [node_from(node_de(1),3),node_from(node_de(1),2)]; %go to first instance of node_from and collect those coordinates
    node_a_lugar = [node_to(node_de,3),node_to(node_de,2)]; %collect coordinates of all corresponding connections to node i
    for j = 1:length(node_de) %for all j connections of node i which resides in its ith polygon, need to check the intersections
            if isempty(setdiff(node_a(j),index_bus_in_poly(:,J(i)))) %if there are no elements in common between this node and the nodes in polygon i.e. this node is not in the same polygon
            [x_intersect,y_intersect] = polyxpoly([node_de_lugar(1),node_a_lugar(j,1)], ... will check for intersection point of polygon and connection
                [node_de_lugar(2), node_a_lugar(j,2)], polygonx(J(i,:)), polygony(J(i,:))); %checks polygon corresponding to node i 
            dist_node_in_poly = deg2km(sqrt((node_de_lugar(1)-x_intersect)^2+(node_de_lugar(2)-y_intersect)^2)); %node i in its corresp poly has j connections, whose dist in poly is recorded here
            connection_index = node_de(j); %returns index of the connection specified from node i to node j. works because node_de contains the indices within node_from, which has the same length as the total number of connections. the jth entry is the jth connection of type 'begins with node i', so it returns appropriately
        
            dist_in_poly(connection_index,2) = dist_node_in_poly; %store distance penalty for final output
            dist_in_poly(connection_index,1) = J(i); %store polygon id for final output
            dist_in_poly(connection_index,4) = 1; %information re: number of nodes in the polygon in question
        else
            %the polygon contains both nodes in the connection e.g. they
            %are both located in the same mountain range 
            dist_nodes_in_poly = deg2km(sqrt((node_de_lugar(1)-node_a_lugar(j,1))^2 + (node_de_lugar(2)-node_a_lugar(j,2))^2));
            connection_index = node_de(j);
            dist_in_poly(connection_index,2) = dist_nodes_in_poly;
            dist_in_poly(connection_index,1) = J(i);
            dist_in_poly(connection_index,4) = 2; %two nodes contained inside the polygon in question
        end 
    end
end



%now need to check the nodes that are not contained within any polygons of
%interest and their connections. to be more computationally efficient, we
%can reject connections that have a low probability of crossing a polygon
%by comparing the proximity of the nodes to the polygon. 
%then need to evaluate the cheapest approach: through the polygon, or
%around the polygon. penalty will be assigned depending on the type of
%obstacle. 


K = setdiff(all_nodes, I); %those nodes that are not in any polygon











if 
dist_in_poly

%then establish node connections
x_cross_first = zeros(length(node_from),length(polygonx));
y_cross_first = zeros(length(node_from),length(polygonx));
x_cross_last = zeros(length(node_from),length(polygonx));
y_cross_last = zeros(length(node_from),length(polygonx));
dist_in_poly = zeros(length(node_from),length(polygonx));
index_line_in_poly = zeros(length(node_from),2);
index_bus_in_poly = zeros(length(node_from),2);
for i = 1:row_location
    for j = 1:row_polygon
        [x_intersect,y_intersect] = polyxpoly([node_from_location(i,2),node_to_location(i,2)], ...
            [node_from_location(i,1), node_to_location(i,1)], polygonx(j,:),polygony(j,:)); % checks for intersections 
        if ~isempty(x_intersect)
            x_cross_first(i,j) = x_intersect(1);
            y_cross_first(i,j) = y_intersect(1);
            if ~inpolygon(node_from_location(i,2),node_from_location(i,1),polygonx(j,:),polygony(j,:)) && ... %will check if both nodes are outside polygon
                    ~inpolygon(node_to_location(i,2),node_to_location(i,1),polygonx(j,:),polygony(j,:))
                x_cross_last(i,j) = x_intersect(2);
                y_cross_last(i,j) = y_intersect(2);
                dist_in_poly(i,j) = deg2km(sqrt((x_cross_first - x_cross_last)^2 + (y_cross_first - y_cross_last)^2)); % distance that connection i spends in polygon j
                index_line_in_poly(i,:) = [i j]; %returns index of the connection i in polygon j 
            elseif index_bus_in_poly(i,j) %inpolygon(node_from_location(i,2),node_from_location(i,1),polygonx(j,:),polygony(j,:)) % if the node from is inpolygon
                dist_in_poly(i,j) = deg2km(sqrt((node_from_location(i,2) - x_intersect(1))^2 ...
                    + (node_from_location(i,1) - y_intersect(1))^2)); 
                index_bus_in_poly(i,:) = [node_from(1,i),j];
                
            elseif ind%inpolygon(node_to_location(i,2),node_to_location(i,1),polygonx(j,:),polygony(j,:)) %if the node to is inpolygon
                dist_in_poly(i,j) = deg2km(sqrt((node_to_location(i,2) - x_intersect(1))^2 ...
                    + (node_to_location(i,1) - y_intersect(1))^2));
                index_bus_in_poly(i,:) = [node_to(1,i),j];

            else %both nodes are in polygon
                dist_in_poly(i,j) = deg2km(sqrt((node_from_location(i,2) - node_to_location(i,2))^2 + ...
                    (node_from_location(i,1) - node_to_location(i,1))^2));
                index_bus_in_poly(i,:) = [1e+6,j]; %using 1e+6 to signify that both nodes are in the polygon
                
            end
        end
    end
end
dist_in_poly = sparse(dist_in_poly);
index_bus_in_poly = sparse(index_bus_in_poly);
index_line_in_poly = sparse(index_line_in_poly);
end
