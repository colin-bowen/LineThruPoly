function [ dist_in_poly, index_bus_in_poly, index_line_in_poly] ...
    = LineThruPoly(polygonx, polygony, node_from_index, node_to_index, node_from_location, node_to_location, polycost)
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
dist_total = zeros(length(node_from_index),length(polygonx));
dist_in_poly = zeros(length(node_from_index),4); %preallocate output matrix. col 1 is the polygon in question in which the principal node
%is located in the connection. col2 is the linear distance through poly
%col3 is the distance around the poly. col4 is the number of nodes in the
%connection contained within a poly 
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

[num_connections,~] = size(node_from_location);
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
%connections_from_size = zeros(length(I),1); %preallocate for speed

%selecting connections to analyze. use linear regression with weighted
%centers of polygons. If the linear correlation is strong enough then it
%will be submitted for checking for the values of the intersection of
%polygon

%firstly, calculate centers of each polygon
polygon_center = zeros(row_polygon,2);
for i = 1:row_polygon
    polygon_center(i,:) = [mean(polygonx(i)),mean(polygony(i))];
end

%calculate correlation coefficients for 4 points. 1. principal node 2.
%terminal node 3. polygon center 4. midpoint of node shortest path

%firstly, calculate midpoint of node shortest path
connection_midpoint = zeros(num_connections,2);
for i = 1:num_connections
    connection_midpoint(i,:) = [mean(connection(i,3), connection(i,6)), mean(connection(i,2),connection(i,5))];
end





%perform correlation check to see if there is a strong correlation between
%node location and polygon center
check_matrix = zeros(num_connections,row_polygon);
for i = 1:num_connections
    for j = 1:row_polygon
        correlation = corr([connection(i,3), connection_midpoint(i,1), connection(i,6), polygon_center(j,1)], ...
            [connection(i,2), connection_midpoint(i,2), connection(i,5), polygon_center(j,2)]);
        if correlation > 0.9 || correlation < -0.9 %if the correlation is strongly positive or strongly negative (since it is a spatial distribution it doesn´t really matter)
            %check to see if nodes are farther away than principal node is from polygon center
            node_distance = sqrt((connection(i,3)-connection(i,6))^2 + (connection(i,2)-connection(i,5))^2);
            node_to_polygon = sqrt((connection(i,3)-polygon_center(j,1))^2 + (connection(i,2)-polygon_center(j,2))^2);
            if node_distance > node_to_polygon
                check_matrix(i,j) = 1; %check matrix will be those connections to check against which polygons
            end
        end
    end
end
check_matrix = sparse(check_matrix);
[M,N] = find(check_matrix); %M is the set of connections to be checked against the N polygons with which they could possibly be incident (m,n) is one such pair

blocking_polys = zeros(size(check_matrix)); %confirms that a poly n will block a connection m
x_intersect = zeros([size(check_matrix) length(polygonx)]);
y_intersect = zeros([size(check_matrix) length(polygony)]);
for i = 1:length(M)
    [x_cross, y_cross] = polyxpoly([connections(M(i),3), connections(M(i),6)],[connections(M(i),2), connections(M(i),5)], ...
        polygonx(N(i)), polygonx(N(i)));
    if ~isempty(x_intersect)
        blocking_polys(M(i),N(i)) = 1;
        x_intersect(M(i),N(i),1:length(x_cross)) = x_cross;
        y_intersect(M(i),N(i),1:length(y_cross)) = y_cross;
    end 
end

x_intersect = sparse(x_intersect);
y_intersect = sparse(y_intersect);
i = 1;
while 1 %chain of if statements to evaluate the properites of the connection with respect to location. i corresponds to the ith connection
    if i > 2 %shortcut to avoid checking if ismember because we already sorted
        while connections(i,1) == connections(i - 1,1) %while principal nodes are the same (because they are in sorted order)
            if ismember(connections(i,4),I)
                k = J(find(connections(i,4),1)); %find the corresp polygon for the terminal node
            %if the terminal node is also in a polygon
                if j == k
                    %if the nodes are in the same poly, then distance is simply
                    %between them times the polycost
                    dist_total(i) = deg2km(sqrt((connections(i,3)-connections(i,6))^2+(connections(i,2)-connections(i,5))^2))*polycost(j); %total distance of connection
                    break
                end
                
                
                if find(blocking_polys(:,1)==i) %if there exists at least one blocking poly corresponding to this connection
                    
                    %the nodes are in different polys but there is a blocking
                    %poly, then the distance is principal node ->
                    %boundary of poly1 -> boundary of poly2 -> other
                    %boundary of poly2 -> ... -> boundary of polyn ->
                    %terminal node. must decide whether it is cheaper to
                    %operate the line around the shortest boundary of the
                    %polygon or through the polygon based on polycost
                else 
                    %the nodes are in in different polys but there is no
                    %blocking poly,  then the distance is principal node ->
                    %boundary of poly1 -> boundary of poly2 -> terminal
                    %node
                    [x_cross, y_cross] = polyxpoly([connections(i,3), connections(i,6)],[connections(i,2), connections(i,5)],...
                        [polygonx(j), nan, polygonx(k)], [polygony(j), nan, polygony(k)]);
                    dist_total(i) = deg2km(sqrt((connections(i,3)-x_cross(1))^2 + ...
                        (connections(i,2)-y_cross(1))^2))*polycost(j) ...
                        + deg2km(sqrt((x_cross(1)-x_cross(2))^2 + ...
                        (y_cross(1)-y_cross(2))^2))...
                        + deg2km(sqrt((x_cross(2) - connections(i,6))^2 + ...
                        (y_cross(2) - connections(i,5))^2))*polycost(k);
                    
                end
            elseif 1
                %the terminal node is not in a poly, but there is a blocking poly
                %j = J(find(I==connections(i,1))); %index 
                %[x_intersect,y_intersect] = polyxpoly(
            else
                %the terminal node is not in a poly and there is not a blocking poly
            
                dist_total(i) = deg2km(sqrt((connections(i,3)-x_cross)^2 + (connections(i,2)-y_cross)^2))*polycost(j) + ...
                    deg2km(sqrt((x_cross-connections(i,6))^2 + (y_cross - connections(i,5))^2)); %total distance is sum of penalized dist and reg dist
            end
        end
    end
    if ismember(connections(i,1),I) 
        %if the principal node is in a polygon
        j = J(find(I==connections(i,1),1)); %corresp polygon to node´s location. only need to find first instance because the corresp poly
        if ismember(connections(i,4),I)
            %if the terminal node is also in a polygon
            if find(I==connections(i,1)) == find(I==connections(i,4)) 
                %if the nodes are in the same poly, then distance is simply
                %between them times the polycost
                dist_total(i) = deg2km(sqrt((connections(i,3)-connections(i,6))^2+(connections(i,2)-connections(i,5))^2))*polycost(j); %total distance of connection
            elseif 1
                %the nodes are in different polys but there is a blocking
                %poly
            else
                %the nodes are in in different polys but there is no
                %blocking poly
            end
        elseif 1
            %the terminal node is not in a poly, but there is a blocking poly
            %j = J(find(I==connections(i,1))); %index 
            %[x_intersect,y_intersect] = polyxpoly(
        else
            %the terminal node is not in a poly and there is not a blocking poly
            
            dist_total(i) = deg2km(sqrt((connections(i,3)-x_cross)^2 + (connections(i,2)-y_cross)^2))*polycost(j) + ...
                deg2km(sqrt((x_cross-connections(i,6))^2 + (y_cross-connections(i,5))^2)); %total distance is sum of penalized dist and reg dist
        end
           
                
    elseif 1
        %principal node not in poly, but is terminal node?
        if 1%principal node not in poly, terminal node is, is there a blocking poly?
        else %principal node not in poly, terminal node is, but no blocking poly
        end
        
    %for j =1:length(row_polygon)
    %[x_intersect,y_intersect] = polyxpoly([connections(i,3), connections(i,6)],[connections(i,2), connections(i,5)], polygonx(
    elseif 1
        %neither node is in a polygon, but are there any blocking polygons?
    else 
        %neither node is in a polygon and there are no intersecting polys
        dist_total(i) = deg2km(sqrt((connections(i,3)-connections(i,6))^2+(connections(i,2)-connections(i,5))^2)); %total distance of connection
    end
end
        
                
                
%for i = 1:length(I) %checking the num of connections for each node i in its corresp polygon
    %connection_indices(i) = find(node_from(:,1)==I(i));
 %   connections_from_size(i) = length(find(node_from(:,1)==I(i))); 
    
%end
%dist_node_in_poly = zeros(length(I),max(connections_from_size)); %preallocate matrix size based on i nodes in corresp polys, and the one with the most connections
%for i =1:length(I)%knowing which nodes are already in a polygon, find the distance of the connection inside the polygon
 %   node_de = find(node_from(:,1)==I(i)); %find index of node inside corresp polygon and all its connections
  %  node_a = node_to(node_de); %will return node indices of node i´s connections
   % node_de_lugar = [node_from(node_de(1),3),node_from(node_de(1),2)]; %go to first instance of node_from and collect those coordinates
   % node_a_lugar = [node_to(node_de,3),node_to(node_de,2)]; %collect coordinates of all corresponding connections to node i
   % for j = 1:length(node_de) %for all j connections of node i which resides in its ith polygon, need to check the intersections
    %    if isempty(setdiff(node_a(j),index_bus_in_poly(:,J(i)))) %if there are no elements in common between this node and the nodes in polygon i.e. this node is not in the same polygon
     %       [x_intersect,y_intersect] = polyxpoly([node_de_lugar(1),node_a_lugar(j,1)], ... will check for intersection point of polygon and connection
      %          [node_de_lugar(2), node_a_lugar(j,2)], polygonx(J(i,:)), polygony(J(i,:))); %checks polygon corresponding to node i 
       %     dist_node_in_poly = deg2km(sqrt((node_de_lugar(1)-x_intersect)^2+(node_de_lugar(2)-y_intersect)^2)); %node i in its corresp poly has j connections, whose dist in poly is recorded here
        %    connection_index = node_de(j); %returns index of the connection specified from node i to node j. works because node_de contains the indices within node_from, which has the same length as the total number of connections. the jth entry is the jth connection of type 'begins with node i', so it returns appropriately
        %
         %   dist_in_poly(connection_index,2) = dist_node_in_poly; %store distance penalty for final output
          %  dist_in_poly(connection_index,1) = J(i); %store polygon id for final output
          %  dist_in_poly(connection_index,4) = 1; %information re: number of nodes in the polygon in question
        %else
            %the polygon contains both nodes in the connection e.g. they
            %are both located in the same mountain range 
        %    dist_nodes_in_poly = deg2km(sqrt((node_de_lugar(1)-node_a_lugar(j,1))^2 + (node_de_lugar(2)-node_a_lugar(j,2))^2));
        %    connection_index = node_de(j);
         %   dist_in_poly(connection_index,2) = dist_nodes_in_poly;
         %   dist_in_poly(connection_index,1) = J(i);
         %   dist_in_poly(connection_index,4) = 2; %two nodes contained inside the polygon in question
       % end 
   % end
%end



%now need to check the nodes that are not contained within any polygons of
%interest and their connections. to be more computationally efficient, we
%can reject connections that have a low probability of crossing a polygon
%by comparing the proximity of the nodes to the polygon. 
%then need to evaluate the cheapest approach: through the polygon, or
%around the polygon. penalty will be assigned depending on the type of
%obstacle. 


%K = setdiff(all_nodes, I); %those nodes that are not in any polygon


%for i = 1:length(K)
 %   connections_K = find(node_from(:,1)==K(i));
  %  for k = 1:length(connections_K)
   %     for j = 1:length(row_polygon)
    %    [x_intersect,y_intersect] = polyxpoly([node_from_location(connections_K(k),2),node_to_location(connections_K(k),2)], ...
     %       [node_from_location(connections_K(k),1), node_to_location(connections_K(k),1)], polygonx(j,:),polygony(j,:)); % checks for intersections 
      %      if ~isempty(x_intersect) % if a connection exists
       %         x_cross_first = x_intersect(1);
        %        y_cross_first = y_intersect(1);
         %       x_cross_last = x_intersect(2);
          %      y_cross_last = y_intersect(2);
           %     dist_in_poly(connections_K(k),j) = deg2km(sqrt((x_cross_first - x_cross_last)^2 + (y_cross_first - y_cross_last)^2));
           % end
        %end
    %end
%end
            



%then establish node connections
%x_cross_first = zeros(length(node_from),length(polygonx));
%y_cross_first = zeros(length(node_from),length(polygonx));
%x_cross_last = zeros(length(node_from),length(polygonx));
%y_cross_last = zeros(length(node_from),length(polygonx));
%dist_in_poly = zeros(length(node_from),length(polygonx));
%index_line_in_poly = zeros(length(node_from),2);
%index_bus_in_poly = zeros(length(node_from),2);
%for i = 1:row_location
 %   for j = 1:row_polygon
  %      [x_intersect,y_intersect] = polyxpoly([node_from_location(i,2),node_to_location(i,2)], ...
   %         [node_from_location(i,1), node_to_location(i,1)], polygonx(j,:),polygony(j,:)); % checks for intersections 
    %    if ~isempty(x_intersect)
     %       x_cross_first(i,j) = x_intersect(1);
      %      y_cross_first(i,j) = y_intersect(1);
       %     if ~inpolygon(node_from_location(i,2),node_from_location(i,1),polygonx(j,:),polygony(j,:)) && ... %will check if both nodes are outside polygon
        %            ~inpolygon(node_to_location(i,2),node_to_location(i,1),polygonx(j,:),polygony(j,:))
         %       x_cross_last(i,j) = x_intersect(2);
          %      y_cross_last(i,j) = y_intersect(2);
           %     dist_in_poly(i,j) = deg2km(sqrt((x_cross_first - x_cross_last)^2 + (y_cross_first - y_cross_last)^2)); % distance that connection i spends in polygon j
            %    index_line_in_poly(i,:) = [i j]; %returns index of the connection i in polygon j 
            %%elseif index_bus_in_poly(i,j) %inpolygon(node_from_location(i,2),node_from_location(i,1),polygonx(j,:),polygony(j,:)) % if the node from is inpolygon
             %   dist_in_poly(i,j) = deg2km(sqrt((node_from_location(i,2) - x_intersect(1))^2 ...
       %             + (node_from_location(i,1) - y_intersect(1))^2)); 
        %        index_bus_in_poly(i,:) = [node_from(1,i),j];
         %       
          %  elseif ind%inpolygon(node_to_location(i,2),node_to_location(i,1),polygonx(j,:),polygony(j,:)) %if the node to is inpolygon
           %     dist_in_poly(i,j) = deg2km(sqrt((node_to_location(i,2) - x_intersect(1))^2 ...
            %        + (node_to_location(i,1) - y_intersect(1))^2));
             %   index_bus_in_poly(i,:) = [node_to(1,i),j];

            %else %both nodes are in polygon
            %    dist_in_poly(i,j) = deg2km(sqrt((node_from_location(i,2) - node_to_location(i,2))^2 + ...
                %    (node_from_location(i,1) - node_to_location(i,1))^2));
 %               index_bus_in_poly(i,:) = [1e+6,j]; %using 1e+6 to signify that both nodes are in the polygon
  %              
   %         end
    %    end
    %end
%end
%dist_in_poly = sparse(dist_in_poly);
%index_bus_in_poly = sparse(index_bus_in_poly);
%index_line_in_poly = sparse(index_line_in_poly);
end


