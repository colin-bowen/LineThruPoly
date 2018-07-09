function [ dist_total] ...
    = LineThruPoly(polygonx, polygony, node_from_index, node_to_index, node_from_location, node_to_location, polycost,resolution, ...
    goodorbad,residual, proximity)
tic
%function will accept vector of transmission lines from, and to, determine
%their Euclidean distance in the plane using weight functions to describe
%passing through "high cost" areas such as lakes or mountains and penalize
%accordingly The function will determine if it
%is cheaper to go 'around' the polygon or through it using a shortest path
%algorithm when it is determined to be necessary via linear correlation.
%the function inputs are node indices corresponding to connection, node
%location in LONGITUDE then LATITUDE (CORRESPONDING TO X AND Y), polycost which is a linear coefficient
%describing the 'cost' of going through a polygon, and resolution creates a
%resolution-by-resolution mesh when the shortest path algorithm is called

%   output is a vector of total distance of each transmission line determined by  
%   euclidean distance and weight of the blocking polygons
%   INPUT LATITUDE THEN LONGITUDE FOR 
%   NODE LOCATION. polygonx is longitude, polygony is latitude. (yes i know
%   that is confusing, that's how the NUTS data was given and how i
%   received location of node data).
dist_total = zeros(length(node_from_index),1); %dist_total returns the total distance that a connection should have to travel, minimized. 
%dist_in_poly = zeros(length(node_from_index),4); %preallocate output matrix. col 1 is the polygon in question in which the principal node
%is located in the connection. col2 is the linear distance through poly
%col3 is the distance around the poly. col4 is the number of nodes in the
%connection contained within a poly 
node_from = [node_from_index, node_from_location];
node_to = [node_to_index, node_to_location];
new_node_from = zeros(length(node_from_index),3);
for i = 1:length(node_from_index) %switch all node_froms to be the minimum node index in the pair
    if node_from(i,1) > node_to(i,1)
        new_node_from(i,:) = node_to(i,:);
        node_to(i,:) = node_from(i,:);
        node_from(i,:) = new_node_from(i,:);
    end
end
connections = [node_from, node_to];
%connections = sortrows(connections,1); %will sort rows based on node_froms        
all_nodes = unique([node_from; node_to],'rows'); %will return the unique rows and reject the ones that are identical i.e. same index and location
all_nodes_x = all_nodes(:,2); %collects longitudes for y data
all_nodes_y = all_nodes(:,3); %collects latitudes for x data
[row_polygon, boundaryx] = size(polygonx);
index_bus_in_poly = zeros(length(all_nodes),row_polygon);

[num_connections,~] = size(node_from_location);
[row_nodes,~] = size(all_nodes);

[~, boundaryy] = size(polygony);



%firstly, calculate centers of each polygon
polygon_center = zeros(row_polygon,2);
for i = 1:row_polygon
    endex = find(isnan(polygonx(i,:)),1)-1; %go to the one before the first nan
    if isempty(endex) %if there are no nans
        polygon_center(i,:) = [mean(polygonx(i,:)),mean(polygony(i,:))];
    else %if there dooby nans
        polygon_center(i,:) = [mean(polygonx(i,1:endex)),mean(polygony(i,1:endex))];
    end
end

%before finding nodes in polygon, need to flip if they're good polygons 
polygonx = [polygonx nan(row_polygon,6)];
polygony = [polygony nan(row_polygon,6)];


for i = 1:row_polygon
    if goodorbad(i) == true 
        max_x = max(polygonx(i,:));
        min_x = min(polygonx(i,:));
        max_y = max(polygony(i,:));
        min_y = min(polygony(i,:));
        new_shape_x = [max_x+1, max_x+1, min_x - 1, min_x-1, max_x+1]; %adds one degree boundary outside the shape
        new_shape_y = [min_y-1, max_y+1, max_y + 1, min_y-1, min_y-1]; %adds one degree boundary outside the shape
        %boop = 1;
        %new_shape_x = zeros(1,boundaryx); %prealloc8
        %new_shape_y = zeros(1,boundaryy); %prealloc9
        
        %while ~isnan(polygonx(i,boop))
         %   checkaroox  = (polygonx(i,boop) - polygon_center(i,1)); %creates the vector from the center of the polygon for each point in the polygon shape definition
          %  checkarooy = (polygony(i,boop) - polygon_center(i,2));
            %checkaroox = checkaroox/(sqrt(checkaroox^2+checkarooy^2));
            %checkarooy = checkarooy/(sqrt(checkaroox^2+checkarooy^2));
           % new_shape_x(boop) = polygon_center(i,1)+(checkaroox * 10);
           % new_shape_y(boop) = polygon_center(i,2)+(checkarooy * 10);
           % boop = boop + 1;
       % end
        
        polygonx(i,:) = [fliplr(polygonx(i,1:boundaryx)) nan new_shape_x]; %create a bad polygon that instead sets "outside" the good polygon as the polygon to be checked.
        polygony(i,:) = [fliplr(polygony(i,1:boundaryy)) nan new_shape_y];
        
    end
end
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



%calculate correlation coefficients for 4 points. 1. principal node 2.
%terminal node 3. polygon center 4. midpoint of node shortest path

%firstly, calculate midpoint of node shortest path
connection_midpoint = zeros(num_connections,2);
for i = 1:num_connections
    connection_midpoint(i,:) = [mean([connections(i,2), connections(i,5)]), mean([connections(i,3),connections(i,6)])];
end





%perform correlation check to see if there is a strong correlation between
%node location and polygon center
check_matrix = zeros(num_connections,row_polygon);
for i = 1:num_connections
    for j = 1:row_polygon
        if goodorbad(j) == 0 %if it's a bad polygon we can check statistical correlation
            [~,S] = polyfit([connections(i,2); connection_midpoint(i,1); connections(i,5); polygon_center(j,1)], ...
                [connections(i,3); connection_midpoint(i,2); connections(i,6); polygon_center(j,2)],1);
            if S.normr < residual %if the correlation is strongly positive or strongly negative (since it is a spatial distribution it doesn´t really matter)
                %check to see if nodes are farther away than principal node is
                %from polygon center i.e. the polygon could be in between them 
                node_distance = sqrt((connections(i,2)-connections(i,5))^2 + (connections(i,3)-connections(i,6))^2);
                node_to_polygon = sqrt((connections(i,2)-polygon_center(j,1))^2 + (connections(i,3)-polygon_center(j,2))^2);
                if node_distance > node_to_polygon || ismember(connections(i,4),I) % double conditioned OR such that if a node_to is the maximum node value, will check if it is in a polygon
                    check_matrix(i,j) = 1; %check matrix will be those connections to check against which polygons
                end
            end
        else %it's a good polygon and we'll check the inner polygon boundary to see if it's sufficiently close to one of the nodes
            counter = 1;
            minimum = 100000;
            while 1
                if isnan(counter) %only want to check the inner boundary of the polygon. If we reach this point, that means the specified minimum proximity was never reached and we can move on without checking this connection against the polygon
                    break
                end
                dist_from_poly = deg2km(sqrt((connections(i,2)-polygonx(counter))^2+(connections(i,3)-polygony(counter))^2)); %over all the shape of the inner polygon
                if dist_from_poly < minimum %if the new distance from the polygon is smaller than the previous distance, set that as the new minimum distance
                    minimum = dist_from_poly;
                end
                if minimum < proximity %if the new minimum value is less than the prespecified proximity, then we add it to the check matrix and move on
                    check_matrix(i,j) = 1;
                    break
                end
                counter = counter+1;
            end
        end
    end
end
check_matrix = sparse(check_matrix);
[M,N] = find(check_matrix); %M is the set of connections to be checked against the N polygons with which they could possibly be incident (m,n) is one such pair

blocking_polys = zeros(size(check_matrix)); %confirms that a poly n will block a connection m
x_intersect = zeros([size(check_matrix) length(polygonx)]); %x_intersect will store in a 3x3 matrix connection m's crossings with polygon n
y_intersect = zeros([size(check_matrix) length(polygony)]); %actually i
%feel like i don't need this because you're just gonna evaluate the
%shortest path anyway, although maybe it would save memory
for i = 1:length(M)
    [x_cross, y_cross] = polyxpoly([connections(M(i),2), connections(M(i),5)],[connections(M(i),3), connections(M(i),6)], ...
        polygonx(N(i),:), polygony(N(i),:));
    if ~isempty(x_cross) %if there are crossings, then we indicate that in the blocking_polys matrix 
        blocking_polys(M(i),N(i)) = 1; %confirms that checked polygon M(i) is intersected by a corresponding polygon N(i)
        x_intersect(M(i),N(i),1:length(x_cross)) = x_cross;
        y_intersect(M(i),N(i),1:length(y_cross)) = y_cross;
    end 
end


blocking_polys = sparse(blocking_polys);




[P,Q] = find(blocking_polys); %stores the indices of a connection p blocked by a poly q 
i = 1;
while 1 %chain of if statements to evaluate the properites of the connection with respect to location. i corresponds to the ith connection
    if i > num_connections
        break
    end
    [connection_blocked, poly_blocker] = find(P==i); %finds whether connection i is an element of the set of blocked connections P, and also returning the id of the blocking polygon
    if ismember(connections(i,1),I) %if the first node is contained within a polygon
        j = J(find(I==connections(i,1),1)); %corresp polygon to node´s location. only need to find first instance because the corresp poly will always be the same for the same principal node
            
            if ismember(connections(i,4),I) %if the terminal node is in a polygon
                k = J(find(I==connections(i,4),1)); %find the corresp polygon for the terminal node
            
                if j == k
                    %if the nodes are in the same poly, then distance is simply
                    %between them times the polycost
                    dist_total(i) = deg2km(sqrt((connections(i,2)-connections(i,5))^2+(connections(i,3)-connections(i,6))^2))*polycost(j); %total distance of connection
                    
                
                
                
                elseif length(connection_blocked)>2 %if there exists at least one blocking poly (not including the poly containing princ. and term. nodes) corresponding to this connection
                    
                    %run shortest path algorithm to determine the total
                    %cost. need to propose the line and blocking polygons
                    %but don't forget to add the distance inside the
                    %principal and terminal nodes
                    
                    %x_crossings = x_intersect(P(connection_blocked),Q(poly_blocker),:); %going to the pth connection in the blocking poly matrix and all of its blocking polygons and collecting the
                    %values of the intersections.
                    %y_crossings = y_intersect(P(connection_blocked),Q(poly_blocker),:);
                    poly_in = find(Q==j); %finding the polygon index in set Q of blocking polys equal to the poly in which the principal node resides NOTE THIS IS PROBABLY NOT FUNCTIONAL RN 
                    poly_out = find(Q==k); %finding the polygon index in set K of blocking polys equal to the poly in which the terminal node resides
                    count = 1;
                    f = 1; % the number of blocking polygons after excluding the polygons in which the principal and terminal nodes reside
                    blockersx = zeros(length(connection_blocked)-2, length(polygonx));
                    blockersy = zeros(length(connection_blocked)-2, length(polygony));
                    polycost_Dijkstra = zeros(length(connection_blocked)-2, 1);
                    while count <= length(connection_blocked)
                        if ismember(connection_blocked(count),poly_in)
                            dist_total(i) = dist_total(i) + deg2km(sqrt((connections(i,2)-x_intersect(i,count, 1))^2 + ...
                                (connections(i,3)-y_intersect(i,count,1))^2))*polycost(j);
                            principal = [x_intersect(i,count,1) y_intersect(i,count,1)]; %store the outbound location as the boundary of the resident polygon
                        elseif ismember(connection_blocked(count),poly_out)
                            dist_total(i) = dist_total(i) + deg2km(sqrt((connections(i,5)-x_intersect(i,count,1))^2 + ...
                                (connections(i,6)-y_intersect(i,count,1))^2))*polycost(k);
                            terminal = [x_intersect(i,count,1) y_intersect(i,count,1)]; %store the final destination as the boundary of the resident node to polygon
                        else
                            blockersx(f,:) = polygonx(Q(connection_blocked(count)),:); %store for the shortest path algo the polygons that you'll be using, associated with the 
                            blockersy(f,:) = polygony(Q(connection_blocked(count)),:);
                            polycost_Dijkstra(f) = polycost(Q(connection_blocked(count)));
                            f = f+1;
                        end
                        count = count + 1;
                    end
                    dist_total(i) = dist_total(i) + Dijkstra(principal, terminal, blockersx, blockersy, polycost_Dijkstra,resolution);
                    %first need to access the intersect vector prestored
                    %from above
                    %x_crossings = x_intersect(P(blocked),Q(blocked),:); %going to the pth connection in the blocking poly matrix and all of its blocking polygons and collecting the
                    %values of the intersections.
                    %y_crossings = y_intersect(P(blocked),Q(blocked),:);
                    %poly_in = find(Q==j); %finding the polygon index in set Q of blocking polys equal to the poly in which the principal node resides
                    %poly_out = find(Q==k); %finding the polygon index in set K of blocking polys equal to the poly in which the terminal node resides
                    
                    %for count = 1:length(blocked) %for loop to sum together all of the polygon distances for this connection
                     %   if  count==poly_in %the principal node's polygon index in Q
                      %      dist_total(i) = dist_total(i) + deg2km(sqrt((connections(i,3)-x_crossing(i,count,1))^2 + ...
                       %         (connections(i,2)-y_crossing(i,count,1))^2))*polycost(j);
                       % elseif count==poly_out %the terminal node's polygon index in Q
                        %    dist_total(i) = dist_total(i) + deg2km(sqrt((connections(i,6)-x_crossing(i,count,1))^2 + ...
                         %       (connections(i,5)-y_crossing(i,count,1))^2))*polycost(k);
                        %else %the other blocking polygons that do not contain the nodes. this will add the distances inside the polygons.
                            %here try to MINIMIZE DISTANCE by either going
                            %around or through the polygon
                         %   dist_total(i) = dist_total(i) + deg2km(sqrt((x_crossings(count,count,1)-x_crossings(count,count,2))^2 + ...
                          %      (y_crossings(count,count,1) - y_crossings(count,count,2))^2))*polycost(Q(count)); %in the countth iteration of the set of blocking polys corresp to this connection
                        %end
                        %now need to sum up interpolygon distances.
                        %if count == 1 %if the first iteration, then go from principal polygon boundary to boundary of next polygon
                         %   dist_total(i) = dist_total(i) + deg2km(sqrt((x_crossings(count,count,1) - x_crossings(count+1,count+1,1))^2 + ...
                          %      (y_crossings(count,count,1) - y_crossings(count+1,count+1,1))^2));
                            
                        %elseif count ~= length(blocked) %if not a principal or terminal polygon, then compare distances between end of the connection's crossing of polygon count and the beginning of connection's crossing of count+1
                            %dist_total(i) = dist_total(i) + deg2km(sqrt((x_crossings(count,count,2) - x_crossings(count+1,count+1,1))^2 + ...
                             %   (y_crossings(count,count,2) - y_crossings(count+1,count+1,1))^2));
                        %else
                            %hi
                        %end
                        
                           % cross = 2;
                            %while cross < length(x_crossings)
                             %   dist_total(i) = dist_total(i) + deg2km(sqrt((x_crossings(i,count,cross)+x_crossings(i,count,cross+1))^2 + ...
                              %      (y_crossings(i,cross,count)+y_crossings(i,cross,count+1))^2))*polycost(blocked(count)); 
                               % cross = cross + 1;
                                %dist_total(i) = dist_total(i) + deg2km(sqrt((x_crossings(i,count,cross)+x_crossings(i,count,cross+1))^2 + ... in between crossings
                                 %   (y_crossings(i)+y_crossing(i+1))^2));
                                %cross = cross + 1;
                            %end
                    %end
                    
                    
                    
                    %over all polygons, sum the blocked polys)
                    
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
                    [x_cross, y_cross] = polyxpoly([connections(i,2), connections(i,5)],[connections(i,3), connections(i,6)],...
                        [polygonx(j,:), nan, polygonx(k,:)], [polygony(j,:), nan, polygony(k,:)]);
                    dist_total(i) = deg2km(sqrt((connections(i,2)-x_cross(1))^2 + ...
                        (connections(i,3)-y_cross(1))^2))*polycost(j) ...
                        + deg2km(sqrt((x_cross(1)-x_cross(2))^2 + ...
                        (y_cross(1)-y_cross(2))^2))...
                        + deg2km(sqrt((x_cross(2) - connections(i,5))^2 + ...
                        (y_cross(2) - connections(i,6))^2))*polycost(k);
                    
                    
                end
            elseif length(connection_blocked)>1
                %the terminal node is not in a poly, but there is a blocking poly
                poly_in = find(Q==j); %finding the polygon index in set Q of blocking polys equal to the poly in which the principal node resides
                count = 1;
                f = 1;
                blockersx = zeros(length(connection_blocked)-1, length(polygonx));
                blockersy = zeros(length(connection_blocked)-1, length(polygony));
                polycost_Dijkstra = zeros(length(connection_blocked)-1, 1);
                while count <= length(connection_blocked) %for loop to sum together all of the polygon distances for this connection
                        if  ismember(connection_blocked(count),poly_in) %the principal node's polygon index in Q
                            dist_total(i) = dist_total(i) + deg2km(sqrt((connections(i,2)-x_intersect(i,count,1))^2 + ...
                                (connections(i,3)-y_intersect(i,count,1))^2))*polycost(j);
                            principal = [x_intersect(i,count,1) y_intersect(i,count,1)];
                        else % store the data for the rest of the polygons
                            blockersx(f,:) = polygonx(Q(connection_blocked(count)),:); %store for the shortest path algo the polygons that you'll be using, associated with the 
                            blockersy(f,:) = polygony(Q(connection_blocked(count)),:);
                            polycost_Dijkstra(f) = polycost(Q(connection_blocked(count)));
                            f = f+1;
                        end
               count = count + 1;
                        %else %the other blocking polygons that do not contain the nodes. this will add the distances inside the polygons.
                            %here try to MINIMIZE DISTANCE by either going
                            %around or through the polygon. Use MATLAB
                            %shortest path algorithm to find this distance
                            %dist_total(i) = dist_total(i) + deg2km(sqrt((x_crossings(count,count,1)-x_crossings(count,count,2))^2 + ...
                               % (y_crossings(count,count,1) - y_crossings(count,count,2))^2))*polycost(Q(count)); %in the countth iteration of the set of blocking polys corresp to this connection
                end
                terminal = [connections(i,5) connections(i,6)];
                dist_total(i) = dist_total(i) + Dijkstra(principal, terminal, blockersx, blockersy, polycost_Dijkstra,resolution);
                        %now need to sum up interpolygon distances.
                        %if count == 1 %if the first iteration, then go from principal polygon boundary to boundary of next polygon
                         %   dist_total(i) = dist_total(i) + deg2km(sqrt((x_crossings(count,count,1) - x_crossings(count+1,count+1,1))^2 + ...
                         %       (y_crossings(count,count,1) - y_crossings(count+1,count+1,1))^2));
                            
                       % elseif count ~= length(blocked) %if not a principal or terminal polygon, then compare distances between end of the connection's crossing of polygon count and the beginning of connection's crossing of count+1
                        %    dist_total(i) = dist_total(i) + deg2km(sqrt((x_crossings(count,count,2) - x_crossings(count+1,count+1,1))^2 + ...
                        %        (y_crossings(count,count,2) - y_crossings(count+1,count+1,1))^2));
                       % else
                            %hi
                       % end
                        
                        
                        
                        
                        %else %the other blocking polygons that do not contain the nodes. this will add the distances inside the polygons.
                             
                            %cross = 1;
                            %while cross < length(x_crossings)
                             %   dist_total(i) = dist_total(i) + deg2km(sqrt((x_crossings(i)+x_crossings(i+1))^2 + ...
                              %      (y_crossings(i)+y_crossings(i+1))^2))*polycost(blocked(count));
                               % i = i + 2;
                            %end
                           
                
                
            else
                %the terminal node is not in a poly and there is not a blocking poly
                x_cross = x_intersect(P(connection_blocked),poly_blocker,1);
                y_cross = y_intersect(P(connection_blocked),poly_blocker,1);
            
                dist_total(i) = deg2km(sqrt((connections(i,2)-x_cross)^2 + (connections(i,3)-y_cross)^2))*polycost(j) + ...
                    deg2km(sqrt((x_cross-connections(i,5))^2 + (y_cross - connections(i,6))^2)); %total distance is sum of penalized dist and reg dist
                
            end
        
       
           
                
    elseif ismember(connections(i,4),I) %if the terminal node is in a polygon when the principal node is not
        k = J(find(I==connections(i,4),1)); %find the corresp polygon for the terminal node
        
        if length(connection_blocked) >= 2 %principal node not in poly, terminal node is, is there a blocking poly
            %x_crossings = x_intersect(P(connection_blocked),Q(poly_blocker),:);
            %y_crossings = y_intersect(P(connection_blocked),Q(poly_blocker),:);
            poly_out = find(Q==k); %finding the polygon index in set Q of blocking polys equal to the poly in which the terminal node resides I DON'T THINK THIS FINNA WORK
            count = 1;
            f = 1;
            blockersx = zeros(length(connection_blocked)-1, length(polygonx));
            blockersy = zeros(length(connection_blocked)-1, length(polygony));
            polycost_Dijkstra = zeros(length(connection_blocked)-1, 1);
            while count <= length(connection_blocked) %over all the blocking polygons
                if ismember(connection_blocked(count), poly_out) %if the blocking poly is the home poly. checks for correspondence between the index of blocked connections and the index of the resident polygon from P and Q
                    dist_total(i) = dist_total(i) + deg2km(sqrt((connections(i,5)-x_intersect(i,count,1))^2 + ...
                        (connections(i,6)-y_intersect(i,count,1))^2))*polycost(k);
                    terminal = [x_intersect(i,count,1) y_intersect(i,count,1)];
                else
                    blockersx(f,:) = polygonx(Q(connection_blocked(count)),:);
                    blockersy(f,:) = polygony(Q(connection_blocked(count)),:);
                    polycost_Dijkstra(f) = polycost(Q(connection_blocked(count)));
                    f = f+1;
                 %   dist_total(i) = dist_total(i) + deg2km(sqrt((x_crossings(count,count,1)-x_crossings(count,count,2))^2 + ...
                  %              (y_crossings(count,count,1) - y_crossings(count,count,2))^2))*polycost(Q(count)); %in the countth iteration of the set of blocking polys corresp to this connection
                end
                
                %now need to sum up interpolygon distances.
                %if count == length(blocked) %if the first iteration, then go from terminal polygon boundary to boundary of next polygon
                 %   dist_total(i) = dist_total(i) + deg2km(sqrt((x_crossings(count,count,1) - x_crossings(count-1,count-1,1))^2 + ...
                  %  (y_crossings(count,count,1) - y_crossings(count-1,count-1,1))^2));
                            
                %elseif count ~= length(blocked) %if not a principal or terminal polygon, then compare distances between end of the connection's crossing of polygon count and the beginning of connection's crossing of count+1
                 %   dist_total(i) = dist_total(i) + deg2km(sqrt((x_crossings(count,count,2) - x_crossings(count-1,count-1,1))^2 + ...
                  %      (y_crossings(count,count,2) - y_crossings(count-1,count-1,1))^2));
                %else
                            %hi
                %end
                count = count + 1;
            end
            principal = [connections(i,2) connections(i,3)];
            dist_total(i) = dist_total(i) + Dijkstra(principal, terminal, blockersx, blockersy, polycost_Dijkstra,resolution); %total distance is distance spent inside terminal polygon + shortest path algo from princip node to boundary of terminal node
        else %principal node not in poly, terminal node is, but no blocking poly. then blocked is only the index of the crossage
            x_cross = x_intersect(P(connection_blocked), poly_blocker, 1);
            y_cross = y_intersect(P(connection_blocked), poly_blocker, 1);
            
            dist_total(i) = deg2km(sqrt((connections(i,5)-x_cross)^2 + (connections(i,6)-y_cross)^2))*polycost(k) + ...
                    deg2km(sqrt((x_cross-connections(i,2))^2 + (y_cross - connections(i,3))^2)); %total distance is sum of penalized dist and reg dist
                
        end
        
    %for j =1:length(row_polygon)
    %[x_intersect,y_intersect] = polyxpoly([connections(i,3), connections(i,6)],[connections(i,2), connections(i,5)], polygonx(
    elseif ~isempty(connection_blocked) %neither node is in a polygon, but  there are blocking polygons
        principal = [connections(i,2) connections(i,3)];
        terminal = [connections(i,5) connections(i,6)];
        count = 1;
        
        blockersx = zeros(length(connection_blocked), length(polygonx));
        blockersy = zeros(length(connection_blocked), length(polygony));
        polycost_Dijkstra = zeros(length(connection_blocked), 1);
        while count <= length(connection_blocked)
            blockersx(count,:) = polygonx(Q(connection_blocked(count)),:); %find the data points for the corresponding blocking polygon
            blockersy(count,:) = polygony(Q(connection_blocked(count)),:);
            polycost_Dijkstra(count) = polycost(Q(connection_blocked(count)));
            count = count + 1; 
        end
        dist_total(i) = dist_total(i) + Dijkstra(principal, terminal, blockersx, blockersy, polycost_Dijkstra,resolution);
        
        %for count = 1:length(blocked)
         %   dist_total(i) = dist_total(i) + deg2km(sqrt((x_crossings(count,count,1)-x_crossings(count,count,2))^2 + ...
          %                      (y_crossings(count,count,1) - y_crossings(count,count,2))^2))*polycost(Q(count)); %in the countth iteration of the set of blocking polys corresp to this connection
            %now need interpolygon distances
           % dist_total(i) = dist_total(i) + deg2km(sqrt((x_crossings(count,count,2)-x_crossings(count+1,count+1,1))^2 + ...
            %    (y_crossings(count,count,2)-x_crossings(count+1,count+1,1))^2));
        %end
        
    else 
        %neither node is in a polygon and there are no intersecting polys
        dist_total(i) = deg2km(sqrt((connections(i,2)-connections(i,5))^2+(connections(i,3)-connections(i,6))^2)); %total distance of connection
    end
    i = i + 1; %indexing!
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
toc
end

