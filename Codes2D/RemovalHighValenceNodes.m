function [EToV,VX,VY,K,Nv] = RemovalHighValenceNodes(EToV,VX,VY,Valence)
% function [EToV,VX,VY,K,Nv] = RemovalHighValenceNodes(EToV,VX,VY)
%
% Removal of high valence nodes from quad mesh
%
% By Allan P. Engsig-Karup, apek@imm.dtu.dk.
%
K = size(EToV,1);
Nv = length(VX);
if Valence<6
    disp('For unstructured quad meshes nodes with')
    disp('valence smaller than 6 can only be removed. Aborting routine.')
    
    return
end
% Remove high valence nodes by doing high valence nodes first then smaller
% etc.
Nfaces = size(EToV,2);
Nnodes = length(VX);
EToVOutput = EToV; % make duplicate to be updated
for n = 1 : Nnodes %TODO: ONLY RUN OVER INTERIOR NODES?? (APEK)
    % Detect valence of current node 
    [i,j] = find(EToV==n);
    M = length(i); % node valence
    if M==Valence
        % extract mesh elements which connect to node
        EToVhigh = EToV(i,:);
        % rotate all vertices so that position of first face is known
        for q = 1 : M
            pos = j(q);
            shift = mod(pos-1,Nfaces);
            shiftvec = mod([0:Nfaces-1]+shift,Nfaces)+1;
            EToVhigh(q,:) = EToVhigh(q,shiftvec);
        end
        % visualize detected node
        quadplot(EToVhigh,VX,VY,'r')
        % Compute connectivity tables
        [EToEhigh,EToFhigh]= tiConnect2Dquad(EToVhigh);
        % make linked list (ring)
        linkedlist = [1]; % element no. 1 in EToEhigh
        for s = 2:M
            connect = linkedlist(s-1);
            linkedlist = [linkedlist; EToEhigh(connect,Nfaces)];
        end
        if length(linkedlist)~=M & unique(linkedlist)<M
            disp('Wrong length of linked list. Aborting current node.')
        else
            % create two new vertices
            Elm2 = linkedlist(2);
            x1 = VX(EToVhigh(Elm2,3));
            y1 = VY(EToVhigh(Elm2,3));
            ElmM = linkedlist(2+round(M/2));
            x2 = VX(EToVhigh(ElmM,3));
            y2 = VY(EToVhigh(ElmM,3));
            % compute two new vertices
            xnew1 = 0.33*(x2-x1)+x1;
            ynew1 = 0.33*(y2-y1)+y1;
            xnew2 = 0.67*(x2-x1)+x1;
            ynew2 = 0.67*(y2-y1)+y1;
            plot(xnew1,ynew1,'ro')
            plot(xnew2,ynew2,'ro')
            % and add to coordinate lists
            Nnodes = length(VX);
            VX = [VX; xnew1; xnew2];
            VY = [VY; ynew1; ynew2];
            % Remove high-valency node
            vh  = n; % high-valency node, global index
            vn1 = Nnodes+1; % new node 1
            vn2 = Nnodes+2; % new node 2
            % modify EToV table and insert new element
            for ii = 1:round(M/2)
                EToVhigh(linkedlist(ii),1) = vn1;
            end
            for ii = round(M/2)+1:M
                EToVhigh(linkedlist(ii),1) = vn2;
            end
            % add new element (always split by 1st element)
            EToVnew = [EToVhigh(linkedlist(1),2) vn1 ...
                EToVhigh(linkedlist(round(M/2)),4) vn2];
            EToVhigh = [EToVhigh; EToVnew];
            quadplot(EToVhigh,VX,VY,'b-')
            % update EToV
            K2 = size(EToVOutput,1);
            EToVOutput(K2+1,:) = zeros(1,Nfaces);
            EToVOutput([i; K2+1],:) = EToVhigh;
        end
    end
end
EToV = EToVOutput;
K  = size(EToV,1);
Nv = length(VX);
return

