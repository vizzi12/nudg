function plotWall()
        hold on
        xc=100; yc=162.5+37.5/2; zc=5;    % coordinated of the center

        X = [0 0 0 0 0 1; 1 0 1 1 1 1; 1 0 1 1 1 1; 0 0 0 0 0 1];
        Y = [0 0 0 0 1 0; 0 1 0 0 1 1; 0 1 1 1 1 1; 0 0 1 1 1 0];
        Z = [0 0 1 0 0 0; 0 0 1 0 0 0; 1 1 1 0 1 1; 1 1 1 0 1 1];
        
        X = 10*(X-0.5) + xc;
        Y = 37.5*(Y-0.5) + yc;
        Z = 10*(Z-0.5) + zc; 
        
        fill3(X,Y,Z,[0.9290 0.6940 0.1250],'Edgecolor','None');    % draw cube

        xc=100; yc=87.5/2; zc=5;    % coordinated of the center

        X = [0 0 0 0 0 1; 1 0 1 1 1 1; 1 0 1 1 1 1; 0 0 0 0 0 1];
        Y = [0 0 0 0 1 0; 0 1 0 0 1 1; 0 1 1 1 1 1; 0 0 1 1 1 0];
        Z = [0 0 1 0 0 0; 0 0 1 0 0 0; 1 1 1 0 1 1; 1 1 1 0 1 1];           
        
        X = 10*(X-0.5) + xc;
        Y = 87.5*(Y-0.5) + yc;
        Z = 10*(Z-0.5) + zc; 
        
        fill3(X,Y,Z,[0.9290 0.6940 0.1250],'Edgecolor','None');    % draw cube
        return