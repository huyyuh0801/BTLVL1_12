clc;
clear all;
mu0 = 4*pi*1e-7;         
I_current = input('Enter I (A): '); 
Radius = input('Enter radius of the loop (m): ');
Constant = mu0/(4*pi) * I_current;   % Useful constant
NGrid = 31;            % Number of grid points for plots
xMax = Radius;             % Limits for graphics
yMax = xMax*1.5;           % Limits for graphics
fprintf('Field plotted from x = %g m to x = %g m\n',-xMax,xMax);
fprintf('Field plotted from y = %g m to y = %g m\n',-yMax,yMax);
for i=1:NGrid
    xObs(i) = -xMax + (i-1)/(NGrid-1)*(2*xMax); % x values to plot
    yObs(i) = -yMax + (i-1)/(NGrid-1)*(2*yMax); % y values to plot
end
%@ Loop over the segments in the current loop in yz plane
NSegments = 50;
for k=1:NSegments
    %@ Compute location of the endpoints of a segment
    theta1 = 2*pi*(k-1)/NSegments;
    x1 = 0;
    y1 = Radius*cos(theta1);
    z1 = Radius*sin(theta1);
    theta2 = 2*pi*k/NSegments;
    x2 = 0;
    y2 = Radius*cos(theta2);
    z2 = Radius*sin(theta2);
    %@ Compute components of segment vector dl
    dlx(k) = x2-x1;
    dly(k) = y2-y1;
    dlz(k) = z2-z1;
    %@ Compute the location of the midpoint of a segment
    xc(k) = (x2+x1)/2;
    yc(k) = (y2+y1)/2;
    zc(k) = (z2+z1)/2;
end
%@ Loop over all grid points and evaluate B(x,y) on grid
for i=1:NGrid
    for j=1:NGrid

        Bx = 0;  By = 0;  % Initialize B to zero
        %@ Loop over the segments in the loop
        for k=1:NSegments

            %@ Compute components of the r vector (vector between
            %% segment on loop and observation point)
            rx = xObs(j) - xc(k);
            ry = yObs(i) - yc(k);
            rz = -zc(k);             % Observation points are in xy plane

            %@ Compute r^3 from r vector
            r3 = sqrt(rx^2 + ry^2 + rz^2)^3;

            %@ Compute x and y components of cross product dl X r
            dlXr_x = dly(k)*rz - dlz(k)*ry;
            dlXr_y = dlz(k)*rx - dlx(k)*rz;

            %@ Increment sum of x and y components of magnetic field
            Bx = Bx + Constant*dlXr_x/r3;
            By = By + Constant*dlXr_y/r3;
        end

        %@ Compute normalized vectors of magnetic field direction
        BMag = sqrt(Bx^2 + By^2);
        BDirx(i,j) = Bx/BMag;
        BDiry(i,j) = By/BMag;

    end
    fprintf('Calculation %g%% complete\n',100*i/NGrid);
end
fprintf('\nThe magnetic field strength at the center of a circular loop  = %g (T)\n',Constant*2*pi/Radius);
clf;
figure(gcf); 
quiver(xObs,yObs,BDirx,BDiry); 
hold on;
plot(0,Radius,'bo');
plot(0,-Radius,'rx');
title('Magnetic field direction');
xlabel('x');  
ylabel('y');
axis([-Radius Radius -Radius*1.2 Radius*1.5]);
grid on;
hold off;
