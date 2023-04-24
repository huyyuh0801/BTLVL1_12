clc;
clear all;
mu0 = 4*pi*1e-7;        %hang so tu  
I_current = input('Enter I (A): '); %cuong do dong dien 
Radius = input('Enter radius of the loop (m): ');   %ban kinh vong day
Constant = mu0/(4*pi) * I_current; 
NGrid = 31;            % Luoi de ve
xMax = Radius;             % Gioi han x
yMax = xMax*1.5;           % Gioi han y
for i=1:NGrid
    xObs(i) = -xMax + (i-1)/(NGrid-1)*(2*xMax); 
    yObs(i) = -yMax + (i-1)/(NGrid-1)*(2*yMax);
end
NSegments = 50; %so doan chia
for k=1:NSegments
    theta1 = 2*pi*(k-1)/NSegments;
    x1 = 0;
    y1 = Radius*cos(theta1);
    z1 = Radius*sin(theta1);
    theta2 = 2*pi*k/NSegments;
    x2 = 0;
    y2 = Radius*cos(theta2);
    z2 = Radius*sin(theta2);
    dlx(k) = x2-x1;
    dly(k) = y2-y1;
    dlz(k) = z2-z1;
    xc(k) = (x2+x1)/2;
    yc(k) = (y2+y1)/2;
    zc(k) = (z2+z1)/2;
end
for i=1:NGrid
    for j=1:NGrid

        Bx = 0;  By = 0;
        for k=1:NSegments
            rx = xObs(j) - xc(k);
            ry = yObs(i) - yc(k);
            rz = -zc(k); 

            r3 = sqrt(rx^2 + ry^2 + rz^2)^3;

            dlXr_x = dly(k)*rz - dlz(k)*ry;
            dlXr_y = dlz(k)*rx - dlx(k)*rz;

            Bx = Bx + Constant*dlXr_x/r3;
            By = By + Constant*dlXr_y/r3;
        end

        BMag = sqrt(Bx^2 + By^2);
        BDirx(i,j) = Bx/BMag;
        BDiry(i,j) = By/BMag;

    end
end
fprintf('\nThe magnetic field strength at the center of a circular loop  = %g (T)\n',Constant*2*pi/Radius);
clf;
figure(gcf); 
quiver(xObs,yObs,BDirx,BDiry);  %Ve duong suc tu
hold on;
plot(0,Radius,'bo');
plot(0,-Radius,'rx');
title('Magnetic field direction');
xlabel('x');  
ylabel('y');
axis([-Radius Radius -Radius*1.2 Radius*1.5]);
grid on;
hold off;
