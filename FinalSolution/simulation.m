%{
2014 3D Mathematics Coursework
by Calum Brown, Euan McMenemin, Zolt�n Tompa


refferences;
% http://uk.mathworks.com/help/matlab/ref/input.html
% http://uk.mathworks.com/help/matlab/matlab_env/save-load-and-delete-workspace-variables.html
% http://uk.mathworks.com/help/matlab/ref/if.html
% http://uk.mathworks.com/help/matlab/ref/try.html
% http://uk.mathworks.com/help/matlab/ref/while.html
% 
%}


% -- initialisation --

clear all;

% --


% --- INITIALISE  VARIABLES --- %

% - preset variables (constants) -

    m = 0.5; %mass of ball (kg)
    d = 1.0; %diameter of the ball (m)
    
    C = 0; %dragCoefficient ()
    p = 0; %density of air ()
    D = 0.07; % air resistance / Drag (Kg/s)

    g = 9.80665; %gravity (m/s^2)
    Wh = 5.0; %wall height (m)
    Wdis = 20.0; %wall distance (m)
    IFH = 1.0; %initial firing altitude (m)
    x = 0.0; %initial X coordinate
    
    CR = 0.58; %coefficient of restitution ()

    Tstep = 0.01; %stepping used for Eulers
    T =0.0; %timer

% - runtime input variables

    Vo = 0.0; %initial velocity (m/s)
    Z = 0.0; %firing angle (deg)

% - calculated variables
    Zrad = 0.0; %firing angle (rad)

    VxWH = 0.0; %velocity of ball on X axis after hitting the wall
    VyWH = 0.0; %velocity of ball on Y axis after hitting the wall
    
    YWH = 0.0; %Y axis value when hitting the wall

    tWH = 0.0; %time when ball hits wall
    tGH = 0.0; %time when ball hits ground
    tAWH = 0.0; %time between ball hits wall and then ground
    
    BdistGHI = 0.0; %ball's distance from initial firing position after hits the ground
    BdistGHW = 0.0; %ball's distance from wall after hits the ground

    
% --- PROGRAM STARTS HERE ---





% -- input + validation --

prompt = 'please enter initial speed of the ball (m/s) ';
Vo = input(prompt);

if Vo<=0
       disp('ERROR - entered speed is invalid!');
       break;
end

prompt2 = 'please enter firing angle (degrees) ';
Z = input(prompt2);

if (Z >= 90) || (Z <= 0)
     disp('ERROR - entered angle is invalid!'); 
     break;
end
    
Zrad = degtorad(Z);




disp ('all COOL');



                   
% -Function - calculate wall hit time

        i=1; %loop iterator
        while x < (Wdis-(d/2)) % 19.5 by deff.
            i=i+1;
            T=(T+Tstep);
            x = ((m*Vo)/D)*cos(Zrad)*(1-exp((-1*D/m)*T));
            Xa(i,1)= x;
            y = (m/D)*(Vo*sin(Zrad)+(m*g/D))*(1-exp(-1*(D/m)*T))-(m*g*T/D);
            Ya(i,1) = y;

            
                   %- FUNCTION - CHECK IF WALL REACHED AT ALL on the X%
                   if y<0
                       disp(' ');
                       out = ['The ball NEVER hit the wall! It landed at ', num2str(x), ' meters away from the firing spot,'];
                       disp(out);
                       out = [num2str(T) , ' seconds after firing.'];
                       disp(out);
                       break;
                   end
    
        end
        tWH = T;
 %-FUNCTION END


%- FUNCTION - CHECK IF WALL REACHED AT ALL on the Y %
YWH = (m/D)*(Vo*sin(Zrad)+(m*g/D))*(1-exp(-1*(D/m)*tWH))-(m*g*tWH/D);
if YWH > Wh
   disp(' ');
   out = ['The ball NEVER hit the wall, it flew over it! '];
   disp(out);
   out = [' Its height when passing the wall was: ',num2str(y)];
   disp(out);
   break;
end

%- FUNCTION END

plot (Xa,Ya);





    