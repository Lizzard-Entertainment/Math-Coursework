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

    Tstep = 0.001; %stepping used for Eulers
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

% -- initialisation --

Zrad = degtorad(Z);
y = IFH;

% --



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
    
disp ('all COOL');

% --

% Function - calculate wall hit

        i=0; %loop iterator
        while x < (Wdis-(d/2))
            x = ((m/D)*Vo*cos(Zrad))*(1-exp((-1*D*T)/m));
            T=T+Tstep;
            i=i+1;
        end
        tWH = T;
% FUNCTION END




    