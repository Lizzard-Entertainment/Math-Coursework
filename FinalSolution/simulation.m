%{
2014 3D Mathematics Coursework
by Calum Brown, Euan McMenemin, Zoltán Tompa


refferences;
% http://uk.mathworks.com/help/matlab/ref/input.html
% http://uk.mathworks.com/help/matlab/matlab_env/save-load-and-delete-workspace-variables.html
% http://uk.mathworks.com/help/matlab/ref/if.html
% http://uk.mathworks.com/help/matlab/ref/try.html

%}


% --- INITIALISE  VARIABLES --- %

% - preset variables (constants) -

    m = 0.5; %mass of ball (kg)
    d = 1; %diameter of the ball (m)
    
    C = 0; %dragCoefficient ()
    p = 0; %density of air ()
    D = 0.07; % air resistance / Drag (multiplier)

    g = 9.80665; %gravity (m/s^2)
    Wh = 5; %wall height (m)
    Wdis = 20; %wall distance (m)
    CR = 0.58; %coefficient of restitution ()

    Tstep = 0.001; %stepping used for Eulers

% - runtime input variables

    Vo = 0; %initial velocity (m/s)
    Z = 0; %fining angle (deg)

% - calculated variables
    VxWH = 0; %velocity of ball on X axis after hitting the wall
    VyWH = 0; %velocity of ball on Y axis after hitting the wall
    
    YWH = 0; %Y axis value when hitting the wall

    tWH = 0; %time when ball hits wall
    tGH = 0; %time when ball hits ground
    tAWH = 0; %time between ball hits wall and then ground
    
    BdistGHI = 0; %ball's distance from initial firing position after hits the ground
    BdistGHW = 0; %ball's distance from wall after hits the ground

    
% --- PROGRAM STARTS HERE ---

% - input + validation -

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


%plot([1 5], [0 4]);
    