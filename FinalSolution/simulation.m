%{
2014 3D Mathematics Coursework
by Calum Brown, Euan McMenemin, Zoltán Tompa


refferences;
% http://uk.mathworks.com/help/matlab/ref/input.html
% http://uk.mathworks.com/help/matlab/matlab_env/save-load-and-delete-workspace-variables.html
% http://uk.mathworks.com/help/matlab/ref/if.html
% http://uk.mathworks.com/help/matlab/ref/try.html
% http://uk.mathworks.com/help/matlab/ref/while.html
% http://uk.mathworks.com/help/matlab/ref/plot.html?searchHighlight=plot
% 
%}%{
2014 3D Mathematics Coursework
by Calum Brown, Euan McMenemin, Zoltán Tompa


refferences;
% http://uk.mathworks.com/help/matlab/ref/input.html
% http://uk.mathworks.com/help/matlab/matlab_env/save-load-and-delete-workspace-variables.html
% http://uk.mathworks.com/help/matlab/ref/if.html
% http://uk.mathworks.com/help/matlab/ref/try.html
% http://uk.mathworks.com/help/matlab/ref/while.html
% http://uk.mathworks.com/help/matlab/ref/plot.html?searchHighlight=plot
% 
%}


%{
    NOTES / KNOWN BUGS / STILL TO BE IMPLEMENTED:
    
    - implement initial firing altitude
    - implement formula for bounce
    - fill new array with after-bounce movement
    - plot new arrays
    - make sure all variables are filled
    - set wall hight back to 5.0

%}

%{
   FORMULAS USED IN CALCULATIONS:

        x = ((m*Vo)/D)*cos(Zrad)*(1-exp((-1*D/m)*T));
        y = (m/D)*(Vo*sin(Zrad)+(m*g/D))*(1-exp(-1*(D/m)*T))-(m*g*T/D);

        Vx = Vo*cos(Zrad)*exp((-1*D/m)*T)
        Vy = Vo*sin(Zrad)+(m*g/D)*(exp(-1*(D/m)*T))-m/D*g)

        angle after bounce = tan^-1(Vwhy/Vwhx)

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
    Wh = 50.0; %wall height (m)
    Wdis = 20.0; %wall distance (m)
    IFH = 0.0; %initial firing altitude (m)
    x = 0.0; %initial X coordinate
    
    CR = 0.58; %coefficient of restitution ()

    Tstep = 0.01; %stepping used for Eulers
    T =0.0; %timer

% - runtime input variables

    Vo = 0.0; %initial velocity (m/s)
    Z = 0.0; %firing angle (deg)

% - calculated variables

    ZW = 0.0; %angle of wall bounce
    ZB = 0.0; %angle after wall bounce

    VxWH = 0.0; %velocity of ball on X axis after hitting the wall - DONE
    VyWH = 0.0; %velocity of ball on Y axis after hitting the wall - DONE
    
    YWH = 0.0; %Y axis value when hitting the wall - DONE
    XWH = 0.0; %X axis value when hitting the wall - DONE

    tWH = 0.0; %time when ball hits wall - DONE
    tGH = 0.0; %time when ball hits ground (from 0 = when wall hit)
    tWHaGH = 0.0; % total fly time
    
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



                   
% -Function - motion before wall hit

        Xb(1,1)= 0;
        Yb(1,1)=IFH+(d/2);
        i=1; %loop iterator
        while x < (Wdis-(d/2)) % 19.5 by deff.
            i=i+1;
            T=(T+Tstep);
            x = ((m*Vo)/D)*cosd(Z)*(1-exp((-1*D/m)*T));
            Xb(i,1)= x;
            y = IFH+(d/2)+(m/D)*(Vo*sind(Z)+(m*g/D))*(1-exp(-1*(D/m)*T))-(m*g*T/D);
            Yb(i,1) = y;

            
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
        XWH = x;
 %-FUNCTION END
 
 
 h = plot (Xb,Yb,'Marker','o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',5);



%- FUNCTION - CHECK IF WALL REACHED AT ALL on the Y %
YWH = y;
if YWH > Wh
   disp(' ');
   out = ['The ball NEVER hit the wall, it flew over it! '];
   disp(out);
   out = [' Its height when passing the wall was: ',num2str(y),' meters'];
   disp(out);
   break;
end

%- FUNCTION END





%- FUNCTION - CALCULATE VELOCITIES WHEN WALL HIT

VxWH = Vo*cosd(Z)*exp((-1*D/m)*tWH); %velocity of ball on X axis after hitting the wall
VyWH = (Vo*sind(Z)+((m/D)*g))*(exp(-1*(D/m)*tWH)-(m/D)*g); %velocity of ball on Y axis after hitting the wall

VW = sqrt((VxWH*cosd(Z)*VxWH*cosd(Z))+((VyWH*sind(Z))*VyWH*sind(Z)));

% FUNCTION TO CALCULATE THE ANGLE OF THE BOUNCE



%  angle at bounce = tan^-1(Vwhy/Vwhx)

ZW = atand((VyWH/VxWH));
ZB = abs(ZW);


 
% -Function - motion after wall hit



        i=1; %loop iterator
        Xa(1,1)= XWH;
        Ya(1,1) = YWH;
        while (y > (d/2) )  % when d/2 it's touching the ground   || (i < 1000)
            
            i=i+1;
            tGH=(tGH+Tstep);
            x = XWH +((m*Vo)/D)*cosd(Z)*(1-exp((-1*D/m)*tGH));
            Xa(i,1)= x;
            %EXPERIMANTAL !!!! NOT ACCTUAL FORMULA (Vo/5) !!!
            %y = (m/D)*(Vo*sind(Z)+(m*g/D))*(1-exp(-1*(D/m)*tGH))-(m*g*tGH/D);
            y = IFH+(d/2)+(m/D)*(Vo*sind(Z)+(m*g/D))*(1-exp(-1*(D/m)*tGH))-(m*g*tGH/D);
            Ya(i,1) = y;
  
                
        end
  
 %-FUNCTION END



% FUNCTION CALCULATE TOTAL FLIGHT TIME
tWHaGH = tWH+tGH;

%function end





% FUNCTION DRAW THE WALL
 
%{
j=1;
isEven = false;
for j = 1:(Wh*10)
    
    if (isEven)
        Wallx(j,1)=(Wdis+0.1)-0.1;
        isEven = false;
    else
       Wallx(j,1)=(Wdis+0.1)+0.1;
       isEven = true;
    end
Wally(j,1)=j*0.1;

end

%}


% draw the calculations

%h = plot (Xa,Ya,'Marker','o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',2);
h = plot (Xb,Yb,Xa,Ya);
%h2 =plot (Xb,Yb);
%createfigure(Xb,Yb);
%{
createfigure(Xa,Ya);

    set(h,'LineWidth',2);
    xlabel('distance (meters)'); ylabel('height (meters)');
    hold on;
    pause(0.001);


%h = plot (Xb,Yb,Wallx,Wally);
%h = plot (Xb,Yb,'Marker','o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',5);
%plot (Xa,Ya,'Marker','o','MarkerFaceColor','red','MarkerEdgeColor','red','MarkerSize',4);


[x,y] = meshgrid(-1.75:.2:3.25);
z = x.*exp(-x.^2-y.^2);

figure
surf(x,y,z)
xlim([-1.75,3.25])
ylim([-1.75,3.25])


%}
 %}
    

