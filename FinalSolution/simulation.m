%%
%{
2014 3D Mathematics Coursework
by Calum Brown, Euan McMenemin, Zolt�n Tompa


refferences;
% http://uk.mathworks.com/help/matlab/ref/input.html
% http://uk.mathworks.com/help/matlab/matlab_env/save-load-and-delete-workspace-variables.html
% http://uk.mathworks.com/help/matlab/ref/if.html
% http://uk.mathworks.com/help/matlab/ref/try.html
% http://uk.mathworks.com/help/matlab/ref/while.html
% http://uk.mathworks.com/help/matlab/ref/plot.html?searchHighlight=plot
% 


    NOTES / KNOWN BUGS / STILL TO BE IMPLEMENTED:
    

    - plot new arrays
    - make sure all variables are filled
    - set wall hight back to 5.0


   FORMULAS USED IN CALCULATIONS:

        x = ((m*Vo)/D)*cos(Zrad)*(1-exp((-1*D/m)*T));
        y = (m/D)*(Vo*sin(Zrad)+(m*g/D))*(1-exp(-1*(D/m)*T))-(m*g*T/D);

        Vx = Vo*cos(Zrad)*exp((-1*D/m)*T)
        Vy = Vo*sin(Zrad)+(m*g/D)*(exp(-1*(D/m)*T))-m/D*g)

        angle after bounce = tan^-1(Vwhy/Vwhx)

%}

%% -- initialisation --

clear all;

% --



%% --- INITIALISE  VARIABLES --- %

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
    CR = 0.58; %coefficient of restitution ()
    Tstep = 0.01; %stepping used for Eulers



% - runtime input variables

    Vo = 0.0; %initial velocity (m/s)
    Z = 0.0; %firing angle (deg)
    
% - runtime input variables 
    x = 0.0; %initial X coordinate
    T =0.0; %timer
    TERMINATE = 0;

% - calculated variables

    ZW = 0.0; %angle of wall bounce - DONE
    ZB = 0.0; %angle at wall bounce - DONE

    VxWH = 0.0; %velocity of ball on X axis after hitting the wall - DONE
    VyWH = 0.0; %velocity of ball on Y axis after hitting the wall - DONE
    
    YWH = 0.0; %Y axis value when hitting the wall - DONE
    XWH = 0.0; %X axis value when hitting the wall - DONE

    tWH = 0.0; %time when ball hits wall - DONE
    tGH = 0.0; %time when ball hits ground (from 0 = when wall hit) - DONE
    tWHaGH = 0.0; % total fly time - DONE
    
    BdistGHI = 0.0; %ball's distance from initial firing position after hits the ground - DONE
    BdistGHW = 0.0; %ball's distance from wall after hits the ground - DONE

    
%% --- PROGRAM STARTS HERE ---

    



%% -- input + validation --

prompt = 'please enter initial speed of the ball (m/s) ';
Vo = input(prompt);

if Vo<=0
       disp('ERROR - entered speed is invalid!');
       TERMINATE=1;
       break;
end

prompt2 = 'please enter firing angle (degrees) ';
Z = input(prompt2);

if (Z >= 90) || (Z <= 0)
     disp('ERROR - entered angle is invalid!');
     TERMINATE=1;
     break;
end
    

                   
%% - Function - motion before wall hit

if (TERMINATE==0)

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
                           TERMINATE=1;
                           break;
                       end

            %% calculate velocities - DEBUG
            %{
            VxWHb(i,1) = Vo*cosd(Z)*exp((-1*D/m)*T); %velocity of ball on X axis after hitting the wall
            VyWHb(i,1) = (Vo*sind(Z)+m/D*g)*(exp(-1*(D/m)*T))-((m/D)*g);
            %}



    end
            tWH = T;
            XWH = x;
    % -FUNCTION END
 
end




%% - FUNCTION - CHECK IF WALL REACHED AT ALL on the Y %
if (TERMINATE==0)

    YWH = y;
    if YWH > Wh
       disp(' ');
       out = ['The ball NEVER hit the wall, it flew over it! '];
       disp(out);
       out = ['Its height when passing the wall was: ',num2str(y),' meters'];
       disp(out);
       TERMINATE=1;
       break;
    end
    
end
%- FUNCTION END





%% - FUNCTION - CALCULATE VELOCITIES WHEN WALL HIT

if (TERMINATE==0)
    
    VxWH = Vo*cosd(Z)*exp((-1*D/m)*tWH);                        %velocity of ball on X axis after hitting the wall
    VyWH = (Vo*sind(Z)+m/D*g)*(exp(-1*(D/m)*tWH))-((m/D)*g);  %velocity of ball on Y axis after hitting the wall

    VW = sqrt((VxWH*cosd(Z)*VxWH*cosd(Z))+((VyWH*sind(Z))*VyWH*sind(Z)));


    %% FUNCTION TO CALCULATE THE ANGLE OF THE BOUNCE

    %  angle at bounce = tan^-1(Vwhy/Vwhx)

    ZW = atand((VyWH/VxWH));
    ZB = abs(ZW);



    %% -Function - motion after wall hit

            i=1; %loop iterator
            Xa(1,1)= XWH;
            Ya(1,1) = YWH;
            while (y > (d/2) )  % when d/2 it's touching the ground   || (i < 1000)

                i=i+1;
                tGH=(tGH+Tstep);
                x = XWH -((m*VxWH)/D)*cosd(ZB)*(1-exp((-1*D/m)*tGH));
                Xa(i,1)= x;
                y=  YWH+(m/D)*((VyWH*CR)*sind(ZB)+(m*g/D))*(1-exp(-1*(D/m)*tGH))-(m*g*tGH/D);
                Ya(i,1) = y;


            end

     %-FUNCTION END



    %% FUNCTION CALCULATE TOTAL FLIGHT TIME
    tWHaGH = tWH+tGH;

    %function end

    %% FUNCTION CALCULATE ball's distances hits the ground

        BdistGHI = x; %ball's distance from initial firing position after hits the ground
        BdistGHW = Wdis-x; %ball's distance from wall after hits the ground

    %function end

    
    %% FUNCTION OUTPUT CALCULATED VALUES
    
       out = ['------------------------------------------------------------'];
       disp(out)
       out = ['The ball hit the wall at the height of ',num2str(YWH),' meters, ',num2str(tWH),' seconds after firing.'];
       disp(out);
       out = ['It bounced back and landed ',num2str(tGH),' seconds later, ',num2str(BdistGHW),' meters away from the wall, and ', num2str(BdistGHI),' meters away from its initial firing position.'];
       disp(out);
       out = [''];
       disp(out);
       out = ['The ball spent ',num2str(tWHaGH),' seconds in the air.' ];
       disp(out);
    
    
    %% FUNCTION DRAW THE WALL

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


    %% draw the calculations

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
    

end