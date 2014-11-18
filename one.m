%input('enter variable name: ',s);

% http://uk.mathworks.com/help/matlab/ref/input.html
% http://uk.mathworks.com/help/matlab/matlab_env/save-load-and-delete-workspace-variables.html
% http://uk.mathworks.com/help/matlab/ref/if.html
% http://uk.mathworks.com/help/matlab/ref/try.html


% - INITIALISE  VARIABLES - %


%mass:
%drag coefficient (C)
%density of air (p)
%initial velocity
%angle
%gravity (g)
%D:

s = 10;
H = zeros(s);


for c = 1:s2
    
    for r = 1:s
        H(r) = sin(c);
    end
end

plot(2s);

prompt = 'please enter initial speed (m/s) ';
V0 = input(prompt);

if V0>0 


prompt2 = 'please enter firing angle (degrees) ';
alpha = input(prompt2);

if (alpha < 90) && (alpha > 0)
    
  disp('values are fine');
    
else
 
     disp('ERROR - entered angle is invalid!'); 
     break;
end
    

else
   disp('ERROR - entered speed is invalid!');
end
    