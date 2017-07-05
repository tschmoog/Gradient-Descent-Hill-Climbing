function GradientAscentEgSimple()

% plot the landscape
%ezmesh(@ComplexLandscape,[-3 7],[-3 7])
%ezmesh(@SimpleLandscape,[-3 7],[-3 7])
% Enter maximum number of iterations of the algorithm, learning rate and
% mutation range
NumSteps=50;
LRate=0.1;
MaxMutate=1;

% DONE: choose a random starting point with x and y in the range, (-2, 2)
%x = (rand(1)*4)-2;
%y = (rand(1)*4)-2;
%StartPt = [x,y]; 


%Grid points covering the landscape
GridPoints=[-2:0.25:2];
    for i = 1:length(GridPoints)
        for j = 1:length(GridPoints)
            StartPt = [GridPoints(i) GridPoints(j)];
            [MaxReached(i,j), Iters(i,j)] = GradAscent(StartPt, NumSteps, LRate);
        end
    end 
    
    
pcolor(GridPoints, GridPoints, MaxReached);
colorbar;

figure; 
pcolor(GridPoints, GridPoints, Iters);
colorbar;

%number of starting points reacing maximum
%numPoints = sum(MaxReached)


%GradAscent(StartPt,NumSteps,LRate);
%HillClimb(StartPt,NumSteps,MaxMutate);

function[MaxReached, Iters] = GradAscent(StartPt,NumSteps,LRate)
%PauseFlag=1;
hold on;
Iters = 0;
MaxReached = 0;
for i = 1:NumSteps
    %Check to see if optimum value is reached
    %if true, loop ends
    if MaxReached == 1;
        break
    end
    
    Iters = Iters + 1;
    % DONE: Calculate the 'height' at StartPt using SimpleLandscape
    height = SimpleLandscape(StartPt(1),StartPt(2));
    % DONE: Plot a point on the landscape in 3D 
    % use plot3(x,y,z,,'r*','MarkerSize',10)
    % to get a marker you can see well
    plot3(StartPt(1),StartPt(2),height,'r*','MarkerSize',10);

    % DONE: Calculate the gradient at StartPt
    grad = SimpleLandscapeGrad(StartPt(1),StartPt(2));
    % DONE: Calculate the new point and update StartPoint
    NewPt = StartPt + LRate * grad;
    StartPt = NewPt;
    % Make sure StartPt is within the specified bounds
    StartPt = max([StartPt;-2 -2]);
    StartPt = min([StartPt;2 2]);
    
    % Pause to view output
    %if(PauseFlag)
        %x=input('Press return to continue\nor 0 and return to stop pausing\n');
        %if(x==0) PauseFlag=0; end;
    %end
    
    maximumSimple = 4;
   if height == maximumSimple;
       MaxReached = 1;
   end
   
end
hold off

function[MaxReached, Iters] = HillClimb(StartPt,NumSteps,MaxMutate)
%PauseFlag=1;
hold on;
Iters = 0;
MaxReached = 0;
for i = 1:NumSteps
    if MaxReached == 1
        break
    end

    Iters = Iters +1
    % DONE: Calculate the height of StartPt and plot it in 3d
    height = SimpleLandscape(StartPt(1),StartPt(2));
    %Plot in 3D
    plot3(StartPt(1),StartPt(2),height,'r*','MarkerSize',10);

    % Mutate StartPt
    NewPt = Mutate(StartPt,MaxMutate);
    
    % DONE: Make sure NewPt is within the specified bounds
    
    NewPt = max([NewPt;-2 -2]);
    NewPt = min([NewPt;2 2]);
    
    
    % DONE: Calculate the height of the new pt 
    
    newHeight = SimpleLandscape(NewPt(1),NewPt(2));
    
    
    % DONE: decide whether to update StartPt or Not or not
    
    if newHeight > height
        StartPt = NewPt;
    end
    

    % DONE: pause to view output
    
   % if(PauseFlag)
        %x=input('Press return to continue\nor 0 and return to stop pausing\n');
        %if(x==0) PauseFlag=0; end;
    %end
    
   maximumSimple = 4;
   if height == maximumSimple;
       MaxReached = 1;
   end
end
hold off

% DONE: Mutation function
% Returns a mutated point given the old point and the range of mutation
function[NewPt] = Mutate(OldPt,MaxMutate)
% DONE: select a random distance to mutate MutDist in the range 
% (-MaxMutate,MaxMutate)

MutDist = (rand(1)*2*MaxMutate) - MaxMutate;

% DONE: randomly choose which element of OldPt to mutate 
% and mutate it by MutDist
%CoinToss returns 1 or 0
CoinToss = round(rand(1));

if CoinToss == 1
    NewPt(1) = OldPt(1) + MutDist;
    NewPt(2) = OldPt(2);
else 
    NewPt(2) = OldPt(2) + MutDist;
    NewPt(1) = OldPt(1);
end

% simple landscape function
% returns 'height' given (x,y) position
function[z] = SimpleLandscape(x,y)
z=max(1-abs(2*x),0)+x+y;

% gradient of simple landscape function
% returns gradient of SimpleLandscape given (x,y) position
function[g] = SimpleLandscapeGrad(x,y)
if(1-abs(2*x) > 0)
    if(x<0) g(1) = 3;
    elseif(x==0) g(1)=0;
    else g(1) = -1;
    end
else g(1)=1;
end
g(2)=1;

% Function which draws a complex landscape
function DrawComplexLandscape
f = ['4*(1-x)^2*exp(-(x^2) - (y+1)^2)' ... 
         '- 15*(x/5 - x^3 - y^5)*exp(-x^2-y^2)' ... 
         '- (1/3)*exp(-(x+1)^2 - y^2)' ...
        '-1*(2*(x-3)^7 -0.3*(y-4)^5+(y-3)^9)*exp(-(x-3)^2-(y-3)^2)'   ];

ezmesh(f,[-3,7])

% complex landscape function
% returns 'height' given (x,y) position
function[f]=ComplexLandscape(x,y)
f=4*(1-x)^2*exp(-(x^2)-(y+1)^2) -15*(x/5 - x^3 - y^5)*exp(-x^2-y^2) -(1/3)*exp(-(x+1)^2 - y^2)-1*(2*(x-3)^7 -0.3*(y-4)^5+(y-3)^9)*exp(-(x-3)^2-(y-3)^2);

% gradient of complex landscape function
% returns gradient of ComplexLandscape given (x,y) position
function[g]=ComplexLandscapeGrad(x,y)
g(1)=-8*exp(-(x^2)-(y+1)^2)*((1-x)+x*(1-x)^2)-15*exp(-x^2-y^2)*((0.2-3*x^2) -2*x*(x/5 - x^3 - y^5)) +(2/3)*(x+1)*exp(-(x+1)^2 - y^2)-1*exp(-(x-3)^2-(y-3)^2)*(14*(x-3)^6-2*(x-3)*(2*(x-3)^7-0.3*(y-4)^5+(y-3)^9));
g(2)=-8*(y+1)*(1-x)^2*exp(-(x^2)-(y+1)^2) -15*exp(-x^2-y^2)*(-5*y^4 -2*y*(x/5 - x^3 - y^5)) +(2/3)*y*exp(-(x+1)^2 - y^2)-1*exp(-(x-3)^2-(y-3)^2)*((-1.5*(y-4)^4+9*(y-3)^8)-2*(y-3)*(2*(x-3)^7-0.3*(y-4)^5+(y-3)^9));