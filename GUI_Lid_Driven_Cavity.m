% Lid-Driven Cavity CFD Code
% Written by: JoshTheEngineer
% YouTube   : https://www.youtube.com/user/joshtheengineer
% Website   : www.joshtheengineer.com
% Started   : 04/28/2013
% Updated   : Unknown

function varargout = GUI_Lid_Driven_Cavity(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Lid_Driven_Cavity_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Lid_Driven_Cavity_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUI_Lid_Driven_Cavity is made visible.
function GUI_Lid_Driven_Cavity_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

% Initialize data
N = 11;                 assignin('base','N',N);
cornerStretch = 0;      assignin('base','cornerStretch',cornerStretch);
r = 1;                  assignin('base','r',r);
Re = 100;               assignin('base','Re',Re);
maxIter = 250;          assignin('base','maxIter',maxIter);
eps = 0;                assignin('base','eps',eps);
CFL = 0.1;              assignin('base','CFL',CFL);

toPlot = 'Streamlines';
toPlotNum = 1;
assignin('base','toPlot',toPlot);
assignin('base','toPlotNum',toPlotNum);

set(handles.editStretchRatio,'Enable','off');  % Grid stretching disabled


% --- Outputs from this function are returned to the command line.
function varargout = GUI_Lid_Driven_Cavity_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% --------------------------- INITIALIZATION ---------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %

% EDIT -------------------- Grid Nodes ---------------------------------- %
function editGridNodes_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% EDIT ------------------ Reynolds Number ------------------------------- %
function editRe_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% EDIT ----------------- Maximum Iterations ----------------------------- %
function editMaxIter_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% EDIT ----------------- Dissipation (epsilon) -------------------------- %
function editEpsilon_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% EDIT ------------------ Stretching Ratio ------------------------------ %
function editStretchRatio_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% EDIT ------------------------ CFL ------------------------------------- %
function editCFL_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% POP ---------------------- What to Plot ------------------------------- %
function popPlot_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------- CALLBACKS ------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %

% EDIT -------------------- Grid Nodes ---------------------------------- %
function editGridNodes_Callback(hObject, eventdata, handles)
N = str2num(get(handles.editGridNodes,'String'));
assignin('base','N',N);

% EDIT ------------------ Reynolds Number ------------------------------- %
function editRe_Callback(hObject, eventdata, handles)
Re = str2num(get(handles.editRe,'String'));
assignin('base','Re',Re);

% EDIT ----------------- Maximum Iterations ----------------------------- %
function editMaxIter_Callback(hObject, eventdata, handles)
maxIter = str2num(get(handles.editMaxIter,'String'));
assignin('base','maxIter',maxIter);

% EDIT ----------------- Dissipation (epsilon) -------------------------- %
function editEpsilon_Callback(hObject, eventdata, handles)
eps = str2num(get(handles.editEpsilon,'String'));
assignin('base','eps',eps);

% EDIT ------------------ Stretching Ratio ------------------------------ %
function editStretchRatio_Callback(hObject, eventdata, handles)
r = str2num(get(handles.editStretchRatio,'String'));
assignin('base','r',r);

% EDIT ------------------------ CFL ------------------------------------- %
function editCFL_Callback(hObject, eventdata, handles)
CFL = str2num(get(handles.editCFL,'String'));
assignin('base','CFL',CFL);

% RADIO ------------------ Grid Stretching ------------------------------ %
function radioGridStretch_Callback(hObject, eventdata, handles)
if (get(hObject,'Value') == get(hObject,'Max'))
    set(handles.editStretchRatio,'Enable','on');   % Grid stretching enabled
    cornerStretch = 1;
else
    set(handles.editStretchRatio,'Enable','off');  % Grid stretching disabled
    cornerStretch = 0;
end
assignin('base','cornerStretch',cornerStretch);

% POP ---------------------- What to Plot ------------------------------- %
function popPlot_Callback(hObject, eventdata, handles)
toPlot    = evalin('base','toPlot');
toPlotNum = evalin('base','toPlotNum');
xx        = evalin('base','xx');
yy        = evalin('base','yy');
x         = evalin('base','x');
y         = evalin('base','y');
Pnew      = evalin('base','Pnew');
unew      = evalin('base','unew');
vnew      = evalin('base','vnew');
uold      = evalin('base','uold');
vold      = evalin('base','vold');
N         = evalin('base','N');
Re        = evalin('base','Re');
r         = evalin('base','r');
eps       = evalin('base','eps');
ErrorTot  = evalin('base','ErrorTot');

popContents = cellstr(get(hObject,'String'));                                 
toPlot = popContents{get(hObject,'Value')};

if (strcmp(toPlot,'Streamlines') == 1)
    toPlotNum = 1;
    set(handles.pushAddStreamlines,'Enable','on');    % Push button enabled
elseif (strcmp(toPlot,'U-Vel Ghia Data') == 1)
    toPlotNum = 2;
    set(handles.pushAddStreamlines,'Enable','off');    % Push button enabled
elseif (strcmp(toPlot,'V-Vel Ghia Data') == 1)
    toPlotNum = 3;
    set(handles.pushAddStreamlines,'Enable','off');    % Push button enabled
elseif (strcmp(toPlot,'Grid') == 1)
    toPlotNum = 4;
    set(handles.pushAddStreamlines,'Enable','off');    % Push button enabled
elseif (strcmp(toPlot,'Error') == 1)
    toPlotNum = 5;
    set(handles.pushAddStreamlines,'Enable','off');    % Push button enabled
elseif (strcmp(toPlot,'Quiver') == 1)
    toPlotNum = 6;
    set(handles.pushAddStreamlines,'Enable','on');    % Push button enabled
end
assignin('base','toPlotNum',toPlotNum);

if (toPlotNum == 1)
    axes(handles.plotData);
    cla;
    hold on;
    % Main vortex (Zone 1)
    s = streamline(xx,yy,unew,vnew,0.88,0.1); set(s,'Color','black');
    s = streamline(xx,yy,unew,vnew,0.02,0.02); set(s,'Color','black');
    s = streamline(xx,yy,unew,vnew,0.7,0.7); set(s,'Color','black');
    s = streamline(xx,yy,unew,vnew,0.97,0.9); set(s,'Color','black');
    s = streamline(xx,yy,unew,vnew,0.97,0.8); set(s,'Color','black');
    s = streamline(xx,yy,unew,vnew,0.97,0.7); set(s,'Color','black');
    s = streamline(xx,yy,unew,vnew,0.97,0.6); set(s,'Color','black');
    s = streamline(xx,yy,unew,vnew,0.97,0.5); set(s,'Color','black');
    s = streamline(xx,yy,unew,vnew,0.97,0.4); set(s,'Color','black');
    s = streamline(xx,yy,unew,vnew,0.97,0.3); set(s,'Color','black');
    s = streamline(xx,yy,unew,vnew,0.97,0.2); set(s,'Color','black');
    s = streamline(xx,yy,unew,vnew,0.97,0.1); set(s,'Color','black');
    % Right vortex (Zone 2)
    s = streamline(xx,yy,unew,vnew,0.97,0.05); set(s,'Color','black');
    s = streamline(xx,yy,unew,vnew,0.93,0.05); set(s,'Color','black');
    s = streamline(xx,yy,unew,vnew,0.94,0.05); set(s,'Color','black');
    % Left vortex (Zone 3)
    s = streamline(xx,yy,unew,vnew,0.03,0.03); set(s,'Color','black');
    s = streamline(xx,yy,unew,vnew,0.04,0.06); set(s,'Color','black');
    s = streamline(xx,yy,unew,vnew,0.02,0.06); set(s,'Color','black');
    s = streamline(xx,yy,unew,vnew,0.02,0.05); set(s,'Color','black');
    title(['Streamline: N = ',num2str(N),', Re = ',num2str(Re),...
                ', r = ',num2str(r),', eps_d_i_s_s = ',num2str(eps)])
    xlabel('X Distance');
    ylabel('Y Distance');
    axis([0 1 0 1]);
elseif (toPlotNum == 2)
    axes(handles.plotData);
    cla;
    hold on;
    mid = ceil(N/2);        % Centerline node index
    xData = fliplr(x);      % Flip the x data left-right
    uData = uold(:,mid);    % Take the centerline values of u velocity
    yData = y;              
    vData = vold(mid,:);    % Take the centerline values of v velocity

    % Table 1 - U Velocity
    valY = [1 .9766 .9688 .9609 .9531 .8516 .7344 .6172 .5 .4531 .2813 ...
                .1719 .1016 .0703 .0625 .0547 0]';
    uy100 = [1 .84123 .78871 .73722 .68717 .23151 .00332 -.13641 -.20581 ...
                -.21090 -.15662 -.10150 -.06434 -.04775 -.04192 -.03717 0]';
    uy400 = [1 0.75837 0.68439 0.61756 0.55892 0.29093 0.16256 0.02135 ...
                -.11477 -.17119 -.32726 -.24299 -.14612 -.10338 -.09266...
                    -.08186 0]';
    uy1000 = [1 .65928 .57492 .51117 .46604 .33304 .18719 .05702 -.06080...
                    -.10648 -.27805 -.38289 -.29730 -.22220 -.20196 -.18109...
                        0]';
    if (Re == 100); plot(valY,uy100,'ro'); end;
    if (Re == 400); plot(valY,uy400,'ro'); end;
    if (Re == 1000); plot(valY,uy1000,'ro'); end;
    plot(flipud(xData),uData,'k-');
    title('U Velocity along Vertical Line at Mid-Cavity');
    xlabel('Cavity Height');
    ylabel('U Velocity');
    axis([0 1 -0.5 1]);
    grid on;
elseif (toPlotNum == 3)
    axes(handles.plotData);
    cla;
    hold on;
    mid = ceil(N/2);        % Centerline node index
    vData = vold(mid,:);    % Take the centerline values of v velocity
    valX = [1 .9688 .9609 .9531 .9453 .9063 .8594 .8047 .5 .2344 .2266 .1563 ...
            .0938 .0781 .0703 .0625 0]';
    vx100 = [0 -.05906 -.07391 -0.08864 -.10313 -.16914 -.22445 -.24533 ...
                    .05454 .17527 .17507 .16077 .12317 .10890 .10091 .09233 0]'; 
    vx400 = [0 -.12146 -.15663 -.19254 -.22847 -.23827 -.44993 -.38598 ...
                .05186 .30174 .30203 .28124 .22965 .20920 .19713 .18360 0]';
    vx1000 = [0 -.21388 -.27669 -.33714 -.39188 -.51550 -.42665 -.31966 ...
                    .02526 .32235 .33075 .37095 .32627 .30353 .29012 .27485 0]';
    if (Re == 100); plot(valX,vx100,'bo'); end;
    if (Re == 400); plot(valX,vx400,'bo'); end;
    if (Re == 1000); plot(valX,vx1000,'bo'); end;
    plot(y,vData,'k-');
    title('V Velocity along Horizontal Line at Mid-Cavity');
    xlabel('Cavity Length');
    ylabel('V Velocity');
    axis([0 1 -0.4 0.4]);
    grid on;
elseif (toPlotNum == 4)
    axes(handles.plotData);
    cla;
    plot(xx,yy,'rx');
    xlabel('X Location');
    ylabel('Y Location');
    axis([0 1 0 1]);
    grid on;
elseif (toPlotNum == 5)
    axes(handles.plotData);
    cla;
    semilogy(ErrorTot);
    title(['N = ',num2str(N),', Re = ',num2str(Re),...
                ', r = ',num2str(r),', eps_d_i_s_s = ',num2str(eps)])
    xlabel('Iterations');
    ylabel('log(Error)');
    axis auto;
    grid on;
elseif (toPlotNum == 6)
    axes(handles.plotData);
    cla;
    quiver(xx,yy,unew,vnew);
    title(['N = ',num2str(N),', Re = ',num2str(Re),...
                ', r = ',num2str(r),', eps_d_i_s_s = ',num2str(eps)])
    xlabel('X Grid Point');
    ylabel('Y Grid Point');
    axis([0 1 0 1]);
    grid on;
end

% PUSH --------------------- Exit the GUI ------------------------------- %
function pushExit_Callback(hObject, eventdata, handles)
clc;                            % Clear command window
evalin('base','clear all');     % Clear all base workspace variables
delete(handles.figure1);        % Delete the figure

% PUSH ------------------- COMPUTE SOLTUION ----------------------------- %
function pushCompute_Callback(hObject, eventdata, handles)
N       = evalin('base','N');
r       = evalin('base','r');
Re      = evalin('base','Re');
maxIter = evalin('base','maxIter');
eps     = evalin('base','eps');
CFL     = evalin('base','CFL');
cornerStretch = evalin('base','cornerStretch');

% Length and Width of Cavity
xlo = 0;                    % Lower x bound
xhi = 1;                    % Upper x bound
ylo = 0;                    % Lower y bound
yhi = 1;                    % Upper y bound

if (cornerStretch == 0)     % If stretching from only bottom left corner
    % X Grid Points
    x = zeros(1,N);             % Initialize to the number of nodes defined
    x(1) = xlo;                 % Set first node equal to lowest value
    i = 2;                      % Absolute index for inner X grid nodes
    for k = 0:1:N-2             % Loop for inner X nodes
        x(i) = x(i-1) + r^k;    % Geometric stretching node definition
        i = i + 1;              % Increment absolute X node counter
    end
    xFactor = x(N)/xhi;     % Normalization factor for grid
    x = x/xFactor;          % Normalize the grid to the maximum X value

    % Y Grid Points
    y = zeros(1,N);             % Initialize to the number of nodes defined
    y(1) = ylo;                 % Set first node equatl to lowest value
    j = 2;                      % Absolute index for inner Y grid nodes
    for k = 0:1:N-2             % Loop for inner Y nodes
        y(j) = y(j-1) + r^k;    % Geometric stretching node definition
        j = j + 1;              % Increment absolute Y node counter
    end
    yFactor = y(N)/yhi;     % Normalization factor for grid
    y = y/yFactor;          % Normalize the grid to the maximum Y value
    
elseif (cornerStretch == 1)
    % X Grid Points
    x = zeros(1,N);             % Initialize all values to zero
    x(1) = xlo;                 % Set the first point to lowest x value
    half = ceil(N/2);           % Midpoint of the number of nodes (use odd)
    k1 = 0;                     % Stretching ratio exponent
    for i = 2:1:N                   % For every node other than first node
        if (i <= half)              % First half of X nodes in array
            x(i) = x(i-1) + r^k1;   % Increasing stretching geometrically
            k1 = k1 + 1;            % Increment stretching exponent
        elseif (i > half)               % Second half of nodes
            x(i) = x(i-1) + r^(k1-1);   % Decrease stretching geometrically
            k1 = k1 - 1;                % Decrement stretching exponent
        end
    end
    xFactor = x(N)/xhi;         % Normalization factor for grid
    x = x/xFactor;              % Normalize the grid to maximum X value

    % Y Grid Points
    y = zeros(1,N);             % Initialize all values to zero
    y(1) = ylo;                 % Set the first point to lowest y value
    k1 = 0;                     % Stretching ratio exponent
    for j = 2:1:N                   % For every node other than first node
        if (j <= half)              % First half of Y nodes in array
            y(j) = y(j-1) + r^k1;   % Increasing stretching geometrically
            k1 = k1 + 1;            % Increment stretching exponent
        elseif (j > half)               % Second half of nodes
            y(j) = y(j-1) + r^(k1-1);   % Decrease stretching geometrically
            k1 = k1 - 1;                % Decrement stretching exponent
        end
    end
    yFactor = y(N)/yhi;         % Normalization factor for grid
    y = y/yFactor;              % Normalize the grid to maximum Y value
end
x = x';     % Transpose X grid values for easier viewing
y = y';     % Transpose Y grid values for easier viewing

% Generate the meshed grid for use later
[xx,yy] = meshgrid(x,y);
assignin('base','xx',xx);
assignin('base','yy',yy);
assignin('base','x',x);
assignin('base','y',y);

% Find smallest dx and dy grid values
for i = 1:1:N-1
    dx(i) = x(i+1) - x(i);
    dy(i) = y(i+1) - y(i);
end
dxMin = min(dx);
dyMin = min(dy);

%% Pre-Allocate Matrices

% Metrics of the Inverse Transformation
x_zeta = zeros(N,1);
y_eta  = zeros(N,1);

% Contravariant velocities
U = zeros(N,N);
V = zeros(N,N);

%% Calculate Metrics

% Defined as unity in computational domain
dZeta = 1;  assignin('base','dZeta',dZeta);
dEta = 1;   assignin('base','dEta',dEta);

% Metrics of Inverse Transformation
for col = 1:1:N
    for row = 1:1:N
        %----------
        if (col == 1)
            x_zeta(row,col) = (-3*x(col)+4*x(col+1)-x(col+2))/(2*dZeta);
        elseif (col == N)
            x_zeta(row,col) = (3*x(col)-4*x(col-1)+x(col-2))/(2*dZeta);
        elseif (col ~= 1 && col ~= N)
            x_zeta(row,col) = (x(col+1)-x(col-1))/(2*dZeta);
        end
        %----------
        if (row == 1)
            y_eta(row,col)  = (-3*y(row)+4*y(row+1)-y(row+2))/(2*dEta);
        elseif (row == N)
            y_eta(row,col)  = (3*y(row)-4*y(row-1)+y(row-2))/(2*dEta);
        elseif (row ~= 1 && row ~= N)
            y_eta(row,col)  = (y(row+1)-y(row-1))/(2*dEta);
        end
    end
end

J = (x_zeta.*y_eta).^-1;    % Jacobian
zeta_x = J.*y_eta;          % Metric of the transformation
eta_y = J.*x_zeta;          % Metric of the transformation
g11 = zeta_x.^2;            % (11) component of metric tensor
g22 = eta_y.^2;             % (22) component of metric tensor

assignin('base','J',J);
assignin('base','zeta_x',zeta_x);
assignin('base','eta_y',eta_y);
assignin('base','g11',g11);
assignin('base','g22',g22);

%% Initial and Boundary Conditions

% Initial Conditions - value of all primitive variables is zero
Pinit = zeros(N,N);
uinit = zeros(N,N);
vinit = zeros(N,N);

% Boundary Condition
velLid = 1;             % Lid velocity [m/s]
uinit(1,:) = velLid;    % Set the u velocity at the lid

assignin('base','velLid',velLid);

%% Time Stepping

Pold = Pinit;   % [25x25] Pressure from initial conditions
uold = uinit;   % [25x25] X-velocity from initial conditions
vold = vinit;   % [25x25] Y-velocity from initial conditions

% Time step variables
dt = dxMin*CFL;            % Time step based on smallest dx
epsError = 1e-8;           % Convergence criterion
assignin('base','dt',dt);
assignin('base','epsError',epsError);

for iter = 1:1:maxIter
    
    % Previous values for use in error calculation
    Pprev = Pold;
    uprev = uold;
    vprev = vold;
    
    % Show iteration in GUI
    set(handles.textIterShow,'String',num2str(iter));
    set(handles.textMaxIterShow,'String',num2str(maxIter));
    drawnow();
    
    % Contravariant velocities
    U = uold.*zeta_x;
    V = vold.*eta_y;
    
    % *** Function ***
    % - Calculate spectral radius of the grid
    % - Returns an [NxN] matrix
    [pA1,pA2] = SPECTRAL_RADIUS3(U,V,g11,g22,J,N);
    
    % *** Function ***
    % - Calculate spectral radius values at the half nodes
    % - pA1_hfnd is [NxN-1] matrix
    % - pA2_hfnd is [N-1xN] matrix
    [pA1_hfnd,pA2_hfnd] = HALFNODE3(pA1,pA2,N);
    
    % *** Function ***
    % - Calculate dissipation values at half nodes
    [D1c,D1u,D1v,D2c,D2u,D2v] = DISSIPATION3(eps,pA1_hfnd,pA2_hfnd,Pold,...
                                                uold,vold,N);
    
    % Runge-Kutta 4th Order Time Stepping
    for l = 4:-1:1
        [RHSc] = CONTINUITY3(N,uold,vold,zeta_x,eta_y,J,dZeta,dEta,...
                                D1c,D2c);
        [RHSu] = XMOMENTUM3(N,uold,vold,Pold,zeta_x,eta_y,J,dZeta,...
                                dEta,Re,D1u,D2u);
        [RHSv] = YMOMENTUM3(N,uold,vold,Pold,zeta_x,eta_y,J,dZeta,...
                                dEta,Re,D1v,D2v);
        P2 = Pold(2:N-1,2:N-1) + (1/l)*dt*RHSc;
        u2 = uold(2:N-1,2:N-1) + (1/l)*dt*RHSu;
        v2 = vold(2:N-1,2:N-1) + (1/l)*dt*RHSv;
        [Pold,uold,vold] = BOUNDARY_CONDITIONS3(P2,u2,v2,velLid,N);
        Pold = Pold - Pold(2,2);
    end
    
    PRestart = Pold;        assignin('base','PRestart',PRestart);
    uRestart = uold;        assignin('base','uRestart',uRestart);
    vRestart = vold;        assignin('base','vRestart',vRestart);
    
    % Check for NaNs
    nanP = isnan(Pold);
    nanu = isnan(uold);
    nanv = isnan(vold);
    checkP = max(max(nanP));
    checku = max(max(nanu));
    checkv = max(max(nanv));
    if (checkP == 1 || checku == 1 || checkv == 1)
        fprintf('Error! -> NaN\nBroke at: %6i iterations\n',iter);
        if (checkP == 1); fprintf('P had NaN\n'); end
        if (checku == 1); fprintf('u had NaN\n'); end
        if (checkv == 1); fprintf('v had NaN\n'); end
        break;
    end 
    
    % Check convergence of velocity field
    errTot = 0;
    erru = (uprev-uold).^2;
    errv = (vprev-vold).^2;
    err  = (erru+errv).^(1/2);
    for row = 1:1:N
        for col = 1:1:N
            errTot = errTot + err(row,col);
        end
    end 
    ErrorTot(iter) = (1/N^2)*errTot;
    ErrorTot = ErrorTot';
    set(handles.textErrorShow,'String',num2str(ErrorTot(iter),'%e'));
    drawnow();
    if (min(ErrorTot) < epsError)
        fprintf('Total Error Calculation\n');
        fprintf('r:    %2.2f\nN:    %3i\nError: %3.7s\n',r,N,min(ErrorTot));
        fprintf('Converged at: %6i iteration\n',iter);
        break;
    end

end

vold = -vold;

Pnew = flipud(Pold);
unew = flipud(uold);
vnew = flipud(vold);

% Restart iteration values
iterRestart = maxIter+1;
totalIters  = maxIter;
assignin('base','iterRestart',iterRestart);
assignin('base','totalIters',totalIters);

% Send final values back into base workspace
assignin('base','ErrorTot',ErrorTot);
assignin('base','Pnew',Pnew);
assignin('base','unew',unew);
assignin('base','vnew',vnew);
assignin('base','uold',uold);
assignin('base','vold',vold);

axes(handles.plotData);
cla;
hold on;
% Main vortex (Zone 1)
s = streamline(xx,yy,unew,vnew,0.88,0.1); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.02,0.02); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.7,0.7); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.97,0.9); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.97,0.8); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.97,0.7); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.97,0.6); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.97,0.5); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.97,0.4); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.97,0.3); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.97,0.2); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.97,0.1); set(s,'Color','black');
% Right vortex (Zone 2)
s = streamline(xx,yy,unew,vnew,0.97,0.05); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.93,0.05); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.94,0.05); set(s,'Color','black');
% Left vortex (Zone 3)
s = streamline(xx,yy,unew,vnew,0.03,0.03); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.04,0.06); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.02,0.06); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.02,0.05); set(s,'Color','black');
title(['Streamline: N = ',num2str(N),', Re = ',num2str(Re),...
            ', r = ',num2str(r),', eps_d_i_s_s = ',num2str(eps)])
xlabel('X Distance');
ylabel('Y Distance');
axis([0 1 0 1]);

% PUSH -------------------- Continue Simulation ------------------------- %
function pushContinueSim_Callback(hObject, eventdata, handles)
xx       = evalin('base','xx');
yy       = evalin('base','yy');
N        = evalin('base','N');
r        = evalin('base','r');
Re       = evalin('base','Re');
maxIter  = evalin('base','maxIter');
eps      = evalin('base','eps');
velLid   = evalin('base','velLid');
uold     = evalin('base','uRestart');
vold     = evalin('base','vRestart');
Pold     = evalin('base','PRestart');
zeta_x   = evalin('base','zeta_x');
eta_y    = evalin('base','eta_y');
g11      = evalin('base','g11');
g22      = evalin('base','g22');
J        = evalin('base','J');
dZeta    = evalin('base','dZeta');
dEta     = evalin('base','dEta');
dt       = evalin('base','dt');
epsError = evalin('base','epsError');
ErrorTot = evalin('base','ErrorTot');
iterRestart = evalin('base','iterRestart'); % 251
totalIters  = evalin('base','totalIters');  % 250

% Update iteration bounds for restart simulation
startIter = iterRestart;            assignin('base','startIter',startIter);
endIter   = totalIters+maxIter;     assignin('base','endIter',endIter);

for iter = startIter:1:endIter
    
    % Previous values for use in error calculation
    Pprev = Pold;
    uprev = uold;
    vprev = vold;
    
    % Show iteration in GUI
    set(handles.textIterShow,'String',num2str(iter));
    set(handles.textMaxIterShow,'String',num2str(endIter));
    drawnow();
    
    % Contravariant velocities
    U = uold.*zeta_x;
    V = vold.*eta_y;
    
    % *** Function ***
    % - Calculate spectral radius of the grid
    % - Returns an [NxN] matrix
    [pA1,pA2] = SPECTRAL_RADIUS3(U,V,g11,g22,J,N);
    
    % *** Function ***
    % - Calculate spectral radius values at the half nodes
    % - pA1_hfnd is [NxN-1] matrix
    % - pA2_hfnd is [N-1xN] matrix
    [pA1_hfnd,pA2_hfnd] = HALFNODE3(pA1,pA2,N);
    
    % *** Function ***
    % - Calculate dissipation values at half nodes
    [D1c,D1u,D1v,D2c,D2u,D2v] = DISSIPATION3(eps,pA1_hfnd,pA2_hfnd,Pold,...
                                                uold,vold,N);
    
    % Runge-Kutta 4th Order Time Stepping
    for l = 4:-1:1
        [RHSc] = CONTINUITY3(N,uold,vold,zeta_x,eta_y,J,dZeta,dEta,...
                                D1c,D2c);
        [RHSu] = XMOMENTUM3(N,uold,vold,Pold,zeta_x,eta_y,J,dZeta,...
                                dEta,Re,D1u,D2u);
        [RHSv] = YMOMENTUM3(N,uold,vold,Pold,zeta_x,eta_y,J,dZeta,...
                                dEta,Re,D1v,D2v);
        P2 = Pold(2:N-1,2:N-1) + (1/l)*dt*RHSc;
        u2 = uold(2:N-1,2:N-1) + (1/l)*dt*RHSu;
        v2 = vold(2:N-1,2:N-1) + (1/l)*dt*RHSv;
        [Pold,uold,vold] = BOUNDARY_CONDITIONS3(P2,u2,v2,velLid,N);
        Pold = Pold - Pold(2,2);
    end
    
    % Check for NaNs
    nanP = isnan(Pold);
    nanu = isnan(uold);
    nanv = isnan(vold);
    checkP = max(max(nanP));
    checku = max(max(nanu));
    checkv = max(max(nanv));
    if (checkP == 1 || checku == 1 || checkv == 1)
        fprintf('Error! -> NaN\nBroke at: %6i iterations\n',iter);
        if (checkP == 1); fprintf('P had NaN\n'); end
        if (checku == 1); fprintf('u had NaN\n'); end
        if (checkv == 1); fprintf('v had NaN\n'); end
        break;
    end 
    
    % Check convergence of velocity field
    errTot = 0;
    erru = (uprev-uold).^2;
    errv = (vprev-vold).^2;
    err  = (erru+errv).^(1/2);
    for row = 1:1:N
        for col = 1:1:N
            errTot = errTot + err(row,col);
        end
    end 
    ErrorTot(iter) = (1/N^2)*errTot;
    ErrorTot = ErrorTot';
    set(handles.textErrorShow,'String',num2str(ErrorTot(iter),'%e'));
    drawnow();
    if (min(ErrorTot) < epsError)
        fprintf('Total Error Calculation\n');
        fprintf('r:    %2.2f\nN:    %3i\nError: %3.7s\n',r,N,min(ErrorTot));
        fprintf('Converged at: %6i iteration\n',iter);
        break;
    end

end

% Untouched restart values
assignin('base','uRestart',uold);
assignin('base','vRestart',vold);
assignin('base','PRestart',Pold);

vold = -vold;

Pnew = flipud(Pold);
unew = flipud(uold);
vnew = flipud(vold);

% Restart iteration values
iterRestart = endIter+1;
totalIters  = endIter;
assignin('base','iterRestart',iterRestart);
assignin('base','totalIters',totalIters);

% Send final values back into base workspace
assignin('base','ErrorTot',ErrorTot);
assignin('base','Pnew',Pnew);
assignin('base','unew',unew);
assignin('base','vnew',vnew);
assignin('base','uold',uold);
assignin('base','vold',vold);

axes(handles.plotData);
cla;
hold on;
% Main vortex (Zone 1)
s = streamline(xx,yy,unew,vnew,0.88,0.1); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.02,0.02); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.7,0.7); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.97,0.9); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.97,0.8); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.97,0.7); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.97,0.6); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.97,0.5); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.97,0.4); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.97,0.3); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.97,0.2); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.97,0.1); set(s,'Color','black');
% Right vortex (Zone 2)
s = streamline(xx,yy,unew,vnew,0.97,0.05); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.93,0.05); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.94,0.05); set(s,'Color','black');
% Left vortex (Zone 3)
s = streamline(xx,yy,unew,vnew,0.03,0.03); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.04,0.06); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.02,0.06); set(s,'Color','black');
s = streamline(xx,yy,unew,vnew,0.02,0.05); set(s,'Color','black');
title(['Streamline: N = ',num2str(N),', Re = ',num2str(Re),...
            ', r = ',num2str(r),', eps_d_i_s_s = ',num2str(eps)])
xlabel('X Distance');
ylabel('Y Distance');
axis([0 1 0 1]);


% PUSH ------------------- Add a Streamline to Plot --------------------- %
function pushAddStreamlines_Callback(hObject, eventdata, handles)
xx = evalin('base','xx');
yy = evalin('base','yy');
unew = evalin('base','unew');
vnew = evalin('base','vnew');

[xSL,ySL] = ginput(10);

numSL = length(xSL);

axes(handles.plotData);
hold on;
for n = 1:1:numSL
    s = streamline(xx,yy,unew,vnew,xSL(n),ySL(n));
    set(s,'Color','black');
end

% PUSH ---------------- Clear all the Streamlines ----------------------- %
function pushClearStreamlines_Callback(hObject, eventdata, handles)
toPlotNum = evalin('base','toPlotNum');
xx        = evalin('base','xx');
yy        = evalin('base','yy');
unew      = evalin('base','unew');
vnew      = evalin('base','vnew');
N         = evalin('base','N');
Re        = evalin('base','Re');
r         = evalin('base','r');
eps       = evalin('base','eps');

if (toPlotNum == 1)
    axes(handles.plotData);
    cla;
elseif (toPlotNum == 6)
    axes(handles.plotData);
    cla;
    quiver(xx,yy,unew,vnew);
    title(['N = ',num2str(N),', Re = ',num2str(Re),...
                ', r = ',num2str(r),', eps_d_i_s_s = ',num2str(eps)])
    xlabel('X Grid Point');
    ylabel('Y Grid Point');
    axis([0 1 0 1]);
    grid on;
end
