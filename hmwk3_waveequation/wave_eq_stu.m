% CERI 7104/8104 Data Analysis in Geophysics
%
%
% solve 1D, unidirectional wave equation
%
% @author thgoebel, CERI U of Memphis
% -script modified from Eric Daub, http://www.ceri.memphis.edu/people/egdaub/teaching.html
clear
clf

%=======================1==================================================
%             TODO: set these parameters
%==========================================================================
% physical parameters: wave speed
c   = ??;

xmax = ??;
tmax = ??;

% boundary conditions
BC_l = ??;
BC_r = ??;
%%%   initial Gaussian pulse
A0     = ??;
loc    = ??;%mean
scale  = ??;%stdDev

% spatial increment
dx  = ??;

%----------------------additional params---------------------
% set up spatial grid
%nx  = xmax/(dx)+1;
a_x = ??; %spatial vector 

% time discretization
dt_dx = 0.5/c; % ratio between time step and spatial grid size
dt    = dt_dx*dx; % fixed value for dt/dx to ensure numerical stability
a_t   = ??;%time vector
%nt    = tmax/dt + 1;
%tmax  = (nt-1)*dt;

%=======================2==================================================
%                  solve PDE
%==========================================================================
% Initial conditions
a_u    = ??;

%tic; %tic toc keep track of computational time

% initialize amplitude matrix
m_u      = zeros(nx, nt);
a_unew = zeros(nx, 1);%amplitude profile at new time step

fig = figure;%create figure and axis out here for updateing plots
clf;
ax = axes;
for i_t=1:nt
    % boundary conditions
    a_u(1)   = ??;
    a_u(nx) = ??;
    for i_x = 2:nx % BD for spatial derivative
        a_unew(??) = a_u(??) - c*(dt_dx)*(a_u(??)-a_u(??));
    end
    a_u          = a_unew; %update amplitude profile
    m_u(:,i_t) = a_u;% save amplitudes at each time step
    %=======================3==============================================
    %                     plots
    %======================================================================
    %  plot amplitudes at each time step
    lines = plot(??, ??);
    xlabel(ax,'Distance' );
    ylabel(ax,'Amplitude' );
    % create title string with five elements
    plottitle = sprintf('%s%s%s%s%s','1D Wave Equation After ',num2str(i_t), 'time steps (t = ',num2str(dt*i_t), ')');
    title( plottitle );
    
    set(ax,'xlim',[0 1],'ylim',[0 A0]);
    set(ax,'xtick',[0 0.2 0.4 0.6 0.8 1],'ytick',[0:.3:A0]); 
    %saveas( fig,'hw2_wave_equation_%i.jpg');
    pause(0.1);
%     
end
%toc;
% F2 - plot amplitude profile at all time steps
fig2 = figure('Position', [1 1 900 600]);
clf;
%---plot last profile
ax2 = subplot(2,1,1);
plot(a_x, a_u);
xlabel(ax2, 'Distance' );
ylabel(ax2, 'Amplitude' );
hold on

%TODO: plot pcolor of amplitude profiles
ax3 = subplot(2,1,2);
??
??
??

axis on
xlabel(ax3,'Distance (km)' );
ylabel(ax3,'Time (s)' );
cbar = colorbar('eastoutside');
colormap(ax3, jet(256));
cbar.Label.String = 'Amplitude';
plottitle = sprintf('%s%s%s','Wave Speed: ',num2str(c), 'km/s');
title(ax3, plottitle );
saveas(fig2,'hw2_imshow_slowness.jpg');



	
