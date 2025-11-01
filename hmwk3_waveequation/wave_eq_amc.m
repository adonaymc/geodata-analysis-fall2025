
% CERI 7104/8104 Data Analysis in Geophysics
%
%
% solve 1D, unidirectional wave equation
%
% @author thgoebel, CERI U of Memphis
% -script modified from Eric Daub, http://www.ceri.memphis.edu/people/egdaub/teaching.html
clear;
close all;

%=======================1==================================================
%             TODO: set these parameters
%==========================================================================
% physical parameters: wave speed
c   = 0.5; %km/s

xmax = 1;
tmax = 1;

% boundary conditions
BC_l = 0;
BC_r = 0;
%%%   initial Gaussian pulse
A0     = 1.5;
loc    = 0.3; %mean
scale  = 0.005; %stdDev
% spatial increment
dx  = 0.01; %km

%----------------------additional params---------------------
% set up spatial grid
nx  = xmax/dx + 1;
a_x = linspace(0, xmax, nx); %spatial vector 

% time discretization
dt_dx = 0.5/c; % ratio between time step and spatial grid size
dt    = dt_dx*dx; % fixed value for dt/dx to ensure numerical stability
a_t   = 0:dt:tmax; %time vector
nt    = length(a_t); %number of time steps

%=======================2==================================================
%                  solve PDE
%==========================================================================
% Initial conditions
a_u    = A0*exp(-(a_x-loc).^2/(2*scale^2)); %initial amplitude profile  

tic; %tic toc keep track of computational time

% initialize amplitude matrix
m_u      = zeros(nx, nt);
a_unew = zeros(nx, 1);  %amplitude profile at new time step

% visualize the initial contitions in a plot before starting time loop
figure;
plot(a_x, a_u, 'LineWidth', 2);
xlabel('Distance' );
ylabel('Amplitude' );
title('1D Wave Equation Initial Conditions');
set(gca,'xlim',[0 1],'ylim',[0 A0]);
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1],'ytick',[0:.3:A0]); 
%saveas( fig,'hw2_wave_equation_initial_conditions.jpg');
pause(1);

fig = figure; %create figure and axis out here for updating plots
clf;
ax = axes;
for i_t=1:nt
    % boundary conditions
    a_u(1)   = BC_l;
    a_u(nx) = BC_r;
    for i_x = 2:nx % BD for spatial derivative
        a_unew(i_x) = a_u(i_x) - c*(dt_dx)*(a_u(i_x)-a_u(i_x-1));
    end
    a_u          = a_unew; %update amplitude profile
    m_u(:,i_t) = a_u;  % save amplitudes at each time step
    %=======================3==============================================
    %                     plots
    %======================================================================
    %  plot amplitudes at each time step
    lines = plot(ax, a_x, a_u, 'LineWidth', 2);
    % set plot properties
    xlabel(ax,'Distance' );
    ylabel(ax,'Amplitude' );
    % create title string with five elements
    plottitle = sprintf('%s%s%s%s%s','1D Wave Equation After ',num2str(i_t), ' time steps (t = ',num2str(dt*i_t), ')');
    title( plottitle );
    
    set(ax,'xlim',[0 1],'ylim',[0 A0]);
    set(ax,'xtick',[0 0.2 0.4 0.6 0.8 1],'ytick',[0:.3:A0]); 
    %saveas( fig,'hw2_wave_equation_%i.jpg');
    drawnow;
    frame = getframe(fig);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if i_t == 1
        imwrite(imind,cm,'hw2_wave_equation.gif','gif', 'Loopcount',inf,'DelayTime',0.1);
    else
        imwrite(imind,cm,'hw2_wave_equation.gif','gif','WriteMode','append','DelayTime',0.1);
    end
    % to slow down the animation, uncomment the next line
    % pause(0.1);
end

toc;

% F2 - plot amplitude profile at all time steps
fig2 = figure('Position', [1 1 900 600]);
clf;
%---plot first, middle, last time step
ax1 = subplot(4,1,1);
plot2_1 = plot(ax1, a_x, m_u(:,1), 'LineWidth', 2);
xlabel(ax1,'Distance' );
ylabel(ax1,'Amplitude' );
title(ax1, 'Amplitude Profile at First Time Step');
set(ax1,'xlim',[0 1],'ylim',[0 A0]);
set(ax1,'xtick',[0 0.2 0.4 0.6 0.8 1],'ytick',[0:.3:A0]);
%---middle time step
ax2 = subplot(4,1,2);
mid_t = round(nt/2);
plot2_2 = plot(ax2, a_x, m_u(:,mid_t), 'LineWidth', 2);
xlabel(ax2,'Distance' );
ylabel(ax2,'Amplitude' );
title(ax2, sprintf('%s%s%s','Amplitude Profile at Middle Time Step (t = ',num2str(dt*mid_t),')'));
set(ax2,'xlim',[0 1],'ylim',[0 A0]);
set(ax2,'xtick',[0 0.2 0.4 0.6 0.8 1],'ytick',[0:.3:A0]);

%---last time step
ax3 = subplot(4,1,3);
plot2_3 = plot(ax3, a_x, m_u(:,nt), 'LineWidth', 2);
xlabel(ax3,'Distance' );
ylabel(ax3,'Amplitude' );
title(ax3, sprintf('%s%s%s','Amplitude Profile at Last Time Step (t = ',num2str(dt*nt),')'));
set(ax3,'xlim',[0 1],'ylim',[0 A0]);
set(ax3,'xtick',[0 0.2 0.4 0.6 0.8 1],'ytick',[0:.3:A0]);

%---all time steps
ax4 = subplot(4,1,4);
hold on;
for i_t=1:nt
    plot2_4 = plot(ax4, a_x, m_u(:,i_t), 'LineWidth', 1, 'Color', [0.5 0.5 0.5]);
end
xlabel(ax4,'Distance' );
ylabel(ax4,'Amplitude' );
title(ax4, 'Amplitude Profile at All Time Steps');
set(ax4,'xlim',[0 1],'ylim',[0 A0]);
set(ax4,'xtick',[0 0.2 0.4 0.6 0.8 1],'ytick',[0:.3:A0]);



%=======================4==================================================

% %TODO: plot pcolor of amplitude profiles
fig3 = figure('Position', [1 1 900 600]);
clf;
[m_x, m_t] = meshgrid( a_x, a_t);
plot3 = pcolor(m_x, m_t, flipud(m_u'));
shading interp; 
axis on
xlabel('Distance (km)' );
ylabel('Time (s)' );
cbar = colorbar('eastoutside');
colormap(jet(256));
cbar.Label.String = 'Amplitude';
plottitle = sprintf('%s%s%s','Wave Speed: ',num2str(c), 'km/s');
title(plottitle);