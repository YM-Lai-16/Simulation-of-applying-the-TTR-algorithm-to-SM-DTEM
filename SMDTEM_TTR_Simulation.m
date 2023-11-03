%%%% Matlab code for implementing TwIST-based tomographic reconstruction (TTR) 
%%%% algorithm of compressed ultrafast tomographic imaging (CUTI) to
%%%% the streak-mode dynamic transmission electron microscopy (SM-DTEM)
%%%% =====================================
%%%% This script requires the package of the TwIST algorithm available at 
%%%% http://www.lx.it.pt/~bioucas/TwIST/TwIST.htm
%%%% =====================================
%%%% The simulation and several customized functions are required. They are available at
%%%% https://github.com/YM-Lai-16/Simulation-of-applying-the-TTR-algorithm-to-SM-DTEM
%%%% The customized functions should be put into the same folder of "TwIST.m"

close all; clear all; clc

Norm = @(x)    (x-min(x(:)))./(max(x(:))-min(x(:))); % Create a normalization function

addpath(genpath('TwIST_v2')); % Load the folder includes the package of TwIST 

%% Simulate the SVCT spatiotemporal projections based on streak imaging
load('SMDTEM_Balls.mat');          % Load the simulated event
Event = Balls;    % Normalize the simulated event
[N_y0, N_x0, N_t] = size(Event);       % Get the size of the event

V_S1 = 1;                % S1: shear 1 pixel per frame along the downward direction
D_shearing_y_S1 = abs(round(V_S1.*N_t)); % Sheared distance at the sweeping speed S1

V_S2 = 0.5;             % S2: shear 0.5 pixel per frame along the downward direction
D_shearing_y_S2 = abs(round(V_S2.*N_t)); % Sheared distance at the sweeping speed S2

V_S3 = -1;             % S3: shear 1 pixel per frame along the upward direction
D_shearing_y_S3 = abs(round(V_S3.*N_t)); % Sheared distance at the sweeping speed S3

V_S4 = -0.5;              % S4: shear 0.5 pixel per frame along the upward direction
D_shearing_y_S4 = abs(round(V_S4.*N_t)); % Sheared distance at the sweeping speed S4

V_S5 = 1;               % S5: shear 1 pixel per frame along the rightward direction
D_shearing_x_S5 = abs(round(V_S5.*N_t)); % Sheared distance at the sweeping speed S5

V_S6 = 0.5;            % S6: shear 0.5 pixel per frame along the rightward direction
D_shearing_x_S6 = abs(round(V_S6.*N_t)); % Sheared distance at the sweeping speed S6

V_S7 = -1;              % S7: shear 1 pixel per frame along the leftward direction
D_shearing_x_S7 = abs(round(V_S7.*N_t)); % Sheared distance at the sweeping speed S7

V_S8 = -0.5;            % S8: shear 0.5 pixel per frame along the leftward direction
D_shearing_x_S8 = abs(round(V_S8.*N_t)); % Sheared distance at the sweeping speed S8

%%%% Pad zeros to the datacube. 
%%%% Downward & rightward: direction of padarray is 'post', which pads after the last array element 
%%%% Upward & leftward: direction of padarray is 'pre', which pads before the first array element 
%%%% If the shearing directions changes, assign the array and 'post' or 'pre' accordingly
for i = 1:N_t
    im_S1 = padarray(Event(:,:,i), [D_shearing_y_S1, 0], 0, 'post');  % Pad zeros to the bottom of the datacube
    im_cube_S1(:,:,i) = circshift(im_S1, [round(V_S1*(i-1)) 0]);    % Shear each frame downward at the sweeping speed S1
    
    im_S2 = padarray(Event(:,:,i), [D_shearing_y_S2, 0], 0, 'post');  % Pad zeros to the bottom of the datacube
    im_cube_S2(:,:,i) = circshift(im_S2, [round(V_S2*(i-1)) 0]);    % Shear each frame downward at the sweeping speed S2
    
    im_S3 = padarray(Event(:,:,i), [D_shearing_y_S3, 0], 0, 'pre');  % Pad zeros to the top of the datacube
    im_cube_S3(:,:,i) = circshift(im_S3, [round(V_S3*(i-1)) 0]);    % Shear each frame downward at the sweeping speed S3
    
    im_S4 = padarray(Event(:,:,i), [D_shearing_y_S4, 0], 0, 'pre');  % Pad zeros to the top of the datacube
    im_cube_S4(:,:,i) = circshift(im_S4, [round(V_S4*(i-1)) 0]);    % Shear each frame upward at the sweeping speed S4
        
    im_S5 = padarray(Event(:,:,i), [0, D_shearing_x_S5], 0, 'post');  % Pad zeros to the right of the datacube
    im_cube_S5(:,:,i) = circshift(im_S5, [0 round(V_S5*(i-1))]);    % Shear each frame upward at the sweeping speed S5
    
    im_S6 = padarray(Event(:,:,i), [0, D_shearing_x_S6], 0, 'post');  % Pad zeros to the right of the datacube
    im_cube_S6(:,:,i) = circshift(im_S6, [0 round(V_S6*(i-1))]);    % Shear each frame upward at the sweeping speed S6
    
    im_S7 = padarray(Event(:,:,i), [0, D_shearing_x_S7], 0, 'pre');  % Pad zeros to the left of the datacube
    im_cube_S7(:,:,i) = circshift(im_S7, [0 round(V_S7*(i-1))]);    % Shear each frame upward at the sweeping speed S7
    
    im_S8 = padarray(Event(:,:,i), [0, D_shearing_x_S8], 0, 'pre');  % Pad zeros to the left of the datacube
    im_cube_S8(:,:,i) = circshift(im_S8, [0 round(V_S8*(i-1))]);    % Shear each frame upward at the sweeping speed S8

end

E_0 = Norm(sum(Event, 3));              % Spatiotemporal integration of the unsheared datacube
E_1 = Norm(sum(im_cube_S1, 3));               % Spatiotemporal integration of datacube sheared at the sweeping speed S1
E_2 = Norm(sum(im_cube_S2, 3));               % Spatiotemporal integration of datacube sheared at the sweeping speed S2
E_3 = Norm(sum(im_cube_S3, 3));               % Spatiotemporal integration of datacube sheared at the sweeping speed S3
E_4 = Norm(sum(im_cube_S4, 3));               % Spatiotemporal integration of datacube sheared at the sweeping speed S4
E_5 = Norm(sum(im_cube_S5, 3));               % Spatiotemporal integration of datacube sheared at the sweeping speed S5
E_6 = Norm(sum(im_cube_S6, 3));               % Spatiotemporal integration of datacube sheared at the sweeping speed S6
E_7 = Norm(sum(im_cube_S7, 3));               % Spatiotemporal integration of datacube sheared at the sweeping speed S7
E_8 = Norm(sum(im_cube_S8, 3));               % Spatiotemporal integration of datacube sheared at the sweeping speed S8

I_frame = ones(N_y0, N_x0);
for i = 1:N_t
    Full_vector(:,i) = I_frame(:);        % Generate an all one matrix the same to the event
end

%% Creat the forward model of E_1 
Size_y_S1 = N_y0 + D_shearing_y_S1; % Size in the y-direction after temporal shearing
Size_x_S1 = N_x0; % Size in the x-direction after temporal shearing

for i = 1:N_t
    if V_S1 >= 0
        shift_y_S1(i) = round((i-1).*V_S1); % The shearing distance of each frame along the downward direction
    else
        shift_y_S1(i) = D_shearing_y_S1+round((i-1).*V_S1); % The shearing distance of each frame along the upward direction
    end
    shift_x_S1(i) = 0; % The shearing distance of each frame in the x-direction
end

count = 1;
for i = 1:N_t        % t  
    for k = 1:N_x0  % x
        for j = 1:N_y0  % y
        y_coordinate(count) = j + (k-1)*Size_y_S1 + shift_y_S1(i)+ ...
            Size_y_S1*(shift_x_S1(i));  % Assign the y-coordinate according to the sensing operator
        x_coordinate(count) = (i-1)*N_y0*N_x0 + (k-1)*N_y0 + j; 
        t_coordinate(count) = Full_vector(j + (k-1)*N_y0, i);
        count = count + 1;
        end
    end
end
S_1=sparse(y_coordinate,x_coordinate,t_coordinate,...
    Size_x_S1*Size_y_S1,N_y0*N_x0*N_t); % Generate the sensing operator at the shearing speed of S1
y_verify_S1 = Norm(reshape(S_1*Event(:),Size_y_S1,Size_x_S1)); % Create a measurement to verify the constructed forward model
figure, imagesc(abs(y_verify_S1 - E_1)); axis equal;axis off;colormap parula, title('Verify S1'); % All the elements should be canceled to be zero

clear x_coordinate;clear y_coordinate;clear t_coordinate; 

%% Creat the forward model of E_2
Size_y_S2 = N_y0 + D_shearing_y_S2; % Size of y after temporal shearing
Size_x_S2 = N_x0; % Size of x after temporal shearing

for i = 1:N_t
    if V_S2 >= 0
        shift_y_S2(i) = round((i-1).*V_S2); % The shearing distance of each frame along the downward direction
    else
        shift_y_S2(i) = D_shearing_y_S2+round((i-1).*V_S2); % The shearing distance of each frame along the upward direction
    end
    shift_x_S2(i) = 0; % The shearing distance of each frame in x-direction
end

count = 1;
for i = 1:N_t        % t  
    for k = 1:N_x0  % x
        for j = 1:N_y0  % y
        y_coordinate(count) = j + (k-1)*Size_y_S2 + shift_y_S2(i)+ ...
            Size_y_S2*(shift_x_S2(i));  % Assign the y-coordinate according to the sensing operator
        x_coordinate(count) = (i-1)*N_y0*N_x0 + (k-1)*N_y0 + j; 
        t_coordinate(count) = Full_vector(j + (k-1)*N_y0, i);
        count = count + 1;
        end
    end
end
S_2=sparse(y_coordinate,x_coordinate,t_coordinate,...
    Size_x_S2*Size_y_S2,N_y0*N_x0*N_t); % Generate the sensing operator at the shearing speed of S2
y_verify_S2 = Norm(reshape(S_2*Event(:),Size_y_S2,Size_x_S2)); % Create a measurement to verify the constructed forward model
figure, imagesc(abs(y_verify_S2 - E_2)); axis equal;axis off;colormap parula, title('Verify S2'); % All the elements should be canceled to be zero

clear x_coordinate;clear y_coordinate;clear t_coordinate; 

%% Creat the forward model of E_3
Size_y_S3 = N_y0 + D_shearing_y_S3; % Size of y after temporal shearing
Size_x_S3 = N_x0; % Size of x after temporal shearing

for i = 1:N_t
    if V_S3 >= 0
        shift_y_S3(i) = round((i-1).*V_S3); % The shearing distance of each frame along the downward direction
    else
        shift_y_S3(i) = D_shearing_y_S3+round((i-1).*V_S3); % The shearing distance of each frame along the upward direction
    end
    shift_x_S3(i) = 0; % The shearing distance of each frame in the x-direction
end

count = 1;
for i = 1:N_t        % t
    for k = 1:N_x0  % x
        for j = 1:N_y0  % y
        y_coordinate(count) = j + (k-1)*Size_y_S3 + shift_y_S3(i)+ ...
            Size_y_S3*(shift_x_S3(i));  % Assign the y-coordinate according to the sensing operator
        x_coordinate(count) = (i-1)*N_y0*N_x0 + (k-1)*N_y0 + j;
        t_coordinate(count) = Full_vector(j + (k-1)*N_y0, i);
        count = count + 1;
        end
    end
end
S_3=sparse(y_coordinate,x_coordinate,t_coordinate,...
    Size_x_S3*Size_y_S3,N_y0*N_x0*N_t); % Generate the sensing operator at the shearing speed of S3
y_verify_S3 = Norm(reshape(S_3*Event(:),Size_y_S3,Size_x_S3)); % Create a measurement to verify the constructed forward model
figure, imagesc(abs(y_verify_S3 - E_3)); axis equal;axis off;colormap parula, title('Verify S3');% All the elements should be canceled to be zero

clear x_coordinate;clear y_coordinate;clear t_coordinate;

%% Creat the forward model of E_4
Size_y_S4 = N_y0 + D_shearing_y_S4; % Size of y after temporal shearing
Size_x_S4 = N_x0; % Size of x after temporal shearing

for i = 1:N_t
    if V_S4 >= 0
        shift_y_S4(i) = round((i-1).*V_S4); % The shearing distance of each frame along the downward direction
    else
        shift_y_S4(i) = D_shearing_y_S4+round((i-1).*V_S4); % The shearing distance of each frame along the upward direction
    end
    shift_x_S4(i) = 0; % The shearing distance of each frame in the x-direction
end

count = 1;
for i = 1:N_t        % t
    for k = 1:N_x0  % x
        for j = 1:N_y0  % y
        y_coordinate(count) = j + (k-1)*Size_y_S4 + shift_y_S4(i)+ ...
            Size_y_S4*(shift_x_S4(i));  % Assign the y-coordinate according to the sensing operator
        x_coordinate(count) = (i-1)*N_y0*N_x0 + (k-1)*N_y0 + j;
        t_coordinate(count) = Full_vector(j + (k-1)*N_y0, i);
        count = count + 1;
        end
    end
end
S_4=sparse(y_coordinate,x_coordinate,t_coordinate,...
    Size_x_S4*Size_y_S4,N_y0*N_x0*N_t); % Generate the sensing operator at the shearing speed of S4
y_verify_S4 = Norm(reshape(S_4*Event(:),Size_y_S4,Size_x_S4)); % Create a measurement to verify the constructed forward model
figure, imagesc(abs(y_verify_S4 - E_4)); axis equal;axis off;colormap parula, title('Verify S4'); % All the elements should be canceled to be zero

clear x_coordinate;clear y_coordinate;clear t_coordinate;

%% Creat the forward model of E_5
Size_y_S5 = N_y0; % Size of y after temporal shearing
Size_x_S5 = N_x0 + D_shearing_x_S5; % Size of x after temporal shearing

for i = 1:N_t
    if V_S5 >= 0
        shift_x_S5(i) = round((i-1).*V_S5); % The shearing distance of each frame along the rightward direction
    else
        shift_x_S5(i) = D_shearing_x_S5+round((i-1).*V_S5); % The shearing distance of each frame along the leftward direction
    end
    shift_y_S5(i) = 0; % The shearing distance of each frame in the y-direction
end

count = 1;
for i = 1:N_t        % t
    for k = 1:N_x0  % x
        for j = 1:N_y0  % y
        y_coordinate(count) = j + (k-1)*Size_y_S5 + shift_y_S5(i)+ ...
            Size_y_S5*(shift_x_S5(i));  % Assign the y-coordinate according to the sensing operator
        x_coordinate(count) = (i-1)*N_y0*N_x0 + (k-1)*N_y0 + j;
        t_coordinate(count) = Full_vector(j + (k-1)*N_y0, i);
        count = count + 1;
        end
    end
end
S_5=sparse(y_coordinate,x_coordinate,t_coordinate,...
    Size_x_S5*Size_y_S5,N_y0*N_x0*N_t); % Generate the sensing operator at the shearing speed of S5
y_verify_S5 = Norm(reshape(S_5*Event(:),Size_y_S5,Size_x_S5)); % Create a measurement to verify the constructed forward model
figure, imagesc(abs(y_verify_S5 - E_5)); axis equal;axis off;colormap parula, title('Verify S5'); % All the elements should be canceled to be zero

clear x_coordinate;clear y_coordinate;clear t_coordinate;

%% Creat the forward model of E_6
Size_y_S6 = N_y0; % Size of y after temporal shearing
Size_x_S6 = N_x0 + D_shearing_x_S6; % Size of x after temporal shearing

for i = 1:N_t
    if V_S6 >= 0
        shift_x_S6(i) = round((i-1).*V_S6); % The shearing distance of each frame along the rightward direction
    else
        shift_x_S6(i) = D_shearing_x_S6+round((i-1).*V_S6); % The shearing distance of each frame along the leftward direction
    end
    shift_y_S6(i) = 0; % The shearing distance of each frame in the y-direction
end

count = 1;
for i = 1:N_t        % t
    for k = 1:N_x0  % x
        for j = 1:N_y0  % y
        y_coordinate(count) = j + (k-1)*Size_y_S6 + shift_y_S6(i)+ ...
            Size_y_S6*(shift_x_S6(i));  % Assign the y-coordinate according to the sensing operator
        x_coordinate(count) = (i-1)*N_y0*N_x0 + (k-1)*N_y0 + j;
        t_coordinate(count) = Full_vector(j + (k-1)*N_y0, i);
        count = count + 1;
        end
    end
end
S_6=sparse(y_coordinate,x_coordinate,t_coordinate,...
    Size_x_S6*Size_y_S6,N_y0*N_x0*N_t); % Generate the sensing operator at the shearing speed of S6
y_verify_S6 = Norm(reshape(S_6*Event(:),Size_y_S6,Size_x_S6)); % Create a measurement to verify the constructed forward model
figure, imagesc(abs(y_verify_S6 - E_6)); axis equal;axis off;colormap parula, title('Verify S6'); % All the elements should be canceled to be zero

clear x_coordinate;clear y_coordinate;clear t_coordinate;

%% Creat the forward model of E_7
Size_y_S7 = N_y0; % Size of y after temporal shearing
Size_x_S7 = N_x0 + D_shearing_x_S7; % Size of x after temporal shearing

for i = 1:N_t
    if V_S7 >= 0
        shift_x_S7(i) = round((i-1).*V_S7); % The shearing distance of each frame along the rightward direction
    else
        shift_x_S7(i) = D_shearing_x_S7+round((i-1).*V_S7); % The shearing distance of each frame along the leftward direction
    end
    shift_y_S7(i) = 0; % The shearing distance of each frame in the y-direction
end

count = 1;
for i = 1:N_t        % t
    for k = 1:N_x0  % x
        for j = 1:N_y0  % y
        y_coordinate(count) = j + (k-1)*Size_y_S7 + shift_y_S7(i)+ ...
            Size_y_S7*(shift_x_S7(i));  % Assign the y-coordinate according to the sensing operator
        x_coordinate(count) = (i-1)*N_y0*N_x0 + (k-1)*N_y0 + j;
        t_coordinate(count) = Full_vector(j + (k-1)*N_y0, i);
        count = count + 1;
        end
    end
end
S_7=sparse(y_coordinate,x_coordinate,t_coordinate,...
    Size_x_S7*Size_y_S7,N_y0*N_x0*N_t); % Generate the sensing operator at the shearing speed of S7
y_verify_S7 = Norm(reshape(S_7*Event(:),Size_y_S7,Size_x_S7)); % Create a measurement to verify the constructed forward model
figure, imagesc(abs(y_verify_S7 - E_7)); axis equal;axis off;colormap parula, title('Verify S7'); % All the elements should be canceled to be zero

clear x_coordinate;clear y_coordinate;clear t_coordinate;

%% Creat the forward model of E_8
Size_y_S8 = N_y0; % Size of y after temporal shearing
Size_x_S8 = N_x0 + D_shearing_x_S8; % Size of x after temporal shearing

for i = 1:N_t
    if V_S8 >= 0
        shift_x_S8(i) = round((i-1).*V_S8); % The shearing distance of each frame along the rightward direction
    else
        shift_x_S8(i) = D_shearing_x_S8+round((i-1).*V_S8); % The shearing distance of each frame along the leftward direction
    end
    shift_y_S8(i) = 0; % The shearing distance of each frame in the y-direction
end

count = 1;
for i = 1:N_t        % t
    for k = 1:N_x0  % x
        for j = 1:N_y0  % y
        y_coordinate(count) = j + (k-1)*Size_y_S8 + shift_y_S8(i)+ ...
            Size_y_S8*(shift_x_S8(i));  % Assign the y-coordinate according to the sensing operator
        x_coordinate(count) = (i-1)*N_y0*N_x0 + (k-1)*N_y0 + j;
        t_coordinate(count) = Full_vector(j + (k-1)*N_y0, i);
        count = count + 1;
        end
    end
end
S_8=sparse(y_coordinate,x_coordinate,t_coordinate,...
    Size_x_S8*Size_y_S8,N_y0*N_x0*N_t); % Generate the sensing operator at the shearing speed of S8
y_verify_S8 = Norm(reshape(S_8*Event(:),Size_y_S8,Size_x_S8)); % Create a measurement to verify the constructed forward model
figure, imagesc(abs(y_verify_S8 - E_8)); axis equal;axis off;colormap parula, title('Verify S8'); % All the elements should be canceled to be zero

clear x_coordinate;clear y_coordinate;clear t_coordinate;

%% Creat the forward model of the static image (E_0)
count = 1;
for i = 1:N_t        
    for k = 1:N_x0  
        for j = 1:N_y0  
        y_coordinate(count) = j + (k-1)*N_y0 + 0;
        x_coordinate(count) = (i-1)*N_y0*N_x0 + (k-1)*N_y0 + j;
        t_coordinate(count) = Full_vector(j + (k-1)*N_y0,i);
        count = count + 1;
        end
    end
end

S_0=sparse(y_coordinate,x_coordinate,t_coordinate,N_x0*N_y0,N_y0*N_x0*N_t);    % Generate the sensing operator of the static image
y_verify_static = Norm(reshape(S_0*Event(:),N_y0,N_x0)); % Create a measurement to verify the constructed forward model
figure, imagesc(abs(y_verify_static - E_0)); axis equal;axis off;colormap parula, title('Static'); % All the elements should be canceled to be zero

clear x_coordinate;clear y_coordinate;clear t_coordinate;

%% Create the sets of streak measurements (spatiotemporal projections) and sensing operators
E = [E_0(:);...
    E_1(:);...
    E_2(:);...
    E_3(:);...
    E_4(:);...
    E_5(:);...
    E_6(:);...
    E_7(:);...
    E_8(:);...
        ];   
S = [S_0;...
    S_1;...
    S_2;...
    S_3;...
    S_4;...
    S_5;...
    S_6;...
    S_7;...
    S_8;...
        ];     

%% Run TwIST
%%%% Note: the default maximum number of iterations is 1000. It can be set inside the function "TwIST.m"
close all;

%%%% TwIST parameters
initial_guess = (S'*E);
lambda = 1e-4;  
tau = 0.6;  % Adjust the regularizer parameter 
tolA = 1e-15;
tv_iters=2;

Psi = @(x,th)  TTR_denoise_func(x,th,N_y0,N_x0,N_t,tv_iters); % Define the regularizer
Phi = @(x) TTR_TVphi(x,N_y0,N_x0,N_t);

[Im_cube,x_debias_twist,obj_twist,...
    times_twist,debias_start_twist,mse]= ...
         TwIST(E,S,tau, ...
         'Lambda', lambda, ...
         'Debias',0,...
         'Monotone',1,...
         'Sparse', 1,...
         'Initialization',initial_guess(:),...
         'Psi',Psi,...
         'Phi',Phi,...
         'StopCriterion',1,...
         'ToleranceA',tolA,...
         'Verbose', 1);

Im_cube = Norm(reshape(Im_cube,N_y0,N_x0,N_t)); % Output from the TwIST algorithm

%% Display the result
bg = mean(Im_cube(end-10:end, 1:10, end), 'all'); % Calculate the background of the result

figure(1), set(gcf,'color','k');
for i=1:1:size(Im_cube, 3)
        im = Im_cube(:,:,i);
        im(im<=bg) = bg;
        im = Norm(im);  % Background subtraction of the result
        
        im_gt = Event(:,:,i); % Read the ground truth

        subplot(1,2,1), imagesc(im_gt,[0 1]); axis off;axis image;colormap('gray');title('Ground truth', 'Color', 'white', 'FontSize', 16);   %Display the ground truth 
        subplot(1,2,2), imagesc(im,[0 1]); axis off;axis image;colormap('gray'); title('Reconstruction', 'Color', 'white', 'FontSize', 16);  %Display the reconstructed video 
        
        pause(0.1);
end

