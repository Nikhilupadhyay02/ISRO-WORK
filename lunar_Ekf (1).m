close all;
clear all;
clc;

  file = 'C:\Users\HP\Desktop\Aman_change\Aman\Kalman_Aman\SATB_L1.csv';
% si=input('Enter starting index number: ');
% ei=input('Enter ending index number: ');
%  file = 'C:\Users\HP\Desktop\Nikhil Upadhyay\Matlab\GroundVehicle1_ECF_position_Velocity.xlsx';
%  read_file=[];
 
%  allData=[];
%  file_length=length(file);

  %read_file = xlsread(file);
  
%   allData = xlsread(file);
allData = xlsread(file);
  
  Start_idx = 1;%input('Enter the starting index for satellite data: ');
  end_idx = 3000;%input('Enter the ending index for satellite data: ');
  
% allData = xlsread(file); 
% sum=0;
  pr_indexes = [18,40,62,84,128,150];
  Ino_indexes=[22,44,66,88,132,154];
  Tropo_indexes=[23,45,67,89,133,155];
  p_x_indexes=[25,47,69,91,135,157];
  p_y_indexes=[26,48,70,92,136,158];
  p_z_indexes=[27,49,71,93,137,159];
  clock_correction_indexes=[32,54,76,98,142,164];

 x_refrence = 1146835.138;
 y_refrence = 5393290.448;
 z_refrence = 3196304.915;
 
%  num_iterations=10;
 %intialize arrays to store Rss errors for each row
 rss_errors = zeros(end_idx-Start_idx+1,1);
 rss_errors_kalman = zeros(end_idx-Start_idx+1,1);
 
 %iterate over each row in the specified range 
 
 for rowIdx=Start_idx:end_idx
     %extract data for the current row
     pr_data=allData(rowIdx,pr_indexes);
     Ino_data=allData(rowIdx,Ino_indexes);
     Tropo_data=allData(rowIdx,Tropo_indexes);
     p_x_data= allData(rowIdx,p_x_indexes);
     p_y_data= allData(rowIdx,p_y_indexes);
     p_z_data= allData(rowIdx,p_z_indexes);
     clock_correction_data = allData(rowIdx,clock_correction_indexes);
     
     Final_Pseudorange =  pr_data-Ino_data- Tropo_data+clock_correction_data;
     
     x_initial=p_x_data;
     y_initial=p_y_data;
     z_initial=p_z_data;
     n=length(x_initial);
     
 XG1=0;
 XY1=0;
 XZ1=0;
 XBU=0;

% %  rss_error = sqrt((x_refrence - XG1)^2+ (y_refrence - XY1)^2 + (z_refrence - XZ1)^2);
% %  disp(rss_error)
% %  end                 
% num_iterations=10;
 for t=1:10
     for i=1:n
     X_initial_diff=-(p_x_data-XG1);
     Y_initial_diff=-(p_y_data-XY1);
     Z_initial_diff=-(p_z_data-XZ1);
     
     distancegu=sqrt(X_initial_diff.^2+Y_initial_diff.^2+Z_initial_diff.^2)+XBU;
     del_rho=-distancegu+Final_Pseudorange;
     inv_distance=1 ./(distancegu-XBU);
     
     %compute gradients
     ax1= X_initial_diff .*inv_distance;
     ay1= Y_initial_diff .*inv_distance;
     az1= Z_initial_diff .*inv_distance;
     
     
   %create matrix for least squares solution
%    h1_matrix=[ax1',ay1',az1',ones(size(del_rho))];
     end
 h1_matrix=[ax1',ay1',az1',(ones(6,1))];
  H_non_kalman=h1_matrix;

   delta_X=(pinv(h1_matrix))*del_rho';
   
   XG1=XG1+delta_X(1);
   XY1=XY1+delta_X(2);
   XZ1=XZ1+delta_X(3);
   XBU=XBU+delta_X(4);
   
  
   
 end
 
  if rowIdx==1
       x_kal_ig = XG1;
       y_kal_ig = XY1;
       z_kal_ig = XZ1;
        H_non_kalman=h1_matrix;     
   end
 %compute the final RSS error for the current row
 rss_error_no_ekf=sqrt((x_refrence-XG1)^2+(y_refrence-XY1)^2+(z_refrence-XZ1)^2);
%  disp(['Rss Error:',num2str(rss_error)]);
%store the RSS error in the array
rss_errors(rowIdx-Start_idx+1)=rss_error_no_ekf;
 
%    disp('RSS Errors for each row:');
%    disp(rss_errors);
 end
   %% kalman
 x=[x_kal_ig;y_kal_ig;z_kal_ig;0;0];
 
 
  sampling_freq=1;
sampling_time=1/sampling_freq;
h0=6.32e-22;
h_1=6.32e-20;
h_2=2.00e-34;
T1=1;
T2=10;
A1 =1/((T1*power(T2,3))-((T2*power(T1,3))));
A2 = [power(T2,3), -power(T1,3);-3*T2, 3*T1];
A3 = [T1, power(T1,2), power(T1,3);T2, power(T2,2), power(T2,3)];
A4 = [h0/2; 2*h_1;(2/3)*(pi)*h_2];
output = A1.*(A2*A3*A4);
S_phi = output(1,1);
S_freq = output(2,1);
C=299792458;
%making of the q matrix
Q1 = C^2.*[(S_phi*sampling_freq)+((power(sampling_freq,3)/3)*S_freq), (power(sampling_freq,2)/2)*S_freq;(power(sampling_freq,2)/2)*S_freq, sampling_freq*S_freq];

% F=eye(5,5);
% F(4,5)=sampling_time;
Q=zeros(5);
Q(4:5,4:5)=Q1;
 
 H_non_kalman = [H_non_kalman,zeros(6,1)];
 Pa=pinv(H_non_kalman'*H_non_kalman);
 
% Pa=pinv(F'*F);
P_est= [Pa(1,1), Pa(1,2), Pa(1,3), Pa(1,4),0;
    Pa(2,1), Pa(2,2), Pa(2,3), Pa(2,4), 0;
    Pa(3,1), Pa(3,2), Pa(3,3), Pa(3,4), 0;
    Pa(4,1), Pa(4,2), Pa(4,3), Pa(4,4), 0;
    0,0,0, 0, 2*Q1(2,2)];

     n=length(x_initial);

R=eye(n)*3*5.22*5.22;

F=eye(5,5);
F(4,5)=sampling_time;


for rowIdx=Start_idx:end_idx
    pr_data=allData(rowIdx,pr_indexes);
     Ino_data=allData(rowIdx,Ino_indexes);
     Tropo_data=allData(rowIdx,Tropo_indexes);
     p_x_data= allData(rowIdx,p_x_indexes);
     p_y_data= allData(rowIdx,p_y_indexes);
     p_z_data= allData(rowIdx,p_z_indexes);
     clock_correction_data = allData(rowIdx,clock_correction_indexes);
     
     Final_Pseudorange =  pr_data-Ino_data- Tropo_data+clock_correction_data;
     
     x_initial=p_x_data;
     y_initial=p_y_data;
     z_initial=p_z_data;
%      n= length(x_initial);
     
%    plot(rss_errors);


%    for i=1:n
%    H_matrix=eye(4);
%     x_est=F*x_est;
%     P_est=F*P_est;

         
     x_initial_diff=-(p_x_data-x(1));
     y_initial_diff=-(p_y_data-x(2));
     z_initial_diff=-(p_z_data-x(3));
    

    distancegu=sqrt(x_initial_diff.^2+ y_initial_diff.^2+z_initial_diff.^2)+x(4);
    
    del_rho=-distancegu+Final_Pseudorange;

    
    
     inv_distance=1./(distancegu-x(4));
     
     %calculate partial derivtive
     
     ax1_k=x_initial_diff .* inv_distance;
     ay1_k=y_initial_diff .* inv_distance;
     az1_k=z_initial_diff .* inv_distance;
     
%      H(i,:)=[x_initial_diff(i)*inv_distance,y_initial_diff(i)*inv_distance,z_initial_diff(i)*inv_distance,-1];
%     end

%     rss_error = sqrt((x_refrence - x(1))^2+ (y_refrence - x(2))^2 + (z_refrence - x(3))^2);
  
% h1_matrix= [ax1',ay1',az1',ones(6,1)];
    H=[ax1_k',ay1_k',az1_k',ones(n,1),zeros(6,1)];
% H=[H,ones(6,1),zeros(6,1)];
   H_matrix =H;
   P_est=F*P_est*F'+Q;
    S=pinv((H_matrix*P_est*H_matrix')+R);
    K=(P_est*H_matrix')*S;
   
           x=x+(K*del_rho');

    P_est=(eye(5)-K*H_matrix)*P_est;

%        
%        Fact=pinv((H_matrix*P*H_matrix')+R);
%        K=(P*H_matrix')*Fact;
      

   

    

 rss_error_ekf = sqrt((x_refrence - x(1))^2+ (y_refrence - x(2))^2 + (z_refrence - x(3))^2);
 
 rss_errors_kalman(rowIdx-Start_idx+1)=rss_error_ekf;
 
%  disp(rss_errors_kalman);
%  disp(rss_errors);
 
%  rss_errors=rand(1,100);
%  rss_errors_kalman=rand(1,100);
end
% low dynamic
  file = 'C:\Users\HP\Desktop\Nikhil Upadhyay\Matlab\GroundVehicle1_ECF_position_Velocity.xlsx';

 
allData = xlsread(file);
 vx_indexes=[4,5,6];
 v_xu_guess = -0.015635633292553;
 v_yu_guess = -0.021385502742435;
 v_zu_guess = -0.024359600504399;
 
  x=[x_kal_ig;y_kal_ig;z_kal_ig;v_xu_guess;v_yu_guess;v_zu_guess;1;0];
% x=[0,0,0,0,0,0,0,0];
  x_inv=x';
 
sampling_freq=1;
sampling_time=1/sampling_freq;
h0=6.32e-22;
h_1=6.32e-20;
h_2=2.00e-34;
T1=1;
T2=10;
A1 =1/((T1*power(T2,3))-((T2*power(T1,3))));
A2 = [power(T2,3), -power(T1,3);-3*T2, 3*T1];
A3 = [T1, power(T1,2), power(T1,3);T2, power(T2,2), power(T2,3)];
A4 = [h0/2; 2*h_1;(2/3)*(pi)*h_2];
output = A1.*(A2*A3*A4);
S_phi = output(1,1);
S_freq = output(2,1);
C=299792458;
%making of the q matrix
Q1 = C^2.*[(S_phi*sampling_freq)+((power(sampling_freq,3)/3)*S_freq), (power(sampling_freq,2)/2)*S_freq;(power(sampling_freq,2)/2)*S_freq, sampling_freq*S_freq];

% F=eye(5,5);
% F(4,5)=sampling_time;
Q=zeros(5);
Q(4:5,4:5)=Q1;
 
 H_non_kalman = [H_non_kalman,zeros(6,1)];
 Pa=pinv(H_non_kalman'*H_non_kalman);
 
% Pa=pinv(F'*F);
% P_est= [Pa(1,1), Pa(1,2), Pa(1,3), Pa(1,4),0;
%     Pa(2,1), Pa(2,2), Pa(2,3), Pa(2,4), 0;
%     Pa(3,1), Pa(3,2), Pa(3,3), Pa(3,4), 0;
%     Pa(4,1), Pa(4,2), Pa(4,3), Pa(4,4), 0;
%     0,0,0, 0, 2*Q1(2,2)];

     n=length(x_initial);

R=eye(n)*3*5.22*5.22;

% F=eye(8,8);
% F(4,5)=sampling_time;

F = [1,sampling_time,0,0,0,0,0,0;...
        0,1,0,0,0,0,0,0;...
        0,0,1,sampling_time,0,0,0,0;...
        0,0,0,1,0,0,0,0;...
        0,0,0,0,1,sampling_time,0,0;...
        0,0,0,0,0,1,0,0;...
        0,0,0,0,0,0,1,sampling_time;...
        0,0,0,0,0,0,0,1];
    Tau =   [0.5*sampling_time^2,0,0,0,0;...
        sampling_time,0,0,0,0;...
        0,0.5*sampling_time^2,0,0,0;...
        0,sampling_time,0,0,0;...
        0,0,0.5*sampling_time^2,0,0;...
        0,0,sampling_time,0,0;...
        0,0,0,1,0;...
        0,0,0,0,1];
     Q = Tau*Q*Tau'
% Pa = pinv(h'*h);  % cov of position errors
% P = [Pa(1,1), Pa(1,2), Pa(1,3), Pa(1,4),Pa(1,5),Pa(1,6),Pa(1,7),0;
%     Pa(2,1), Pa(2,2), Pa(2,3), Pa(2,4),Pa(2,5),Pa(2,6),Pa(2,7), 0;
%     Pa(3,1), Pa(3,2), Pa(3,3), Pa(3,4),Pa(3,5),Pa(3,6),Pa(3,7), 0;
%     Pa(4,1), Pa(4,2), Pa(4,3), Pa(4,4),Pa(4,5),Pa(4,6),Pa(4,7), 0;
%     Pa(5,1), Pa(5,2), Pa(5,3), Pa(5,4),Pa(5,5),Pa(5,6),Pa(5,7), 0;
%     Pa(6,1), Pa(6,2), Pa(6,3), Pa(6,4),Pa(6,5),Pa(6,6),Pa(6,7), 0;
%     Pa(7,1), Pa(7,2), Pa(7,3), Pa(7,4),Pa(7,5),Pa(7,6),Pa(7,7), 0;
%     0,0,0, 0, 0, 0,0, 2*Q1(2,2)]

%  P = horzcat(A,zero);
%  rP = [P;0 0 0 0 A(4,4)];
%  Pa=pinv(H_non_kalman'*H_non_kalman);
 P = [Pa(1,1) 0 Pa(1,2) 0 Pa(1,3) 0 Pa(1,4) 0;...
                0 Pa(1,1) 0 Pa(1,2) 0 Pa(1,3) 0      Pa(1,4);...
                Pa(2,1) 0 Pa(2,2) 0 Pa(2,3) 0 Pa(2,4) 0;...
                0 Pa(2,1) 0 Pa(2,2) 0 Pa(2,3) 0      Pa(2,4);...
                Pa(3,1) 0 Pa(3,2) 0 Pa(3,3) 0 Pa(3,4) 0;...
                0 Pa(3,1) 0 Pa(3,2) 0 Pa(3,3) 0      Pa(3,4);...
                Pa(4,1) 0 Pa(4,2) 0 Pa(4,3) 0 Pa(4,4) 0;...
                0      0 0      0 0      0 0 2*Q1(2,2)];

            
for rowIdx=Start_idx:end_idx
     pr_data=allData(rowIdx,pr_indexes);
     Ino_data=allData(rowIdx,Ino_indexes);
     Tropo_data=allData(rowIdx,Tropo_indexes);
     p_x_data= allData(rowIdx,p_x_indexes);
     p_y_data= allData(rowIdx,p_y_indexes);
     p_z_data= allData(rowIdx,p_z_indexes);
     clock_correction_data = allData(rowIdx,clock_correction_indexes);
     
     Final_Pseudorange =  pr_data-Ino_data- Tropo_data+clock_correction_data;
     
     x_initial=p_x_data;
     y_initial=p_y_data;
     z_initial=p_z_data;
%      n= length(x_initial);
     
%    plot(rss_errors);


%    for i=1:n
%    H_matrix=eye(4);
%     x_est=F*x_est;
%     P=F*P;

         
     x_initial_diff=-(p_x_data-x(1));
     y_initial_diff=-(p_y_data-x(2));
     z_initial_diff=-(p_z_data-x(3));
    

    distancegu=sqrt(x_initial_diff.^2+ y_initial_diff.^2+z_initial_diff.^2)+x(4);
    
    del_rho=-distancegu+Final_Pseudorange;

    
    
     inv_distance=1./(distancegu-x(4));
     
     %calculate partial derivtive
     
     ax1_k=x_initial_diff .* inv_distance;
     ay1_k=y_initial_diff .* inv_distance;
     az1_k=z_initial_diff .* inv_distance;
     
%      H(i,:)=[x_initial_diff(i)*inv_distance,y_initial_diff(i)*inv_distance,z_initial_diff(i)*inv_distance,-1];
%     end

%     rss_error = sqrt((x_refrence - x(1))^2+ (y_refrence - x(2))^2 + (z_refrence - x(3))^2);
  
% h1_matrix= [ax1',ay1',az1',ones(6,1)];
%     H=[ax1_k',zeros(6,1),ay1_k',zeros(6,1),az1_k',zeros(6,1),ones(n,1),zeros(6,1)];
H=[ax1_k',ay1_k',az1_k',zeros(6,1),zeros(6,1),zeros(6,1),ones(n,1),zeros(6,1)];

    % H=[H,ones(6,1),zeros(6,1)];
   H_matrix =H;
   P=F*P*F'+Q;
    S=pinv((H_matrix*P*H_matrix')+R);
    K=(P*H_matrix')*S;
%    k_inv=(K*del_rho');
            x=x+(K*del_rho');

    P=(eye(8)-K*H_matrix)*P;

%        
%        Fact=pinv((H_matrix*P*H_matrix')+R);
%        K=(P*H_matrix')*Fact;
      

   

    

 rss_error_ekf = sqrt((x_refrence - x(1))^2+ (y_refrence - x(2))^2 + (z_refrence - x(3))^2);
 
 rss_errors_kalman_lowdynamic(rowIdx-Start_idx+1)=rss_error_ekf;
 
%  disp(rss_errors_kalman);
%  disp(rss_errors);
 
%  rss_errors=rand(1,100);
%  rss_errors_kalman=rand(1,100);
 end
 figure;
 hold on 
 plot(rss_errors,'r-','LineWidth',1.5);
%  hold on ;

 plot(rss_errors_kalman_lowdynamic,'b--','LineWidth',1.5);
 
 xlabel('Index');
 ylabel('Error Value');
 title('Comparison of RSS Errors and Kalman Filter Errors lowdynamic');
 legend('RSS Errors','Kalman Filter Errors');
 grid on;
 
 hold off;
