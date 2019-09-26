%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% .m script of web demo "LDPC EXIT Chart"
% Written for GNU Octave
%
% INSTITUTE OF TELECOMMUNICATIONS
% University of Stuttgart
% www.inue-uni-stuttgart.de
% author: Sherif Abdulatif
% date: 01.09.2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ldpcexit_flex(ch,Gdv,Gdc,r1,hr,cont,uploadfile,filename)
% Function ldpcExit used to plot EXIT Chart of an LDPC Code for different
% degree profiles (Regular and Irregular)
% *ch is a controller of the main function as follows:
%        1- A given degree profiles (Gdv,Gdc) by the user
%        2- Wifi standard code
%        3- DVB-S2 standard code
%        4- Wimax standard code
% *uploadfile is activated only for hr=1 and it's the path of the .txt H-matrix
% *Gdv and Gdc are activated only for ch=1 and they are the degree profiles
% given by the user
% *r1 is a controller to choose between different channels as follows:
%        0- Binary Errasure Channel (BEC) Model with errasure probability
%           given as variable *bec
%        1- Additive White Gaussian Noise (AWGN) Model with Eb/No given as
%           a variable *EbNo
% *hr is a controller to control whether a matrix is uploaded or not
% *cont is a controller to control whether BER contour plot is displayed or not
% *uploadfile link to upladed file
% *filename file to save in
% Additional parameters that can be found in the code:
% *R is the rate of the code given as (1-dv/dc) the higher the rate, the
% less the redundant bits
% *IAVND is the aprior mutual information of variable nodes input part
% (Plotted on x-axis)
% *IEVND is the extrinsic mutual information is the extrinsic or output
% mutual information of variable nodes (Plotted on y-axis)
% *IACND is the aprior mutual information of check nodes input part
% (Plotted on y-axis)
% *IECND is the extrinsic mutual information is the extrinsic or output
% mutual information of check nodes (Plotted on x-axis)
% Additional information: the only constraint in the code is that dv must
% never be greater than dc in that case the rate will be negative and code
% will not accept it
% Any other information is commented in the code, in addition to extra 4
% functions at the end
fh = figure(1);
if hr == 1
    try
        H = load(uploadfile);       % Load the parity check matrix (H) from a given path
    catch
        ldpcExit(ch,Gdv,Gdc,r1,0,uploadfile,filename);
    end
    if all(all((H==1|H==0)))
        H = load(uploadfile);           % Load the parity check matrix (H) from a given path
        dvi = sum(H,1);           % <-- Following 6 lines to compute the average dv and dc as it might be an irregular code
        dci = sum(H,2)';
        [adv,bdv] = hist(dvi,unique(dvi));
        [adc,bdc] = hist(dci,unique(dci));
        fractionEdgesV = (adv.*bdv)/sum(adv.*bdv);
        fractionEdgesC = (adc.*bdc)/sum(adc.*bdc);
        dv = sum((adv./sum(adv)).*bdv);
        dc = sum((adc./sum(adc)).*bdc);
        eL = strcat('Exit Curve for LDPC of specified H-Matrix with Rate: ',num2str(1-dv/dc));
    else
        ldpcExit(ch,Gdv,Gdc,r1,0,uploadfile,filename);
    end
else
    if ch == 1                     % Certain degree profile
        bdv = Gdv;                 % Directly place given degree profiles to code
        bdc = Gdc;
        dv = bdv;
        dc = bdc;
        fractionEdgesV = 1;
        fractionEdgesC = 1;
        eL = strcat('Exit Curve for LDPC of certain degree profile with Rate : ',num2str(1-dv/dc));
    elseif ch == 2               % Wifi Standard
        H = load('HWifi.mat');           % Load the parity check matrix (H) from a given path
        H = H.H;
        dvi = sum(H,1);           % <-- Following 6 lines to compute the average dv and dc as it might be an irregular code
        dci = sum(H,2)';
        [adv,bdv] = hist(dvi,unique(dvi));
        [adc,bdc] = hist(dci,unique(dci));
        fractionEdgesV = (adv.*bdv)/sum(adv.*bdv);
        fractionEdgesC = (adc.*bdc)/sum(adc.*bdc);
        dv = sum((adv./sum(adv)).*bdv);
        dc = sum((adc./sum(adc)).*bdc);
        eL = strcat('Exit Curve for LDPC of Wifi IEEE 802.11n with Rate : ',num2str(1-dv/dc));
    elseif ch == 3             % DVB-S2 Standard
        H = load('HDVB_910.mat');           % Load the parity check matrix (H) from a given path
        H = H.H;
        H = full(H);
        dvi = sum(H,1);           % <-- Following 6 lines to compute the average dv and dc as it might be an irregular code
        dci = sum(H,2)';
        [adv,bdv] = hist(dvi,unique(dvi));
        [adc,bdc] = hist(dci,unique(dci));
        fractionEdgesV = (adv.*bdv)/sum(adv.*bdv);
        fractionEdgesC = (adc.*bdc)/sum(adc.*bdc);
        dv = sum((adv./sum(adv)).*bdv);
        dc = sum((adc./sum(adc)).*bdc);
        eL = strcat('Exit Curve for LDPC of DVB-S2 with Rate : ',num2str(1-dv/dc));
    elseif ch == 4            % Wimax Standard
        H = load('HWimax.mat');           % Load the parity check matrix (H) from a given path
        H = H.H;
        dvi = sum(H,1);           % <-- Following 6 lines to compute the average dv and dc as it might be an irregular code
        dci = sum(H,2)';
        [adv,bdv] = hist(dvi,unique(dvi));
        [adc,bdc] = hist(dci,unique(dci));
        fractionEdgesV = (adv.*bdv)/sum(adv.*bdv);
        fractionEdgesC = (adc.*bdc)/sum(adc.*bdc);
        dv = sum((adv./sum(adv)).*bdv);
        dc = sum((adc./sum(adc)).*bdc);
        eL = strcat('Exit Curve for LDPC of WiMax IEEE 802.16e with Rate : ',num2str(1-dv/dc));
    end
end
R = 1-(dv/dc);
IAVND = 0:0.001:1;            % Initialize aprior mutual information for variable and check nodes as IAVND and IACND respectively
IACND = 0:0.001:1;
leg={};
col = [1 0 0; 0 0 0; 0 1 1; 1 0 1; 1 1 0; 0.5 0.5 0.5];
if r1
    IECNDi=[];                                                         % Calculate extrinsic check node mutual information according to AWGN model
    for i = 1:length(bdc)
        IECNDi = [IECNDi ; fractionEdgesC(i)*(1-JFunction(sqrt(bdc(i)-1)*JInverseFunction(1-IACND)))];  %AWGN the .^2 ????
    end
    IECND = sum(IECNDi,1);       %AWGN
    plot(IECND,IACND,'Color',[0,0.8,0],'linewidth',2);                % Plot check nodes (Swapped axis)
    leg = [leg, ['dc_{average} = ', num2str(dc)]];
    hold on;
    minEbNodB = 0;                                                     % Initial minimum EbNodB
    for EbNodB = 0:0.1:5                                                % Loop over different EbNodB values
        EbNo = 10^(EbNodB/10);                                          % Convert to linear scale
        [x1,y1,iter] = convergence(bdv,bdc,fractionEdgesV,fractionEdgesC,R,EbNo,1);   % Function convergenceAWGN return stair case converging plot for AWGN model
        if(~isempty(x1) && ~isempty(y1) && x1(end) > 0.95 && y1(end) > 0.95)    % Check for convergence mutual information reached accepted value
            minEbNodB = EbNodB;                                           % Set the minimum
            break;                                                      % Break the loop when minimum is found
        end
    end
    p = [];                                                    % Markers for different EbNodB values
    for i = 1:length(bdv)
        p = [p ; fractionEdgesV(i)*(JFunction(sqrt(((bdv(i)-1)*(JInverseFunction(0.04).^2))+(8*R*10^(minEbNodB/10)))))]; %AWGN
    end
    p = sum(p,1);  %AWGN
    text(0.04,p,...
            strcat('<--------------------------- Minimum E_b/N_o for convergence = ',num2str(minEbNodB),' dB'),'FontSize',8,'color',col(2,:));
    index = 1;
    for EbNodB = minEbNodB-1:minEbNodB+4                                                 % Loop over different EbNodB values
        EbNo = 10^(EbNodB/10);                                           % Convert to linear scale
        hold on;
        IEVNDi = [];                                                    % Compute extrinsic mutula information for AWGN formula
        for i = 1:length(bdv)
            IEVNDi = [IEVNDi ; fractionEdgesV(i)*(JFunction(sqrt(((bdv(i)-1)*(JInverseFunction(IAVND).^2))+(8*R*EbNo))))]; %AWGN
        end
        IEVND = sum(IEVNDi,1);  %AWGN
        plot(IAVND,IEVND,'color',col(index,:),'linewidth',2);
        if index == 1
        leg = [leg, ['dv_{average} = ', num2str(dv), ' for E_b/N_o = ' , num2str(EbNodB) , ' dB']];
        else 
        leg = [leg, ['E_b/N_o = ' , num2str(EbNodB) , ' dB']];        
        end
        index = index +1;
    end
    plot(x1,y1);                                             % Plot stair case converging plot
    leg=[leg,['Staircase with ', num2str(iter),' iterations']];
else
	 if ch == 3
	 ldpcExit(3,0,0,1,0,0,'',filename);
	 else
    IECNDi=[];                                                           % Calculate extrinsic check node mutual information according to BEC model
    for i = 1:length(bdc)
        IECNDi = [IECNDi ; fractionEdgesC(i)*((IACND).^(bdc(i)-1))];   %BEC
    end
    IECND = sum(IECNDi,1);       %BEC
    plot(IECND,IACND,'Color',[0,0.8,0],'linewidth',2);                % Plot check nodes (Swapped axis)
    leg = [leg, ['dc_{average} = ', num2str(dc)]];
    hold on;
    maxErrasure = 0;                                                     % Initial minimum EbNodB
    for bec = 0.5:-0.05:0.1                                               % Loop over different EbNodB values
        [x1,y1,iter] = convergence(bdv,bdc,fractionEdgesV,fractionEdgesC,R,bec,0);   % Function convergenceAWGN return stair case converging plot for AWGN model
        if(~isempty(x1) && ~isempty(y1) && x1(end) > 0.95 && y1(end) > 0.95)    % Check for convergence mutual information reached accepted value
            maxErrasure = bec;                                           % Set the minimum
            break;                                                      % Break the loop when minimum is found
        end
    end
    p = [];                                              % Markers for different errasure values
    for i = 1:length(bdv)
         p = [p ; fractionEdgesV(i)*(1-(maxErrasure*((1-0.03).^(bdv(i)-1))))];  %BEC
    end
    p = sum(p,1);  %BEC
    text(0.032,p,...
                strcat('<--------------------------- Maximum \epsilon for convergence = ',num2str(maxErrasure)),'FontSize',8,'color',col(2,:));
    loopTerm = maxErrasure-0.4;
    if loopTerm < 0 
    loopTerm = 0;
    end
    index = 1;
    for bec = maxErrasure+0.1:-0.1:loopTerm                                               % Loop over different errasure values
        hold on;
        IEVNDi = [];                                                        % Compute extrinsic mutula information for BEC formula
        for i = 1:length(bdv)
            IEVNDi = [IEVNDi ; fractionEdgesV(i)*(1-(bec*((1-IAVND).^(bdv(i)-1))))];  %BEC
        end
        IEVND = sum(IEVNDi,1);  %BEC
        plot(IAVND,IEVND,'color',col(index,:),'linewidth',2);
        if index == 1
        if bec == 0
        leg = [leg, ['dv_{average} = ', num2str(dv), ' for \epsilon = 0.0']];
        else
        leg = [leg, ['dv_{average} = ', num2str(dv), ' for \epsilon = ' , num2str(bec)]];
        end
        else 
        if bec == 0 || bec < 1e-3
        leg = [leg, ['\epsilon = 0.0']];
        else
        leg = [leg, ['\epsilon = ' , num2str(bec)]];
        end
        end
        index = index +1;
    end
    plot(x1,y1);                                             % Plot stair case converging plot
    leg=[leg,['Staircase with ', num2str(iter),' iterations']];
    end
end
% Adapt figure axis, labels, legends and titles
axis([0 1 0 1]);
xlabel('I_{A,VND},I_{E,CND}');
ylabel('I_{E,VND},I_{A,CND}');
legend(leg,'location','SouthEast');
title(eL);
pbaspect([1 1 1]);
if cont && r1
ie=0:0.005:1;
ebno = 10^(minEbNodB/10);
index = 1;
for ia = 0:0.005:1
pb(:,index) = 0.5*erfc((sqrt((8*R*ebno)+(JInverseFunction(ie).^2)+(JInverseFunction(ia).^2)))/(2*sqrt(2)));
index = index + 1;
end
hold on;
clabel(contour(0:0.005:1,0:0.005:1,pb,[0.1 0.08 0.06 0.04 0.03 0.02 0.01 0.005 0.001]),'FontSize',10,'Color','k','Rotation',0,'labelspacing',700);
end
% Create an image file. This png image is displayed on the webpage. Its size
% must be 640x480.
print(filename, "-dgif","-S640,480", "-r0");
end

function [x,y,iteration] = convergence(bdv,bdc,fractionEdgesV,fractionEdgesC,R,conv,tch)
% Function convergence is used to plot the stair case converging,
% (x1,y1) (x2,y2) (x3,y3) are x-y iterative values in the tunnel
% between the check node and variable node curves and return the
% stair case plot as x-y
iteration = 0;
x1 = 0;
y1 = 0;
x_temp=[];
y_temp=[];
while true
    if tch                                                       % Used in AWGN case (Setting y2 and x3 to AWGN model values)
        y2 = [];
        for i = 1:length(bdv)
            y2 = [y2 ; fractionEdgesV(i)*(JFunction(sqrt(((bdv(i)-1)*(JInverseFunction(x1).^2))+(8*R*conv))))];
        end
        y2 = sum(y2,1);
        x3=[];
        for i = 1:length(bdc)
            x3 = [x3 ; fractionEdgesC(i)*(1-JFunction(sqrt(bdc(i)-1)*JInverseFunction(1-y2)))];  %AWGN
        end
        x3 = sum(x3,1);       %AWGN
    else                                                         % Used in BEC case (Setting y2 and x3 to BEC model values)
        y2 = [];
        for i = 1:length(bdv)
            y2 = [y2 ; fractionEdgesV(i)*(1-(conv*((1-x1).^(bdv(i)-1))))];
        end
        y2 = sum(y2,1);
        x3=[];
        for i = 1:length(bdc)
            x3 = [x3 ; fractionEdgesC(i)*((y2).^(bdc(i)-1))];
        end
        x3 = sum(x3,1);
    end
    x2 = x1;                                                     % Update x2 with x1
    if y2 - y1 < 10e-5                                           % Stop when the difference between y2 and y1 is too small
        break;
    end
    y3 = y2;                                                     % Update y3 with y2
    if x3 -x1 < 10e-5                                            % Stop when the difference between x3 and x1 is too small
        break;
    end
    x_temp = [x_temp x1 x2 x3];                                  % Fill returned arrays
    y_temp = [y_temp y1 y2 y3];
    x1 = x3;                                                     % Update x1 and y1 for next iteration
    y1 = y3;
    iteration = iteration + 1;
end
x = x_temp;                                                      % Return the converging curve
y = y_temp;
end

function [J] = JFunction(alpha)
% JFunction is an approximated function for EXIT curves of AWGN
% channel model, it takes alpha vector and returns a corresponding
% JFunction vector.
% Values are arbitary that fullfill required function
alpha_tp = 1.6363;
a_J1 = -0.0421061;
b_J1 = 0.209252;
c_J1 = -0.00640081;
a_J2 = 0.00181491;
b_J2 = -0.142675;
c_J2 = -0.0822054;
d_J2 = 0.0549608;
J_temp = zeros(1,length(alpha));

J_temp(alpha >= 0 & alpha <= alpha_tp) = (a_J1*(alpha(alpha >= 0 & alpha <= alpha_tp).^3))+(b_J1*(alpha(alpha >= 0 & alpha <= alpha_tp).^2))+(c_J1*(alpha(alpha >= 0 & alpha <= alpha_tp)));
J_temp(alpha > alpha_tp & alpha < 10) = 1-exp((a_J2*(alpha(alpha > alpha_tp & alpha < 10).^3))+(b_J2*(alpha(alpha > alpha_tp & alpha < 10).^2))+(c_J2*(alpha(alpha > alpha_tp & alpha < 10)))+(d_J2));
J_temp(alpha >= 10) = 1;

J = J_temp;
end

function [JI] = JInverseFunction(I)
% JInverseFunction is an approximated function for EXIT curves of
% BEC channel model, it takes mutual information vector and returns
% a corresponding JInverseFunction vector.
% Values are arbitary that fullfill inverse of Jfunction above
I_tp = 0.3646;
a_alpha1 = 1.09542;
b_alpha1 = 0.214217;
c_alpha1 = 2.33727;
a_alpha2 = 0.706692;
b_alpha2 = 0.386013;
c_alpha2 = -1.75017;
JI_temp = zeros(1,length(I));

JI_temp(I >= 0 & I <= I_tp) = (a_alpha1*(I(I >= 0 & I <= I_tp).^2))+(b_alpha1*(I(I >= 0 & I <= I_tp)))+(c_alpha1*(I(I >= 0 & I <= I_tp).^(1/2)));
JI_temp(I > I_tp & I <= 1) = (-a_alpha2*log(-b_alpha2*((I(I > I_tp & I <= 1))-1)))-(c_alpha2*(I(I > I_tp & I <= 1)));

JI = JI_temp;
end
