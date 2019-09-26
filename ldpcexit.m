function ldpcexit(H)
   

        dvi = sum(H,1);           % <-- Following 6 lines to compute the average dv and dc as it might be an irregular code
        dci = sum(H,2)';
        [adv,bdv] = hist(dvi,unique(dvi));
        [adc,bdc] = hist(dci,unique(dci));
        fractionEdgesV = (adv.*bdv)/sum(adv.*bdv);
        fractionEdgesC = (adc.*bdc)/sum(adc.*bdc);
        dv = sum((adv./sum(adv)).*bdv);
        dc = sum((adc./sum(adc)).*bdc);
        eL = strcat('Exit-Chart under AWGNCN for LDPC code  with Rate : ',num2str(1-dv/dc));
  
R = 1-(dv/dc);
IAVND = 0:0.001:1;            % Initialize aprior mutual information for variable and check nodes as IAVND and IACND respectively
IACND = 0:0.001:1;
leg={};
col = [1 0 0; 0 0 0; 0 1 1; 1 0 1; 1 1 0; 0.5 0.5 0.5];

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
            strcat('   <---------------------------------- Minimum E_b/N_o for convergence = ',num2str(minEbNodB),' dB'),'FontSize',11,'color',[0 0 0]);
    
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
    leg=[leg,['Converge after ', num2str(iter),' iterations']];



% Adapt figure axis, labels, legends and titles
axis([0 1 0 1]);
xlabel('I_{A,VND},I_{E,CND}');
ylabel('I_{E,VND},I_{A,CND}');
legend(leg,'location','SouthEast');
title(eL);
pbaspect([1 1 1]);
%if cont && r1
ie=0:0.005:1;
ebno = 10^(minEbNodB/10);
index = 1;
for ia = 0:0.005:1
pb(:,index) = 0.5*erfc((sqrt((8*R*ebno)+(JInverseFunction(ie).^2)+(JInverseFunction(ia).^2)))/(2*sqrt(2)));
index = index + 1;
end
hold on;
%clabel(contour(0:0.005:1,0:0.005:1,pb,[0.1 0.08 0.06 0.04 0.03 0.02 0.01 0.005 0.001]),'FontSize',10,'Color','k','Rotation',0);
end
% Create an image file. This png image is displayed on the webpage. Its size
% must be 640x480.



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
c_alpha1 = 2.33727;
a_alpha2 = 0.706692;
b_alpha1 = 0.214217;
b_alpha2 = 0.386013;
c_alpha2 = -1.75017;
JI_temp = zeros(1,length(I));

JI_temp(I >= 0 & I <= I_tp) = (a_alpha1*(I(I >= 0 & I <= I_tp).^2))+(b_alpha1*(I(I >= 0 & I <= I_tp)))+(c_alpha1*(I(I >= 0 & I <= I_tp).^(1/2)));
JI_temp(I > I_tp & I <= 1) = (-a_alpha2*log(-b_alpha2*((I(I > I_tp & I <= 1))-1)))-(c_alpha2*(I(I > I_tp & I <= 1)));

JI = JI_temp;
end