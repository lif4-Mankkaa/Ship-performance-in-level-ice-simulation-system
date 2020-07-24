function [y1] = NN1(x1)
% Neural network simulation function for stress calculation
%
% ALL PARAMETERS ARE NON-DIMENSIONAL
%
% [y1] = NN1(x1) takes these arguments:
%   x = Qx5 matrix, input #1, sequence of inputs: [theta (deg), b/lc, c/lc, mu1,mu2]
% and returns:
%   y = Qx1 matrix, output #1, y = maximum princ stress / pressure on contact area
% where Q is the number of samples.
% 

%#ok<*RPMT0>
lc = (50e6*0.03^3/(12*0.91*1025*9.81))^0.25*1000;
x1(2:3) = x1(2:3)*lc;
x1(4:5) = x1(4:5)*1000;

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [10;6;15;50;50];
x1_step1.gain = [0.0133333333333333;0.00336700336700337;0.00493827160493827;0.0133333333333333;0.0133333333333333];
x1_step1.ymin = -1;

% Layer 1
b1 = [2.49213365350507;-1.78852461376759;-0.237760868117081;0.482093835055585;1.23929446572885;1.39181876428233;-1.02584511213037;-2.15291006542418;1.72678245881189;-4.03319426672703];
IW1_1 = [0.202459488611851,-0.219729798408599,1.38245467005947,-0.385876119203726,-0.0285611714723076;0.659207303546072,-1.72103419052228,0.263509931068996,-0.847299221509922,0.556850059309043;-0.197176152684097,0.284267977873888,0.706775315864802,0.391391795441077,0.0146023980749638;0.346810653635543,1.00247123837193,0.483923179400450,-0.210806615315506,-0.0155882152905635;0.317384516119153,-0.140057262528129,-0.597533564898498,0.583379841265598,0.0315133806032032;0.322317555073020,0.0282908619626646,-0.362946121577700,-0.532398687343130,-0.0359380765077768;-0.258290776625427,-1.13691296343142,-0.940593975064239,0.0512330224085031,0.0130774211372146;-0.777087262224913,0.443731451773981,-0.513725586580340,0.303635552845485,0.0766238952044174;0.604520341878937,1.00126732949522,-0.872621942401602,0.981729168601916,1.33764195886737;-0.506769873769450,-2.79972597191456,0.373641358335232,0.429147313103166,0.0253847828626398];
% Layer 2
b2 = -1.38457456564762;
LW2_1 = [1.48325129247911,0.00837637882607949,0.486792458057967,-0.421602160189142,0.514129712311897,-1.02849604444290,-0.393474290623295,0.912243272727443,0.00400052095972890,-0.933764986489913];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 388.207725917801;
y1_step1.xoffset = 3.35456970788073e-05;

% ===== SIMULATION ========

% Dimensions
Q = size(x1,1); % samples

% Input 1
x1 = x1';
xp1 = mapminmax_apply(x1,x1_step1);

% Layer 1
a1 = tansig_apply(repmat(b1,1,Q) + IW1_1*xp1);

% Layer 2
a2 = repmat(b2,1,Q) + LW2_1*a1;

% Output 1
y1 = mapminmax_reverse(a2,y1_step1);
y1 = y1'/0.001;
end

% ===== MODULE FUNCTIONS ========

% Map Minimum and Maximum Input Processing Function
function y = mapminmax_apply(x,settings)
y = bsxfun(@minus,x,settings.xoffset);
y = bsxfun(@times,y,settings.gain);
y = bsxfun(@plus,y,settings.ymin);
end

% Sigmoid Symmetric Transfer Function
function a = tansig_apply(n,~)
a = 2 ./ (1 + exp(-2*n)) - 1;
end

% Map Minimum and Maximum Output Reverse-Processing Function
function x = mapminmax_reverse(y,settings)
x = bsxfun(@minus,y,settings.ymin);
x = bsxfun(@rdivide,x,settings.gain);
x = bsxfun(@plus,x,settings.xoffset);
end