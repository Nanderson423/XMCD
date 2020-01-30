function [y_background,varargout] = background(Energy,XAS,varargin)
%[Y_BACKGROUND, step1, step2] = background(ENERGY,XAS,[step1 Erange],[step2
%Erange]) OR [Y_BACKGROUND, m, b] = background(ENERGY,XAS,[Linear Range])
%Background - Creates the background for the input data
%Has 2 required inputs and 2 optional inputs
%Only 2 inputs gives a linear background. 4 inputs will give the step
%background.
%Required inputs:
%   Energy - The energy vector corresponding to the XAS
%   XAS - The intensity of the XAS in a vector
%Optional inputs:
%   The extra two inputs are required for the step background.
%   First varargin: The range of the first edge as a vector (i.e. [700 715])
%   Second varagin: The range of the second edge as a vector (i.e. [715
%   730])
%Outputs:
%   First output is the intensity of the background (y_background)
%   Optional outputs are    1) m of linear or first step location
%                           2) b of linear of second step location

switch nargin
    
    case 3 %For Linear Background
    [~,i1] = min(abs(Energy-varargin{1}(1)));
    [~,i2] = min(abs(Energy-varargin{1}(2)));
    m = (XAS(i2)-XAS(i1))/(Energy(i2)-Energy(i1)); %Slope of the backgorund
    b = XAS(1) -  m*Energy(1); %Y-intercept of the background
    y_background =  m.*Energy + b; %Gets y values for each energy of the background
    varargout{1} = m;
    varargout{2} = b;
    
    
    case 5 %For Steps Background
    %STEPS
    Step1 = varargin{1};
    Step2 = varargin{2};
    Scale = varargin{3};
    
    [~,bLow] = min(abs(Energy-(Step2(2)-.5)));
    bHigh = length(Energy);
    [bMin,bMinI] = min(XAS(bLow:bHigh)); % For adjusting end backround height (min/max/average)
    Step2(2) = Energy(bMinI+bLow-1);
    
    
    %First Step
    [~,bLow] = min(abs(Energy-Step1(1)));
    [~,bHigh] = min(abs(Energy-Step1(2)));
    [L3peak,L3peakI] = max(XAS(bLow:bHigh));
    [~,L3peakL] = min(abs(XAS(bLow:L3peakI+bLow)-L3peak/2-bMin/6));
    L3peakL = L3peakL + bLow-1;
    [~,L3peakH] = min(abs(XAS(L3peakI+bLow:bHigh)-L3peak/2-bMin/6));
    L3peakH = L3peakH + bLow + L3peakI - 1;
    step1 = round((L3peakH+L3peakL)/2);
    
    
    %Second Step
    [~,bLow] = min(abs(Energy-Step2(1)));
    [~,bHigh] = min(abs(Energy-Step2(2)));
    [L2peak,L2peakI] = max(XAS(bLow:bHigh));  
    [~,L2peakL] = min(abs(XAS(bLow:L2peakI+bLow)-L2peak/2-(5/12)*bMin));
    L2peakL = L2peakL + bLow - 1;
    [~,L2peakH] = min(abs(XAS(L2peakI+bLow:bHigh)-L2peak/2-(5/12)*bMin));
    L2peakH = L2peakH + bLow + L2peakI - 1;
    step2 = round((L2peakH+L2peakL)/2);
    
    varargout{1} = step1;
    varargout{2} = step2;
    
    y_background = zeros(size(XAS));
    y_background(step1:step2-1) = bMin*Scale;
    y_background(step2:end) = bMin;
    
    %disp('end')
end