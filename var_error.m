function [ error_top, error_bottom ] = var_error(variance, Ncycles, cycle_length, conf_interval)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Generating the statisitcal error
%conf_interval = 95; % Confidence interval in percentage
alpha = 1 - (conf_interval/100);
DF = Ncycles-1; % degrees of freedom of the chi-squared distribution

% calculate the chi-squared values from table:
chi_left = chi2inv(1-(alpha/2),DF);
chi_right = chi2inv((alpha/2),DF);

%sumx=zeros(1,cycle_length-1);
%sumx2=zeros(1,cycle_length-1);
a=zeros(cycle_length-1,1);
ci_bottom=zeros(cycle_length-1,1);
ci_top=zeros(cycle_length-1,1);
error_bottom=zeros(cycle_length-1,1);
error_top=zeros(cycle_length-1,1);
for j=1:cycle_length-1
    %sumx(j) = sum(newb(j,:));
    %sumx2(j) = sum(newb(j,:).^2);
    a(j) = DF*(variance(j));
    ci_bottom(j) = a(j)/chi_left;
    ci_top(j) = a(j)/chi_right;
    error_top(j) = ci_top(j)-variance(j);
    error_bottom(j) = variance(j) - ci_bottom(j);
end
end

