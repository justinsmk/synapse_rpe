% RPE Model
% t = 0
% v(t) = 0            % Current estimate of the value of a state of an environment
% v(t-1) = 0          % Weighted sum of previous rewards
% alpha = 0.8         % Learning rate constant
% r(t) = 1            % Most recently acquired reward
% RPE = r(t) - v(t-1) %Reward prediction error (RPE)

function RPE()
R_length = 100;
V = zeros (1,R_length);
R = ones(1,R_length);
a = 0.21;

for i=2:R_length
    V(i)= V(i-1)+a*(R(i)-V(i-1));
end
figure(1);hold on;plot(V);
end
    