%% load 
% load('..\results\APAS_result.mat')

%% Extract pressure magnitude and phase (quick access of single simulation)
% a_index = 1;
% r_index = 1;
% z_index = 1;
% f_index = 1;

% P_level = 20*log10(abs(freq_spectrum{a_index,r_index,z_index,f_index}));
% P_phase = angle(freq_spectrum{a_index,r_index,z_index,f_index});

%% Extract pressure magnitude and phase of z and f simulations
% pressure = squeeze(cell2mat(freq_spectrum(1,1,:,:)));
% % plot frequency spectrum for first z-distance
% % plot(f_vec,20*log10(abs(pressure(:,1))))
% figure;plot(f_vec,abs(pressure(:,1)))

%% Extract pressure magnitude and phase of a and r simulations
% pressure = squeeze(cell2mat(pressure_freq_spectrum(:,:,1,1)));

%% Investigate transmission coeffiscient
% eta = 0:1e-2:2500;
% T_coeff = transducerFunctions.TransmissionCoeff(eta);
% figure
% plot(etaPlot,real(T_coeff(etaPlot)));
%% Investigate integrand
eta1 = 0:1e-2:2500;
% eta2 = 1936.716:1e-6:1936.722;
% eta2 = 0:1e-2:2500;

y_vec1 = integrand_org(eta1);
% y_vec2 = integrand_parts(eta2);
% y_vec2 = integrand_org(eta2);

% figure
% plot(eta1,imag(integrand_org(eta1)));

% save test.mat eta1 eta2 y_vec1 y_vec2 C_INT
% save test.mat eta1 y_vec1 y_vec2 C_INT



