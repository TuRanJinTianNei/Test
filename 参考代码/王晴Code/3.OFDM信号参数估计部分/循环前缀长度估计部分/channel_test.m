function my_channel = channel_test(N_CP)
    delay = round([0 0.04 0.09 0.25 0.36 0.81]*N_CP);
    gains_dB = [0 -2 -3 -6 -8 -10];
    gains = 10.^(gains_dB/10);
    my_channel = zeros(1, delay(end) + 1);
    my_channel(delay + 1) = sqrt(gains/2) .* (randn(1,length(delay)) + 1j * randn(1,length(delay)));
end