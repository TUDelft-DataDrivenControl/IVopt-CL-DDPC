function Re = snr2Re(plant,SNR_goal,Ru,Re,Nsim)
SNR_0 = calc_snr(plant,Re=Re,Ru=Ru,Nsim=Nsim);
SNR_ratio = (SNR_0-1)/(SNR_goal-1);
Re = Re*SNR_ratio;
end