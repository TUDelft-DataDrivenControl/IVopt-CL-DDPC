function yr = make_reference(Nbar_ref,ny)
% This function replicates a reference trajectory from
%
% [1] A. Chiuso, M. Fabris, V. Breschi, and S. Formentin, “Harnessing
%     uncertainty for a separation principle in direct data-driven
%     predictive control,” Automatica, vol. 173, p. 112070, Mar. 2025,
%     doi: 10.1016/j.automatica.2024.112070.
%
% using code that acompanies [1] that can be found at
%     https://github.com/marcofabris92/a-separation-principle-in-d3pc

arguments
   Nbar_ref (1,1) double
   ny (1,1) double
end

%% create original reference
yr = zeros(1,500);
for t = 1:10
    yr(t) = -1;
end
for t = 11:20
    yr(t) = 10;
end
for t = 21:30
    yr(t) = -1;
end
for t = 31:40
    yr(t) = 1;
end
for t = 41:180
    yr(t) = -10;
end
for t = 181:320
    yr(t) = 10;
end
for t = 321:460
    yr(t) = -1;
end
for t = 461:470
    yr(t) = 1;
end
for t = 471:480
    yr(t) = -10;
end
for t = 481:490
    yr(t) = 1;
end
for t = 491:500
    yr(t) = -1;
end

%% replication of the original reference
reps = ceil(Nbar_ref/500);
yr = repmat(yr,ny,reps);
yr = yr(:,1:Nbar_ref);

end