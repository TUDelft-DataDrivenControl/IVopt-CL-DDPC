function yr = make_reference(Nbar_ref,ny)
arguments
   Nbar_ref (1,1) double
   ny (1,1) double
end
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
reps = ceil(Nbar_ref/500);
yr = repmat(yr,ny,reps);
yr = yr(:,1:Nbar_ref);

end