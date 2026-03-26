function [A,B,C,D,K] = plant2ABCDK(plant)
    [A,BK,C,DI] = ssdata(plant);
    ny = size(C,1);
    B = BK(:,1:end-ny);
    K = BK(:,end-ny+1:end);
    nu = size(B,2);
    D = DI(:,1:nu);
end