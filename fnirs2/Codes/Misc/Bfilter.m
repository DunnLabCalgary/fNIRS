function filt=Bfilter(Xres)
for j=1:Xres
    rx=(j-Xres/2);
    f(j)=1/(1+exp((rx-0.25*Xres)*(24/Xres)));
end
for j=1:Xres
    filt(j)=f(Xres-j+1);
    filt(j+Xres)=f(j);
end

% plot(filt)
end