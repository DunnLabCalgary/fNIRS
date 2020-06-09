function [d] = inf_resolution_d (d)

%%%%%%%%%%%remove infs%%%%%%%%%%%
[a,b]=find(isinf(d));
if isempty(a)==0
    for i=1:length(a)
        d(a(i),b(i))=nan;
        %d(:,b(i))=inpaint_nans(d(:,b(i)));
        if a(i)~=1
            d(a(i),b(i))=d(a(i)-1,b(i));
        else
            d(a(i),b(i))=d(a(i)+1,b(i));
        end
    end
end
d =inpaint_nans(d);
%ending inf resolution

end