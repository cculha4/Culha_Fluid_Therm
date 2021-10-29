%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Changing Colors for Contour%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mycolormap = colormein(M,clrange1,clrange2,setrangeval,sz)

%M is your plotting matrix
%clrange is the colors you want
if ~exist('sz','var')
    sz = 100;
end
%Do you define where the colors will be set?
if exist('setrangeval','var')

    setrangevalNO = numel(setrangeval)-1; %identifying how many color ranges that need to be identified
    maxM          = max(setrangeval); minM = min(setrangeval); %what are the max and min values of M
%Do you want the code to define the color range
else
    setrangevalNO = 2;
    minM          = min(min(M)); maxM = max(max(M)); medM = mean(mean(M));
    setrangeval   = [ maxM, medM, minM];
end

breaks = zeros(1,setrangevalNO);
mycolormap = zeros(sz,3);
for l = 1:setrangevalNO
    breaks(l)       = [setrangeval(l)-setrangeval(l+1)];
    bks(l)    = floor(abs(breaks(l))/(maxM-minM)*sz);
   if bks(l) == 0
       l = sprintf('the set range is not wide enough between value %d and %d', l, l+1)
   end
end
if sum(bks) ~=sz
    bks(end) = bks(end) + (sz - sum(bks));
end
for l = 1:setrangevalNO
    bksNO = sum(bks(1:l));
    bksR  = bksNO;
    bksL  = bksNO - bks(l) + 1;
    a = linspace(clrange1(l,1),clrange2(l,1),bks(l))';
    b = linspace(clrange1(l,2),clrange2(l,2),bks(l))';
    c = linspace(clrange1(l,3),clrange2(l,3),bks(l))';
    if isnan(bksL) || isnan(bksR)
        mycolormap = clrange1(1,:).*ones(size(mycolormap));
    else
        mycolormap(bksL:bksR,:) = [a b c ]; 
    end
end



