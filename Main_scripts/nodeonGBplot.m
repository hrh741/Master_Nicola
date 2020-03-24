
n0 = []; n1=[]; np=[]; np1=[]; np2=[];
figure
hold on;

 for ng = 1:ngr
    n0 = nodeonGB{ng,1};
    n1 = n0(:,1);
    np = p(n1,:);
    for i = 1:size(n1,1)
        np1 = np(i,:);
        plot(np1(1), np1(2),'ro');
        axis([-1 5 -1 3]);
        title('nodeonGB')
        pause
    end
    
%     plot(np(:,1), np(:,2),'ro');
end

% nGB=[]; nGB1=[]; nGBp=[];
% for ng = 1:ngr
%     nGB = nodeonGB{ng,1};
%     nGB1 = nGB(:,1);
%     for i = 1:size(nGB,1)
%         nGBp = p(nGB1,i);
%         plot(nGBp(:,1), nGBp(:,2))
%         pause
%     end
% end
