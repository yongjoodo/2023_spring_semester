function [x,resid]=qrsolve(A,b,mtd)

    if(mtd=="clqrgrsch")
        n=length(A);
        [Q,R]=clqrgrsch(A);
        tempb=Q'*b;
        x=zeros(n,1);
        x(n,1)=tempb(n,1)/R(n,n);
        for i=1:n-1
            x(n-i,1)=(tempb(n-i,1)-R(n-i,n-i+1:n)*tempb(n-i+1:n,1))/R(n-i,n-i);
        end
        resid=norm(A*x-b,2);
        fprintf("Selected method: clqrgrsch");
    end

    if(mtd=="modqrgrsch")
        n=length(A);
        [Q,R]=modqrgrsch(A);
        tempb=Q'*b;
        x=zeros(n,1);
        x(n,1)=tempb(n,1)/R(n,n);
        for i=1:n-1
            x(n-i,1)=(tempb(n-i,1)-R(n-i,n-i+1:n)*tempb(n-i+1:n,1))/R(n-i,n-i);
        end
        resid=norm(A*x-b,2);
        fprintf("Selected method: modqrgrsch");
    end

    if(mtd=="qr")
        n=length(A);
        [Q,R]=qr(A);
        tempb=Q'*b;
        x=zeros(n,1);
        x(n,1)=tempb(n,1)/R(n,n);
        for i=1:n-1
            x(n-i,1)=(tempb(n-i,1)-R(n-i,n-i+1:n)*tempb(n-i+1:n,1))/R(n-i,n-i);
        end
        resid=norm(A*x-b,2);
        fprintf("Selected method: qr");
    end

end


