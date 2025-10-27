cell_1= -[14;-2];
cell_2=-[12;8;6;4;2;0];

cell_1= -depth2;
cell_2= -depth1;

oerlap = [];

if size(cell_1,1) > 1 && size(cell_2,1) > 1
    for i1=1:size(cell_1,1)-1
        i2=1;
        a = max(0, - max(cell_1(i1,1), cell_2(i2,1)) + min(cell_1(i1+1,1), cell_2(i2+1,1)));
        
        while a <= 0 && i2 < size(cell_2,1)-1
            i2 = i2+1;
            a = max(0, - max(cell_1(i1,1), cell_2(i2,1)) + min(cell_1(i1+1,1), cell_2(i2+1,1)));
        end
        if a>0
            oerlap = [oerlap;  [i1  i2 a (-cell_1(i1,1)-cell_1(i1+1,1))/2 (-cell_2(i2,1)-cell_2(i2+1,1))/2 min(-cell_1(i1,1), -cell_2(i2,1))-a./2]];
        end
        
        i2_start = i2;
        while a > 0 && i2 < size(cell_2,1)-1
            
            i2 = i2+1;
            a = max(0, - max(cell_1(i1,1), cell_2(i2,1)) + min(cell_1(i1+1,1), cell_2(i2+1,1)));
            if a>0
                oerlap = [oerlap;  [i1  i2 a (-cell_1(i1,1)-cell_1(i1+1,1))/2 (-cell_2(i2,1)-cell_2(i2+1,1))/2 min(-cell_1(i1,1), -cell_2(i2,1))-a./2]];
            end
        end
        i2 = i2_start;
    end
end