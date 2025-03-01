clear

cd ..
cd ..
cd Datasets\MANC\Synapses

files = dir("*.txt");

filenames = [];
for i=1:size(files,1)
    disp(i)
    filename = files(i).name;
    filename = filename(1:end-4);
    filenames = [filenames;string(filename)];
end

Weighted_A = zeros(length(filenames),length(filenames));

for i=1:length(filenames)
    disp(i)
    subtable = readtable(strcat(filenames(i),".txt"));
    partners = int64(double(subtable{:,3}));
    [Lia,Lib] = ismember(string(partners),filenames);
    for j=1:length(Lib)
        if Lib(j)~=0
            Weighted_A(i,Lib(j)) = Weighted_A(i,Lib(j)) + 1;
            %Weighted_A(Lib(j),i) = Weighted_A(Lib(j),i) + 1;
        end
    end
end

Weighted_A2 = sparse(Weighted_A);

[row,col,v] = find(Weighted_A);

v = v(randperm(numel(v)));

Shuffled_A = sparse(row, col, v, size(Weighted_A, 1), size(Weighted_A, 2));

Unweighted_A = sparse(row, col, ones(length(v),1), size(Weighted_A, 1), size(Weighted_A, 2));

Degree_List = horzcat(ones(size(Weighted_A,1),1),full(sum(Unweighted_A,2)),full(sum(Weighted_A,2)));

cd ..
cd ..
cd ..
cd ProcessedData/MANC

save("MANC_Degree_Data.mat",'Degree_List');