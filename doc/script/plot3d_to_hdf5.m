clear

%input parameters
gridfile_fortran='sweptwing_coase_973.xyz';
flowfile_fortran='sweptwing_coase_973.fun';
gridfile_hdf5='sweptwing_coase_973_grid.h5';
flowfile_hdf5='sweptwing_coase_973_flow.h5';
curvafile_fortran='curva.dat';
curvafile_ascii = 'curva.plt';
dims = 2;


%read fortran files
fid = fopen(gridfile_fortran, 'rb');
hr = fread(fid, 1, 'int32');
nblk = fread(fid, 1, 'int32');
hr = fread(fid, 1, 'int32');
hr = fread(fid, 1, 'int32');
ishape = fread(fid, dims, 'int32');
hr = fread(fid, 1, 'int32');
ishape= [ishape; dims];
hr = fread(fid, 1, 'int32');
grid = fread(fid, prod(ishape), 'float64');
hr = fread(fid, 1, 'int32');
fclose(fid);
grid = reshape(grid, ishape');
ishape_grid = ishape;

fid = fopen(flowfile_fortran, 'rb');
hr = fread(fid, 1, 'int32');
nblk = fread(fid, 1, 'int32');
hr = fread(fid, 1, 'int32');
hr = fread(fid, 1, 'int32');
ishapes = fread(fid, dims+1, 'int32');
hr = fread(fid, 1, 'int32');
ishapes;
hr = fread(fid, 1, 'int32');
flow = fread(fid, prod(ishapes), 'float64');
hr = fread(fid, 1, 'int32');
fclose(fid);
flow = reshape(flow, ishapes');
ishape_flow = ishapes;

fid = fopen(curvafile_fortran, 'rb');
hr = fread(fid, 1, 'int32');
nblk = fread(fid, 1, 'int32');
hr = fread(fid, 1, 'int32');
hr = fread(fid, 1, 'int32');
ishape = fread(fid, dims-1, 'int32');
hr = fread(fid, 1, 'int32');
hr = fread(fid, 1, 'int32');
ishape = [ishape; 2];
curva = fread(fid, prod(ishape), 'float64');
hr = fread(fid, 1, 'int32');
fclose(fid);
curva = reshape(curva, ishape');
ishape_curva = ishape;

%write hdf5 files
grid_hdf5 = zeros(2, ishape_grid(3), ishape_grid(1), ishape_grid(2));
flow_hdf5 = zeros(2, ishape_flow(3), ishape_flow(1), ishape_flow(2));
for i=1:ishape_grid(3)
    grid_hdf5(1, i, :, :) = grid(:, :, i);
end
for i=1:ishape_flow(3)
    flow_hdf5(1, i, :, :) = flow(:, :, i);
end
if ~exist(gridfile_hdf5,'file')
    h5create(gridfile_hdf5,'/grid', [2, ishape_grid(3), ishape_grid(1), ishape_grid(2)]);
end
h5write(gridfile_hdf5, '/grid', grid_hdf5);
if ~exist(flowfile_hdf5,'file')
    h5create(flowfile_hdf5,'/baseflow', [2, ishape_flow(3), ishape_flow(1), ishape_flow(2)]);
end
h5write(flowfile_hdf5, '/baseflow', flow_hdf5);

%write curva file
fid=fopen(curvafile_ascii,'w');
fprintf(fid,'%.i\n', 1);
if(dims == 2)
    fprintf(fid,'%.i\n', ishape_curva(1));
    for i=1:ishape_curva(1)
        fprintf(fid,'%.7f   %.7f\n', curva(i, 1), curva(i, 2));
    end
elseif(dims ==3)
    fprintf(fid,'%.i %.i\n', ishape_curva(1), ishape_curva(2));
    for k=1:ishape_curva(2)
     for i=1:ishape_curva(1)
        fprintf(fid,'%.7f   %.7f\n', curva(i, k, 1), curva(i, k,  2));
     end
    end
end
