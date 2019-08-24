function plot_local_XYZ_bloch_spheres(local_XYZ_series,axis_labels,f)
    if nargin==2
        f = 0.5; % scale factor for shadow sphere
    elseif nargin==1
        f = 0.5;
        axis_labels = {'Sx','Sy','Sz'};
    end

    figure('Position',[500 400 1000 400])    
    N = size(local_XYZ_series,1) / 3;
    nRows = floor(sqrt(N));
    nCols = ceil(N/nRows);
    colormap gray
    
    if max(abs(imag(local_XYZ_series(:)))) > 1000*eps
        warning('local_XYZ_series has imaginary parts > 1000*eps')
    end
    local_XYZ_series = real(local_XYZ_series);
    
    for k=1:N
        subplot(nRows,nCols,k)
        hold on
        
        xidx = (k-1)*3 + 1;
        yidx = (k-1)*3 + 2;
        zidx = (k-1)*3 + 3;
        
        plot3(local_XYZ_series(xidx,:),...
              local_XYZ_series(yidx,:),...
              local_XYZ_series(zidx,:))
        view(60,30)
        
        % indicate initial state with a blob
        init_coords = local_XYZ_series([xidx,yidx,zidx],1);
        plot3(init_coords(1),init_coords(2),init_coords(3),'r*')
        
        % plot translucent sphere
        [x,y,z] = sphere(128);
        h = surfl(f*x, f*y, f*z); 
        set(h, 'FaceAlpha', 0.1, 'FaceColor', 'k', 'edgecolor', 'none')
        shading interp
        xlabel(axis_labels{1})
        ylabel(axis_labels{2})
        zlabel(axis_labels{3})
        if N~=1, title(sprintf('n=%i',k)), end
    end   
end

